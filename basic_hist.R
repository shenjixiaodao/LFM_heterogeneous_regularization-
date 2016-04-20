library(dplyr)

#1，为了节约计算量，只考虑已确定关系的 gene 和 phene
#2，已确定关系的 gene 和 phene中， gene～gene不存在相似值将被去掉

#Human Gene~Phene
Human_GenePhene = tbl_df(read.csv("GenePhene_Human.csv", header = FALSE, sep = ","))
Human_GenePhene = Human_GenePhene %>% mutate(V1 = as.character(V1), V2 = as.character(V2))
colnames(Human_GenePhene) = c("gene","phene","value")
#gene~gene
GeneGene = tbl_df(read.csv("GeneGene.csv",header = FALSE, sep = ","))
GeneGene = GeneGene %>% mutate(V1 = as.character(V1), V2 = as.character(V2))
GeneGene = GeneGene %>% filter(V1 != V2)
colnames(GeneGene) = c("gene1","gene2","value")
#Phene~Phene
PhenePhene = tbl_df(read.csv("PhenePhene.csv", header = FALSE, sep = ","))
PhenePhene = PhenePhene %>% mutate(V1 = as.character(V1), V2 = as.character(V2))
PhenePhene = PhenePhene %>% filter(V1 != V2)
colnames(PhenePhene) = c("phene1","phene2","value")
# 有些基因在 Human_GenePhene中存在， 但在GeneGene中没有
  #只需要保留 已知的人类基因和表现型
for(i in 1:4){
  genes = intersect(Human_GenePhene$gene, union(GeneGene$gene1, GeneGene$gene2))
  phenes = intersect(Human_GenePhene$phene, union(PhenePhene$phene1, PhenePhene$phene2))
  
  #只需要保留 GeneGene中的基因, 和 PhenePhene中的表现型
  Human_GenePhene = Human_GenePhene %>% filter(gene %in% genes, phene %in% phenes)
    #过滤掉未知 gene 和 phene
  PhenePhene = PhenePhene %>% filter(phene1 %in% phenes, phene2 %in% phenes)
  GeneGene = GeneGene %>% filter(gene1 %in% genes, gene2 %in% genes)
}


  #生成全矩阵
full_Human_GenePene = xtabs(value ~ gene + phene, data = Human_GenePhene)

full_GeneGene = xtabs(value ~ gene1 + gene2, data = GeneGene)


full_PhenePhene = xtabs(value ~ phene1 + phene2, data = PhenePhene)



#读入非人类 gene~phene
if(FALSE){
NonHuman_GenePhene = list()
for(i in 1:8){
  filename = paste0("GenePhene_S",i,".csv")
  NonHuman_GenePhene[[i]] = tbl_df(read.csv(filename, header = FALSE, sep = ",", colClasses = "character"))
  colnames(NonHuman_GenePhene[[i]]) = c("gene","phene","value")
  NonHuman_GenePhene[[i]] = NonHuman_GenePhene[[i]] %>% filter(gene %in% genes)
}
NonFullHuman_GenePhene = sapply(NonHuman_GenePhene, function(x){
  x$value = as.numeric(x$value)
  xtabs(value ~ gene + phene, data = x)
})
#将NonFullHuman_GenePhene 矩阵按列合并, 列编号递增
lens = sapply(NonFullHuman_GenePhene, function(x){
  ncol(x)
})
nonFullHuman_GenePhene = matrix(rep(0, length(genes) * sum(lens)), nrow = length(genes))
rownames(nonFullHuman_GenePhene) = sort(genes)

for(gene in rownames(nonFullHuman_GenePhene)){
  sim = c()
  for(i in 1:length(NonFullHuman_GenePhene)){
    index = which(rownames(NonFullHuman_GenePhene[[i]])==gene)
    if(length(index) == 0)
      sim = c(sim, rep(0, ncol(NonFullHuman_GenePhene[[i]])))
    else
      sim = c(sim,NonFullHuman_GenePhene[[i]][index,])
  }
  nonFullHuman_GenePhene[gene,] = sim
}

#计算 gene 相似值
gene1 = c()
gene2 = c()
value = c()
len_index = c()#用于计算相似值的向量维度
for(i in 1:(nrow(nonFullHuman_GenePhene)-1)){
  for(j in (i+1):nrow(nonFullHuman_GenePhene)){
    index = which(nonFullHuman_GenePhene[i,]>0 & nonFullHuman_GenePhene[j,]>0)
    #判断需要多少个 index
    #增加一项 缩放因子 factor=length(index)/max(index), 代表index个数越多越可信
    if(length(index) > 1){
      len_index = c(len_index, length(index))
      a = nonFullHuman_GenePhene[i,index]
      b = nonFullHuman_GenePhene[j,index]
      sim = sum(a * b)
      sim = sim / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
      
      gene1 = c(gene1,i)
      gene2 = c(gene2,j)
      value = c(value,sim)
    }
  }
}

nonHuman_GeneGene = data.frame(gene1=gene1, gene2=gene2, value=value, len_index=len_index)
#factor 设为90%的分位点
nonHuman_GeneGene$factor = nonHuman_GeneGene$len_index/quantile(nonHuman_GeneGene$len_index,
                                        probs = c(seq(0,1,by = 0.1)))[10]
temp = sapply(nonHuman_GeneGene$factor, function(x){
  if(x>1)
    1
  else
    x
})
nonHuman_GeneGene$factor = temp
nonHuman_GeneGene$value = nonHuman_GeneGene$value * nonHuman_GeneGene$factor

#生成使用其它物种的gene和phene之间生成的gene相似矩阵
nonFull_GeneGene = full_GeneGene
nonFull_GeneGene[which(nonFull_GeneGene>0)] = 0

sapply(unique(nonHuman_GeneGene$gene1), function(x){
  nHGG = nonHuman_GeneGene[which(nonHuman_GeneGene$gene1==x),]
  nonFull_GeneGene[x, nHGG$gene2] <<- nHGG$value
  ""
})

for(i in 1:nrow(nonFull_GeneGene))
  nonFull_GeneGene[i,i] = 1

}

# 1/3 作为测试集， 2/3 作为训练集
# 选取测试集时，不可让某个元素的 行和列 全为0
set.seed(123456789)
temp_fHGP = full_Human_GenePene
test_gene=c()
test_phene=c()
train_gene=c()
train_phene=c()
for(r in 1:nrow(temp_fHGP)){
  for(c in which(temp_fHGP[r,] == 1)){
    # 行和列 不全为0 的正样本作为测试对象, 90%的概率置为预测样本
    temp_row = sum(temp_fHGP[r,])
    temp_col = sum(temp_fHGP[,c])
    temp_sum = temp_row + temp_col 
    if(FALSE){
      if(temp_sum == 3 & sample(c(1,0),size = 1, prob = c(0.4,0.6))){
        test_gene = c(test_gene,r)
        test_phene = c(test_phene,c)
        temp_fHGP[r, c] = 0
      }else if(temp_sum > 3 & temp_sum <= 5 & sample(c(1,0),size = 1, prob = c(0.7,0.3)) == 1){
        test_gene = c(test_gene,r)
        test_phene = c(test_phene,c)
        temp_fHGP[r, c] = 0
      }else if(temp_sum > 5 & sample(c(1,0),size = 1, prob = c(0.8,0.2)) == 1){
        test_gene = c(test_gene,r)
        test_phene = c(test_phene,c)
        temp_fHGP[r, c] = 0
      }else{
        train_gene = c(train_gene,r)
        train_phene = c(train_phene,c)
      }
    }
    if(sample(c(1,0),size = 1, prob = c(0.4,0.6)) & length(which(test_gene==r))==0 &length(which(test_phene==c))==0){
      test_gene = c(test_gene,r)
      test_phene = c(test_phene,c)
      
    }else{
      train_gene = c(train_gene,r)
      train_phene = c(train_phene,c)
    }
    
    
  }
}
test_set = data.frame(gene=test_gene, phene=test_phene, value = 0)
train_set = data.frame(gene=train_gene, phene=train_phene, value = 1)
  #将 full_Human_GenePhene 的测试对象置为 0
for(i in 1:nrow(test_set)){
  ts = test_set[i,]
  full_Human_GenePene[ts[,1],ts[,2]] = 0
}

if(FALSE){
  #统计 行和列 中1的个数
sum_row = sapply(1:nrow(full_Human_GenePene), function(x){
  sum(full_Human_GenePene[x,])
})
sum_col = sapply(1:ncol(full_Human_GenePene), function(x){
  sum(full_Human_GenePene[,x])
})
}


#1， 预测pathMatrix值小的 gene～phene对
#选取 负  样本
#将full_GeneGene归一化
nor_full_GeneGene = (full_GeneGene - min(full_GeneGene[which(full_GeneGene>0)])) / (max(full_GeneGene) - min(full_GeneGene[which(full_GeneGene>0)]))
#将对角线置为  1
for(i in 1:nrow(full_GeneGene))
  full_GeneGene[i,i] = 1  
#将对角线置为  1
for(i in 1:nrow(full_PhenePhene))
  full_PhenePhene[i,i] = 1
#PathMatrix = nonFull_GeneGene %*% full_Human_GenePene %*% full_PhenePhene
PathMatrix = full_GeneGene %*% full_Human_GenePene %*% full_PhenePhene
  #将值为0的设为 负 样本
#将短格式转成长格式  ,速度更快
require(reshape2)
PathMatrix_Frame = as.data.frame(PathMatrix)
colnames(PathMatrix_Frame) = c(1:ncol(PathMatrix_Frame))
PathMatrix_Frame$gene = c(1:nrow(PathMatrix_Frame))#rownames(PathMatrix_Frame)
PathMatrix_Frame = tbl_df(melt(PathMatrix_Frame, "gene"))
colnames(PathMatrix_Frame) = c("gene","phene","value")
PathMatrix_Frame$phene = as.integer(PathMatrix_Frame$phene)
#（负样本在选取时， 同时满足该样本元素的 行和列 不全为0？）   
NegativeSample = PathMatrix_Frame %>% filter(value < 0.77)
#非零数据 且有序排列
NonZeroSample = PathMatrix_Frame %>% filter(value > 0.000000005) %>% arrange(!desc(value))
#预测 有信息支撑的gene～phene对
PredictSample = PathMatrix_Frame %>% filter(value > 0.4)

if(FALSE){
  
  temp = PathMatrix_Frame %>% filter(value > 0.9)
  
  index = sapply(1:nrow(train_set), function(x){
    hgp = train_set[x,]
    which(temp$gene==hgp$gene & temp$phene==hgp$phene)
  })
  SpecialSample = temp[-index,]
}

hgp_values = sapply(1:nrow(train_set), function(x){
  hgp = train_set[x,]
  PathMatrix_Frame[which(PathMatrix_Frame$gene==hgp$gene & PathMatrix_Frame$phene==hgp$phene),]$value
})

test_hgp_values = sapply(1:nrow(test_set), function(x){
  hgp = test_set[x,]
  PathMatrix_Frame[which(PathMatrix_Frame$gene==hgp$gene & PathMatrix_Frame$phene==hgp$phene),]$value
})



#2， 预测hetesim相似值小的 gene～phene对

for(i in 1:nrow(test_set)){
  ts = as.data.frame(test_set[i,])
  test_set[i,]$sim = PathMatrix[ts[,"gene"], ts[,"phene"]]
}

temp = which(full_Human_GenePene==1 & PathMatrix==1)

rows = rownames(PathMatrix)[temp%%nrow(PathMatrix)]
cols = colnames(PathMatrix)[as.integer(temp/nrow(PathMatrix)) + 1]

GP_H = full_GeneGene %*% full_Human_GenePene
H_GP = full_Human_GenePene %*% full_PhenePhene
temp = (GP_H + H_GP) / 2
#gene="1502", phene="1586"是异常值， 强令为1
temp["1502","1586"] = 1

PathMatrix_Frame = as.data.frame(temp)
colnames(PathMatrix_Frame) = c(1:ncol(PathMatrix_Frame))
PathMatrix_Frame$gene = c(1:nrow(PathMatrix_Frame))#rownames(PathMatrix_Frame)
PathMatrix_Frame = tbl_df(melt(PathMatrix_Frame, "gene"))
colnames(PathMatrix_Frame) = c("gene","phene","value")
PathMatrix_Frame$phene = as.integer(PathMatrix_Frame$phene)
hgp_values = sapply(1:nrow(train_set), function(x){
  hgp = train_set[x,]
  print(x)
  PathMatrix_Frame[which(PathMatrix_Frame$gene==hgp$gene & PathMatrix_Frame$phene==hgp$phene),]$value
})

test_hgp_values = sapply(1:nrow(test_set), function(x){
  hgp = test_set[x,]
  PathMatrix_Frame[which(PathMatrix_Frame$gene==hgp$gene & PathMatrix_Frame$phene==hgp$phene),]$value
})

nozero_test_index = unlist(sapply(1:nrow(test_set), function(x){
  # > 1
  if(sum(full_Human_GenePene[,test_set[x,]$phene])>0)
    x
}))

single_test_index = unlist(sapply(1:nrow(test_set), function(x){
  # > 1
  if(sum(full_Human_GenePene[,test_set[x,]$phene])==1)
    x
}))


tmp = test_set[nozero_test_index,]
tmp_PMF = PathMatrix_Frame[order(PathMatrix_Frame$value, decreasing = TRUE),]
sapply(1:100, function(i){
  count = 0
  for(j in 1:nrow(tmp)){
    ts = tmp[j,]
    #PathMatrix_Frame[which(PathMatrix_Frame$gene==ts$gene),]
    if(length(which(tmp_PMF[which(tmp_PMF$gene==ts$gene),]$phene[1:i]==ts$phene))>0)
      count = count + 1
  }
  #print(i)
  count
})
