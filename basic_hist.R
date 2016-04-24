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
    if(TRUE){
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
    if(TRUE){
      if(sample(c(1,0),size = 1, prob = c(0.4,0.6)) & length(which(test_gene==r))==0 &length(which(test_phene==c))==0){
        test_gene = c(test_gene,r)
        test_phene = c(test_phene,c)
        
      }else{
        train_gene = c(train_gene,r)
        train_phene = c(train_phene,c)
      }
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


#1， 预测pathMatrix值小的 gene～phene对
#选取 负  样本
#将full_GeneGene归一化
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

#SpecialSample 不包含已知为1的，随机游走概率较大的样例
temp = PathMatrix_Frame %>% filter(value > 0.9)
index = sapply(1:nrow(train_set), function(x){
  hgp = train_set[x,]
  which(temp$gene==hgp$gene & temp$phene==hgp$phene)
})
SpecialSample = temp[-index,]



if(FALSE){
#多链接实例子
nozero_test_index = unlist(sapply(1:nrow(test_set), function(x){
  # > 1
  if(sum(full_Human_GenePene[,test_set[x,]$phene])>0)
    x
}))
}

