library(dplyr)

#计算gene-gene相似

#读入非人类 gene~phene
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
        
        if(TRUE){
          #pcc计算
          sim = cor(a,b,method = "pearson")
          #可能协方差为0
          sim = ifelse( is.na(sim), 1/2, (sim +1)/2 )
        }
        if(FALSE){
        #vss计算
        #sim = (cor(a,b,method = "pearson") + 1) / 2#sum(a * b)
        #sim = sim / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
        }
        
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
  nonFull_GeneGene_pcc = full_GeneGene
  nonFull_GeneGene_pcc[which(nonFull_GeneGene_pcc>0)] = 0
  
  sapply(unique(nonHuman_GeneGene$gene1), function(x){
    nHGG = nonHuman_GeneGene[which(nonHuman_GeneGene$gene1==x),]
    nonFull_GeneGene_pcc[x, nHGG$gene2] <<- nHGG$value
    ""
  })
  
  for(i in 1:nrow(nonFull_GeneGene_pcc))
    nonFull_GeneGene_pcc[i,i] = 0
  

save(nonFull_GeneGene_pcc, file = "nonFull_GeneGene_pcc.RData")




