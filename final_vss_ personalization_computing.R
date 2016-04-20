#归一化 HumanNet 的gene-gene相似矩阵
nor_full_GeneGene = (full_GeneGene - min(full_GeneGene[which(full_GeneGene>0)])) / (max(full_GeneGene) - min(full_GeneGene[which(full_GeneGene>0)]))

set.seed(123456789)
#nnumber : 负样本比率
# NegativeSample : 计算得到的大概率负样本实例子， SpecialSample ：计算得到的小概率负样本实例子
# train_set ：数据中已知的正样本训练数据
TS = as.data.frame(rbind(train_set,NegativeSample[sample(1:nrow(NegativeSample),
          nrow(train_set) * nnumber),], SpecialSample[sample(1:nrow(SpecialSample), nrow(train_set) * 0.1),]))
maxIter = 200 #最大迭代次数
theta = 0.001 #步长
lambda = 0.001 #
gthreshold = 0.00001#
pthreshold = 0.00001#
for(nnumber in c(2)){
  
  #将训练数据TS每条记录装入list
  nozero_inv = list()
  for(i in 1:nrow(TS)){
    nozero_inv[[length(nozero_inv)+1]] = list(gene = TS[i,1],
                        phene = TS[i,2], value = TS[i,3])
  }

  for(ratio in c(1.5)){
    for(fs in c(10,20,30,40,50,60,70,80)){
      for(be in c(0.001)){
        for(al in c(0.5)){
          #--------------------------------
          #-----------------------------------------------------------------------
          alpha = al#0.03 # 基因相似约束因子
          alpha2 = ratio * al
          beta = be # 表现型相似约束因子
          features = fs
          gene_nbs = list()#基因相似矩阵1邻居缓存
          gene_nbs2 = list()#基因相似矩阵2邻居缓存
          phene_nbs = list()#表型相似矩阵邻居缓存
          
          #存放模型数据文件
          tempfile = paste0("temp_data/features/temp_alpha",al,"_beta",be,"_features",features,"_number",nnumber,".RData")

          #----------------------- 随机初始特征因子矩阵 ----------------------------
          gene_features_dist_reg = matrix(runif(nrow(srcMatrix) * features, min = 0, max = 1), 
                                          nrow = nrow(srcMatrix) , ncol = features)
          phene_features_dist_reg = matrix(runif(ncol(srcMatrix) * features , min = 0, max = 1),
                                           nrow = ncol(srcMatrix) , ncol = features)
          
          #----------------------------------
          
          for(j in 1:features){
            tmp_gradient = 10
            for(iter in 1:maxIter){
              #if(abs(rmse_last - rmse) < min_improvement){
              #  break
              #}
              sq = 0
              for(i in 1:length(nozero_inv)){
                inv = nozero_inv[[i]]
                
                psum = 1 - 1/(1 - exp(-sum((gene_features_dist_reg[inv$gene,] - phene_features_dist_reg[inv$phene,])^2)/2 + lambda))
                error = inv$value + (1 - inv$value) * psum
                
                
                temp1 = gene_features_dist_reg[inv$gene,j]#暂存，用于更新表型特征
                temp2 = phene_features_dist_reg[inv$phene,j]#暂存，用于更新基因特征
                
                #基因特征, 若已加入内存，则读取内存数据
                if(is.null(gene_nbs[inv$gene][[1]])){
                  # 相似度非0
                  order = which(nonFull_GeneGene[inv$gene,]>gthreshold)
                  if(length(order) >=min_nbs){
                    sim = nonFull_GeneGene[inv$gene,order]
                    sim = sim #/ sum(sim)#相似度 归一化
                    nbs = data.frame(order=order, sim=sim)
                  }else
                    nbs = data.frame(order=c(), sim=c())
                  gene_nbs[[inv$gene]] = nbs
                }else
                  nbs = gene_nbs[[inv$gene]]
                gradient = error * (temp2 - temp1)
                if(!is.null(nbs$order) & length(nbs$order)>0)
                  gradient = gradient - alpha * sum(nbs$sim * (temp1 - gene_features_dist_reg[nbs$order,j])) #约束项
                
                if(TRUE){
                  #基因特征2, 若已加入内存，则读取内存数据
                  if(is.null(gene_nbs2[inv$gene][[1]])){
                    # 相似度非0 或 大于指定值
                    order = which(nor_full_GeneGene[inv$gene,]>gthreshold)
                    if(length(order) >=min_nbs){
                      sim = nor_full_GeneGene[inv$gene,order]
                      sim = sim #/ sum(sim)#相似度 归一化
                      nbs = data.frame(order=order, sim=sim)
                    }else
                      nbs = data.frame(order=c(), sim=c())
                    gene_nbs2[[inv$gene]] = nbs
                  }else
                    nbs = gene_nbs2[[inv$gene]]
                  if(!is.null(nbs$order) & length(nbs$order)>0)
                    gradient = gradient - alpha2 * sum(nbs$sim * (temp1 - gene_features_dist_reg[nbs$order,j])) #约束项
                }
                #基因特征更新
                gene_features_dist_reg[inv$gene,j] = temp1 + theta * gradient
                #
                #表现型特征
                if(is.null(phene_nbs[inv$phene][[1]])){
                  # 相似度非0 或 大于指定值
                  order = which(full_PhenePhene[inv$phene,]>pthreshold)
                  if(length(order)>=min_nbs){
                    sim = full_PhenePhene[inv$phene,order]
                    sim = sim #/ sum(sim)#相似度 归一化
                    nbs = data.frame(order=order, sim=sim)
                  }else
                    nbs = data.frame(order=c(), sim=c())
                  phene_nbs[[inv$phene]] = nbs
                }else
                  nbs = phene_nbs[[inv$phene]]
                gradient = error * (temp1 - temp2)
                if(!is.null(nbs$order) & length(nbs$order)>0)
                  gradient = gradient - beta * sum(nbs$sim * (temp2 - phene_features_dist_reg[nbs$order,j]))
                phene_features_dist_reg[inv$phene,j] = temp2 + theta * gradient
              }
              
            }
          }
          
          #计算模型预测值
          temp = sapply(1:nrow(gene_features_dist_reg), function(x){
            sapply(1:nrow(phene_features_dist_reg), function(y){
              1/exp(sum((gene_features_dist_reg[x,] - phene_features_dist_reg[y,])^2)/2 + lambda)
            })
          })
          save(temp,file = tempfile)#保存模型结果

        }
      }
    }
  }


}





