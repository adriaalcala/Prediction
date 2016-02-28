library(AligNet)
library(igraph)
load("Alignments/global-cel-dme.RData")

cel = AligNet::read.network.edges("Nets/cel.tab", cols = c(1,2),sep="\t")
dme = AligNet::read.network.edges("Nets/dme.tab", cols = c(1,2),sep="\t")
testpair = function(p1,p2,net1,net2,aligns){
p1 = as.character(p1)
p2 = as.character(p2)
belongs.p1 = unlist(lapply(aligns, function(i) p1%in%names(i)))
belongs.p2 = unlist(lapply(aligns, function(i) p2%in%names(i)))
bp1.p2 = which(belongs.p1*belongs.p2==1)[1]
if(is.na(bp1.p2)){
  test = c(-10)
}
else{
  test = c(bp1.p2)
}

bp1.p2 = which(belongs.p1+belongs.p2==1)[1]
if(is.na(bp1.p2)){
  test = c(test,-10)
}
else{
  test = c(test,bp1.p2)
}
common.neighs = length(intersect(names(neighbors(cel,p1)),names(neighbors(cel,p2))))
test  = c(test,common.neighs)
for(align in aligns){
  if(p1 %in% names(align)){
    if (p2 %in% names(align)){
      image.p1 = align[p1]
      image.p2 = align[p2]
      if(are_adjacent(net2,image.p1,image.p2)){
        test2 = c(test2,1)
      }
      else{
        test = c(test, -1)
      }
    }
    else{
      test = c(test,0)
    }
  }
  else{
    test = c(test,0)
  }
  }
test = c(test, as.numeric(are_adjacent(net1, p1,p2)))
return(test)
}


data2 = data.table(expand.grid(V(cel),V(cel)$))
test.list = lapply(1:dim(data2)[1], function(i) testpair(data2[i,]$Var1, data2[i,]$Var2,cel, dme, global[[1]]))
data = rbindlist(lapply(test.list, function(i) as.list(i)))
save(data, file="Data-Cel-Dme.RData")
#test.list = lapply(V(cel)$name, function(i) testpair("ce3956", i, cel, dme, global[[1]]))
#test.list2 =lapply(V(cel)$name, function(i) testpair("ce1014", i, cel, dme, global[[1]]))
#test.list3 =lapply(V(cel)$name, function(i) testpair("ce412", i, cel, dme, global[[1]]))
#test.list4 =lapply(V(cel)$name, function(i) testpair("ce3674", i, cel, dme, global[[1]]))
 
#data1 = rbindlist(lapply(test.list, function(i) as.list(i)))
#data2 = rbindlist(lapply(test.list2, function(i) as.list(i)))
#data3 = rbindlist(lapply(test.list3, function(i) as.list(i)))
#data4 = rbindlist(lapply(test.list4, function(i) as.list(i)))

#data = rbindlist(list(data1,data2,data3,data4))
# Set a seed
#set.seed(500)

# Train-test random splitting for linear model
# index <- sample(1:nrow(data),round(0.75*nrow(data)))
# train <- data[index,]
# test <- data[-index,]
# 
# #-------------------------------------------------------------------------------
# # Neural net fitting
# 
# # NN training
# library(neuralnet)
# n <- names(train)
# f <- as.formula(paste("V25 ~", paste(n[!n %in% "V25"], collapse = " + ")))
# nn <- neuralnet(f,data=train,hidden=c(100,50,25,10,5), lifesign = "full", rep=100)
# 
# # Visual plot of the model
# plot(nn)
# 
# # Predict
# pr.nn <- compute(nn,test[,1:24, with=FALSE])
# 
# # Results from NN are normalized (scaled)
# # Descaling for comparison
# pr.nn_ <- pr.nn$net.result
# test.r <- test$V25
# 
# # Calculating MSE
# MSE.nn <- sum((test.r - pr.nn_)^2)/nrow(test)
# 
# 
# # Plot predictions
# 
# plot(test$V25,pr.nn_,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
# abline(0,1,lwd=2)
# legend('bottomright',legend='NN',pch=18,col='red', bty='n')
# 
# P = length(which(test$V25==1))
# TP = lapply(seq(0,1,0.001), function(i) length(intersect(which(pr.nn_>i), which(test$V25==1))))
# tp = unlist(TP)/P
# N = length(which(test$V25==0))
# TN = lapply(seq(0,1,0.001), function(i) length(intersect(which(pr.nn_<=i), which(test$V25==0))))
# tn = unlist(TN)/N
# FP = lapply(seq(0,1,0.001), function(i) length(intersect(which(pr.nn_>i), which(test$V25==0))))
# fp = unlist(FP)/N
# FN = lapply(seq(0,1,0.001), function(i) length(intersect(which(pr.nn_<=i), which(test$V25==1))))
# fn = unlist(FN)/P
# TP = unlist(TP)
# TN = unlist(TN)
# FP = unlist(FP)
# FN = unlist(FN)
# ACC  = (TP+TN)/(P+N)
# PPV = (TP)/(TP+FP)
# NPV = (TN)/(TN+FN)
# 
# LR1 = tp/fp
# LR2 = fn/tn
# 
# DOR = LR1/LR2
# 
# plot(NPV, type = "l")
# lines(ACC,type="l", col = "green")
# lines(PPV, type="l")
# lines(NPV, type="l")
# 
# plot(tp,type = 'l')
# lines(tn, col = "green")
# lines(fp, col = "blue")
# lines(fn, col = "red")
# #-------------------------------------------------------------------------------
# # Cross validating
# 
# library(boot)
# set.seed(200)
# 
# # Neural net cross validation
# set.seed(450)
# cv.error <- NULL
# k <- 10
# 
# # Initialize progress bar
# library(plyr) 
# pbar <- create_progress_bar('text')
# pbar$init(k)
# 
# for(i in 1:k){
#   index <- sample(1:nrow(data),round(0.9*nrow(data)))
#   train.cv <- scaled[index,]
#   test.cv <- scaled[-index,]
#   
#   nn <- neuralnet(f,data=train.cv,hidden=c(5,2),linear.output=T)
#   
#   pr.nn <- compute(nn,test.cv[,1:13])
#   pr.nn <- pr.nn$net.result*(max(data$V25)-min(data$V25))+min(data$V25)
#   
#   test.cv.r <- (test.cv$V25)*(max(data$V25)-min(data$V25))+min(data$V25)
#   
#   cv.error[i] <- sum((test.cv.r - pr.nn)^2)/nrow(test.cv)
#   
#   pbar$step()
# }
# 
# # Average MSE
# mean(cv.error)
# 
# # MSE vector from CV
# cv.error
# 
# # Visual plot of CV results
# boxplot(cv.error,xlab='MSE CV',col='cyan',
#         border='blue',names='CV error (MSE)',
#         main='CV error (MSE) for NN',horizontal=TRUE)