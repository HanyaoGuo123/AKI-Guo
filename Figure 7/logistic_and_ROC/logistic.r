library(pROC)
data <- read.table("D:\\肾损伤数据\\人的数据\\GSE30718\\data_for_ROC.txt",header = TRUE,row.names = 1)

#=============================================Univariate=====================================
#Calculate logistic regression model
fit1<-glm(OS~CBFB,data=data,family=binomial())
fit2<-glm(OS~COL1A1,data=data,family=binomial())
fit3<-glm(OS~EGF,data=data,family=binomial())
summary(fit)

#CBFB
fitted.prob<-predict(fit1, newdata = data, type = "response")
data$pred1<-fit1$fitted.values
roc_multivar_1<-roc(data$OS,data[,"pred1"],smooth = TRUE)
#COL1A1
fitted.prob<-predict(fit2, newdata = data, type = "response")
data$pred2<-fit2$fitted.values
roc_multivar_2<-roc(data$OS,data[,"pred2"],smooth = TRUE)
#EGF
fitted.prob<-predict(fit3, newdata = data, type = "response")
data$pred3<-fit3$fitted.values
roc_multivar_3<-roc(data$OS,data[,"pred3"],smooth = TRUE)

pdf("GSE30718_ROC.pdf",width=6,height=6)
plot.roc(roc_multivar_1,col="#6475e3fe",print.auc=TRUE) 
plot.roc(roc_multivar_2,col="#e36464fe",print.auc=TRUE,add = TRUE) 
plot.roc(roc_multivar_3,col="#64e386fe",print.auc=TRUE,add = TRUE) 
dev.off()

#=======================================Multivariate (pairwise combination)=========================================
#Calculate logistic regression model
fit1<-glm(OS~CBFB+COL1A1,data=data,family=binomial())
fit2<-glm(OS~COL1A1+EGF,data=data,family=binomial())
fit3<-glm(OS~CBFB+EGF,data=data,family=binomial())
summary(fit)

#CBFB + COL1A1：
fitted.prob<-predict(fit1, newdata = data, type = "response")
data$pred1<-fit1$fitted.values
roc_multivar_1<-roc(data$OS,data[,"pred1"],smooth = TRUE)
#COL1A1 + EGF：
fitted.prob<-predict(fit2, newdata = data, type = "response")
data$pred2<-fit2$fitted.values
roc_multivar_2<-roc(data$OS,data[,"pred2"],smooth = TRUE)
#CBFB + EGF：
fitted.prob<-predict(fit3, newdata = data, type = "response")
data$pred3<-fit3$fitted.values
roc_multivar_3<-roc(data$OS,data[,"pred3"],smooth = TRUE)

pdf("GSE30718_ROC_22.pdf",width=6,height=6)
plot.roc(roc_multivar_1,col="#6475e3fe",print.auc=TRUE) 
plot.roc(roc_multivar_2,col="#e36464fe",print.auc=TRUE,add = TRUE) 
plot.roc(roc_multivar_3,col="#64e386fe",print.auc=TRUE,add = TRUE) 
dev.off()


#=======================Multivariate (all)================================
fit<-glm(OS~CBFB+COL1A1+EGF,data=data,family=binomial())
fitted.prob<-predict(fit, newdata = data, type = "response")
data$pred<-fit$fitted.values
roc_multivar_1<-roc(data$OS,data[,"pred"],smooth = TRUE)
pdf("GSE30718_ROC_all.pdf",width=6,height=6)
plot.roc(roc_multivar_1,col="#000000fe",print.auc=TRUE) 
dev.off()


#Convert to forest plot data
formatFit<-function(fit){
  p<-summary(fit)$coefficients[,4]
  wald<-summary(fit)$coefficients[,3]^2
  valueB<-coef(fit)
  valueOR<-exp(coef(fit))
  confitOR<-exp(confint(fit))
  data.frame(
    B=round(valueB,3),
    Wald=round(wald,3),
    OR_with_CI=paste(round(valueOR,3),"(",
               round(confitOR[,1],3),"~",round(confitOR[,2],3),")",sep=""),
    P=format.pval(p,digits = 3,eps=0.001)
  )
}

result = formatFit(fit3)
result
write.csv(result,'D:\\AKI文章\\ROC曲线\\logistic回归分析\\GSE139061_logistics_result.csv')

