install.packages("gemtc")  #调用JAGS程序来进行贝叶斯框架下的分析以及相关图形绘制
install.packages("igraph")  #网状证据图绘制
install.packages("ggplot2")  #绘制排序折线图

library(gemtc)
library(igraph)
library(ggplot2)
library(readr)
setwd('/Users/Desktop/BRAF-CRC/data/second-ORR')
data=read.table('NMAORR.txt',sep='\t',header=T)
treatments=read.table('treatment.txt',sep='\t',header=T)
net_gemtc <- mtc.network(data = data, treatments = treatments)
net_gemtc
summary(net_gemtc)

plot(net_gemtc)
#2. 生成igraph对象
treatments <- data.frame(id = net_gemtc$treatment$id,description = net_gemtc$treatments$description) 
studies <- data.frame(from = summary(net_gemtc)$`Studies per treatment comparison`$t1, 
                      to = summary(net_gemtc)$`Studies per treatment comparison`$t2, 
                      weight = summary(net_gemtc)$`Studies per treatment comparison`$nr)
plot_igraph <- graph_from_data_frame(d = studies, vertices = treatments, directed = F)

#计算每种方案总人数
size_list <-  lapply(split(net_gemtc$data.ab$sampleSize, net_gemtc$data.ab$treatment), sum) 
treatment_size <- vector(mode = "numeric",length = 0) 
for(i in size_list){ 
  treatment_size = append(treatment_size,i) 
}
#输出图像设置
#l = layout_as_star
l = layout_as_star (plot_igraph, order = c(3, 5, 7, 4, 2, 6, 1))
#l = layout_in_circle
#l = layout_in_circle (plot_igraph, order = c(1, 2, 7, 3, 4, 5, 6, 8))
l = layout_as_star (plot_igraph, center = V(plot_igraph)[3], order = c(3, 1, 5, 2, 8, 4, 7, 6))

#网状证据图设置
pdf("network plot.pdf") 
plot(plot_igraph,  
     layout = l, #设置布局,使用了前面定义的l，删除这行会使用默认格式 
     vertex.size = treatment_size^(3/5), #圈大小根据样本量进行调整
     edge.width = E(plot_igraph)$weight, #线粗细
     vertex.label = paste(V(plot_igraph)$description,"\nn =",treatment_size), #药物标签名 
     edge.label = E(plot_igraph)$weight, #研究数量标签
     vertex.color = "#B92D2E",edge.color = "grey",vertex.label.color = "black",edge.label.color = "#184A7F", #颜色，分别为节点圆圈，节点间连接，药物标签颜色，数字颜色，RGB颜色查询
     vertex.label.cex = 1.2, edge.label.cex = 1.2) #字号大小
dev.off()

#网状森林图
model <- mtc.model(net_gemtc, type="consistency", n.chain=4, likelihood="binom", link="logit", linearModel='random') #一致性模型, OR 
#model <- mtc.model(net_gemtc, type="ume", n.chain=4, likelihood="binom", link="logit", linearModel='random') #不一致性模型, OR 
#model <- mtc.model(net_gemtc, type="consistency", n.chain=4, likelihood="binom", link="log", linearModel='random') #一致性模型，RR
#cat(model$code) 显示模型建立的code
results <- mtc.run(model, sampler="JAGS", n.adapt=5000, n.iter=20000, thin=1) 
summary(results)

pdf("forest plot.pdf", width = 8, height = 2.5)
for(i in net_gemtc$treatments$id){ 
  # 连续性变量, xlim为横坐标范围，需根据实际情况手动调整，无需对称 
  #forest(relative.effect(results, i), use.description = TRUE, xlim=c(-700,700)) 
  # 二分类变量, xlim为横坐标范围，需根据实际情况手动调整，无需对称 
  forest(relative.effect(results, i), use.description = TRUE, xlim=c(log(0.001),log(1000)), digits=3)
} 
dev.off()

#联赛表
tb<- relative.effect.table(results)
print(tb)
tb1<-round(exp(tb),2)
tb1
write.csv(tb1,"relative_effect_table_or.csv")
#write.csv(tb1,"relative_effect_table_rr.csv")

#收敛性诊断 Convergence Diagnostics 
gelman.diag(results)

#收敛诊断图
pdf("gelman plot.pdf")
gelman.plot(results)
dev.off()

#轨迹 Trace 和密度图 Density
pdf("plot.pdf")
plot(results)

#排序表
ranks <- rank.probability(results,preferredDirection = 1) #if higher values are preferred
#ranks <- rank.probability(results,preferredDirection = -1) #if lower values are preferred
print(ranks) 
write.csv(ranks,"rank.csv")

# 单个排序等级图
pdf("ranko.pdf") 
plot(ranks,ylim=c(0,1), beside=TRUE) 
dev.off()

#折线图
rankmerge <- read.csv("data/binary variable/rank_merge.csv", header = T, dec = ".")
rankmerge

pdf("rank_merge.pdf")
p<-ggplot(rankmerge, aes(x=rank, y=prob, colour=Treatment)) + geom_line(size=1, lty="solid")
p1<-p+geom_point(size=1) + labs(x="Rank", y="Probability") + theme(axis.text = element_text(size=10),axis.title = element_text(size=12),panel.background=element_blank(),axis.line=element_line(color="black"),legend.title=element_text(size=12, color="black"),legend.position="bottom")
plot(p1)
dev.off()

#累积概率图 (堆积排序图)
pdf("rank_cmu.pdf") 
plot(ranks,xlab='treatment',ylab='cumulative probability') 
dev.off()

cumrank.prob <- apply(t(ranks), 2, cumsum)
sucra <- round(colMeans(cumrank.prob[-nrow(cumrank.prob),]),4)
print(sucra)
write.csv(sucra, "SUCRA.csv")

#或者
sucra(ranks)


#异质性检验 Heterogeneity Analysis
result.anohe <- mtc.anohe(net_gemtc, n.adapt=1000, n.iter=5000)
#result.anohe <- mtc.anohe(net_gemtc, n.adapt = 1000, n.iter =5000, thin = 1, n.chain=4,likelihood="binom",link="logit",linearModel="random")
summary.anohe <- summary(result.anohe)

pdf("anohe.pdf", width=8, height=15)
plot(summary.anohe, xlim=log(c(0.001, 1000)))
a<-summary(result.anohe)
print(a)
dev.off()

#不一致性检验 Node-Splitting Plot
result.node <- mtc.nodesplit(net_gemtc, sampler ="JAGS", thin = 50)
#result.node <-mtc.nodesplit(net_gemtc, n.adapt = 20000, n.iter= 50000, thin = 1, n.chain=4,likelihood="binom",link="logit",linearModel="random")

pdf("node-splitting.pdf") 
#plot(summary(result.node))
plot(summary(result.node),xlim=log(c(0.001,15)), digits=3)
b<-summary(result.node)
print(b)
dev.off()

#回归分析 Regression Analysis
# Load the data and create an mtc.network
data=read.table('DCR.txt',sep='\t',header=T)
study=read.table('study.txt',sep='\t',header=T)
networkreg <- mtc.network(data=data, studies=study)
# Random effect meta-regression
model_reg <- mtc.model(networkreg, type="regression", 
                       regressor=list(coefficient='shared',#只有一个系数，多个系数选unrelated
                                      variable='size',
                                      control='1'))
# Run regression model
sink("regression.csv")
result_reg <- mtc.run(model_reg, sampler ="JAGS", n.adapt=1000, n.iter=5000, thin=1)
summary(result_reg)
sink(NULL)
dev.off()

#二分类协变量meta分回归：协变量应变为0/1
model <- mtc.model(networkreg, type="regression",  #默认为OR
                   regressor=list(coefficient='shared',
                                  variable='study_design',
                                  control='1'))
#run regression model
results <- mtc.run(model)
summary(results)

#二分类协变量meta分回归：协变量应变为0/1
model <- mtc.model(networkreg, type="regression",  #默认为OR
                   regressor=list(coefficient='shared',
                                  variable='data_source',
                                  control='1'))
#run regression model
results <- mtc.run(model)
summary(results)

#二分类协变量meta分回归：协变量应变为0/1
model <- mtc.model(networkreg, type="regression",  #默认为OR
                   regressor=list(coefficient='shared',
                                  variable='prior_therapy',
                                  control='1'))
#run regression model
results <- mtc.run(model)
summary(results)
