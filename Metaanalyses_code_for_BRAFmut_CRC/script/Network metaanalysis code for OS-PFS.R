setwd('/Users/Desktop/BRAF-CRC/data/second-PFS')
data=read.table('pfs.txt',sep='\t',header=T)
treatments=read.table('treatment.txt',sep='\t',header=T)
study=read.table('study.txt',sep='\t',header=T)
library(gemtc)
library(ggplot2)
library(igraph)
net_gemtc <- mtc.network(data.re = data, treatments = treatments) #\item{data.re}{Relative effect data}
net_gemtc

summary(net_gemtc)

#生成igraph输入
treatments <- data.frame(id = net_gemtc$treatment$id,description = net_gemtc$treatments$description) 
studies <- data.frame(from = summary(net_gemtc)$`Studies per treatment comparison`$t1, 
                      to = summary(net_gemtc)$`Studies per treatment comparison`$t2, 
                      weight = summary(net_gemtc)$`Studies per treatment comparison`$nr)
#生成igraph对象
plot_igraph <- graph_from_data_frame(d = studies, vertices = treatments, directed = F)

#计算每种方案总人数
size_list <-  lapply(split(net_gemtc$data.re$sampleSize, net_gemtc$data.re$treatment), sum) 
treatment_size <- vector(mode = "numeric",length = 0) 
for(i in size_list){ 
  treatment_size = append(treatment_size,i) 
}

#输出图像设置
#l = layout_as_star
l = layout_as_star (plot_igraph, center = V(plot_igraph)[3], order = c(3, 1, 5, 2, 8, 4, 7, 6))
#l = layout_in_circle
#l = layout_in_circle (plot_igraph, order = c(1, 2, 7, 3, 4, 5, 6, 8))

#网状证据图设置
pdf("network plot.pdf") 
plot(plot_igraph,  
     layout = l, #设置布局,使用了前面定义的l，删除这行会使用默认格式 
     vertex.size = treatment_size^(3/5), #圈大小根据样本量进行调整
     edge.width = E(plot_igraph)$weight, #线粗细
     vertex.label = paste(V(plot_igraph)$description,"\nn =",treatment_size), #药物标签名 
     edge.label = E(plot_igraph)$weight, #研究数量标签
     vertex.color="#870087",edge.color="#808080",vertex.label.color="#000000",edge.label.color="#5F005F", #颜色，分别为节点圆圈，节点间连接，药物标签颜色，数字颜色，RGB颜色查询
     vertex.label.cex = 0.8, edge.label.cex = 0.8)
dev.off()

#网状森林图
model <- mtc.model(net_gemtc, type="consistency", n.chain=4, likelihood="binom", link="cloglog", linearModel='random') #一致性模型, HR
#model <- mtc.model(net_gemtc, type="ume", n.chain=4, likelihood="binom", link="cloglog", linearModel='random') #不一致性模型, HR
#cat(model$code) 显示模型建立的code
results <- mtc.run(model, sampler="JAGS", n.adapt=5000, n.iter=20000, thin=1) 
summary(results)

pdf("forest plot.pdf", width = 8, height = 5)
for(i in net_gemtc$treatments$id){ 
  forest(relative.effect(results, i), use.description = TRUE, xlim=c(log(0.01),log(100)), digits=3)
} 
dev.off()



#联赛表
tb<- relative.effect.table(results)
print(tb)
tb1<-round(exp(tb),2)
tb1
write.csv(tb1,"relative_effect_table_hr.csv")

#收敛图
gelman.diag(results)
pdf("gelman plot.pdf")
gelman.plot(results)
dev.off()
#3. 轨迹 Trace 和密度图 Density
pdf("plot.pdf")
plot(results)
dev.off()
#排序表     
#ranks <- rank.probability(results,preferredDirection = 1) #if higher values are preferred
ranks <- rank.probability(results,preferredDirection = -1) #if lower values are preferred
print(ranks) 
write.csv(ranks,"rank.csv")
#单个排序等级图
pdf("ranko.pdf") 
plot(ranks,ylim=c(0,1), beside=TRUE) 
dev.off()
#折线图
rankmerge <- read.csv("rank_merge.csv", header = T, dec = ".")
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
#累积概率图 (堆积排序图)
cumrank.prob <- apply(t(ranks), 2, cumsum)
sucra <- round(colMeans(cumrank.prob[-nrow(cumrank.prob),]),4)
print(sucra)
write.csv(sucra, "SUCRA.csv")
# or
sucra(ranks)
# 异质性检验
result.anohe <- mtc.anohe(net_gemtc, likelihood="binom", link="cloglog", linearModel='random')
summary.anohe <- summary(result.anohe)

pdf("anohe.pdf", width=8, height=15)
plot(summary.anohe, xlim=log(c(0.1, 50)))
a<-summary(result.anohe)
print(a)
dev.off()
#不一致检验
result.node <- mtc.nodesplit(net_gemtc, likelihood="binom", link="cloglog", linearModel='random')

pdf("node-splitting.pdf") 
plot(summary(result.node))
b<-summary(result.node)
print(b)
dev.off()



#meta回归分析

# 清除环境并加载必要的库
# rm(list = ls())
library(gemtc)
library(rjags)

# 读取数据
data <- read.table('OS.txt', sep='\t', header=TRUE)
study <- read.table('study.txt', sep='\t', header=TRUE)

# 创建网络
net_gemtc <- mtc.network(data.re = data,studies = study)

# 标准随机效应模型
# model_standard <- mtc.model(net_gemtc, 
#                             linearModel = "random",
#                             n.chain = 4)
model_standard <- mtc.model(net_gemtc, 
                            likelihood = "normal",
                            link = "identity",
                            linearModel = "random",
                            n.chain = 4)

# 运行标准模型
results_standard <- mtc.run(model_standard, 
                            n.adapt = 5000, 
                            n.iter = 50000, 
                            thin = 10)

# 1. (size)的Meta回归
# 准备研究级数据
studies_data <- data.frame(
  study = as.character(unique(study$study)),
  size = study$size,
  design = as.numeric(study$study.design == "RCT"),
  source = as.numeric(study$data.source == "SUB"),
  type = as.numeric(study$mut.type == "Nonmixed")
)

# 标准化size变量
studies_data$size_std <- scale(studies_data$size)[,1]

# 创建包含协变量的网络
net_size <- mtc.network(data.re = data, studies = studies_data)

# 关键修改：为回归器添加control参数
# 为规模变量定义回归器，使用治疗1作为控制组
regressor_size <- list(
  coefficient = "shared",  # 共享系数
  variable = "size_std",   # 使用标准化的规模变量
  control = "1"            # 指定对照组 - 重要！
)

# 创建规模回归模型
model_size <- mtc.model(net_size, 
                        linearModel = "random", 
                        likelihood = "normal",
                        link = "identity",
                        type = "regression",
                        regressor = regressor_size,
                        n.chain = 4)

# 运行规模回归
results_size <- mtc.run(model_size, 
                        n.adapt = 5000, 
                        n.iter = 50000, 
                        thin = 10)

print(summary(results_size))
# 2. (design)的Meta回归
regressor_design <- list(
  coefficient = "shared",
  variable = "design",
  control = "1"            # 指定对照组
)

model_design <- mtc.model(net_size, 
                          linearModel = "random", 
                          likelihood = "normal",
                          link = "identity",
                          type = "regression",
                          regressor = regressor_design,
                          n.chain = 4)

results_design <- mtc.run(model_design, 
                          n.adapt = 5000, 
                          n.iter = 50000, 
                          thin = 10)

print(summary(results_design))

# 4. (source)的Meta回归
regressor_source <- list(
  coefficient = "shared",
  variable = "source",
  control = "1"            # 指定对照组
)

model_source <- mtc.model(net_size, 
                          linearModel = "random", 
                          likelihood = "normal",
                          link = "identity",
                          type = "regression",
                          regressor = regressor_source,
                          n.chain = 4)

results_source <- mtc.run(model_source, 
                          n.adapt = 5000, 
                          n.iter = 50000, 
                          thin = 10)
print(summary(results_source))

# 4. (type)的Meta回归
regressor_type <- list(
  coefficient = "shared",
  variable = "type",
  control = "1"            # 指定对照组
)

model_type <- mtc.model(net_size, 
                        linearModel = "random", 
                        likelihood = "normal",
                        link = "identity",
                        type = "regression",
                        regressor = regressor_type,
                        n.chain = 4)

results_type <- mtc.run(model_type, 
                        n.adapt = 5000, 
                        n.iter = 50000, 
                        thin = 10)
print(summary(results_type))     
