setwd('/Users/Desktop/BRAF-CRC/data/second-ORR')
data=read.table('pairwiseORR.txt',sep='\t',header=T)
#ORR-Pairwise分析
library(meta)
library(dplyr)

corrected_data <- data %>%
  mutate(
    across(c(event1, event2), ~ if_else(.x == 0 | .x == total1, .x + 0.5, .x)),
    across(c(total1, total2), ~ if_else(event1 == 0 | event1 == total1, .x + 1, .x))
  )

# Meta分析
meta_res <- metabin(event.e = event1, n.e = total1,
                    event.c = event2, n.c = total2,
                    studlab = paste("Study", study),
                    data = corrected_data,
                    subgroup = subgroup,
                    sm = "OR",
                    method = "MH",
                    common = TRUE,
                    random = TRUE,
                    method.tau = "REML",
                    hakn = TRUE,
                    title = "Subgroup Analysis of Odds Ratios")

# BMJ格式森林图
forest(meta_res,
       layout = "BMJ",
       # 颜色设置
       col.diamond.fixed = "#2171B5",    # BMJ标准蓝色
       col.diamond.random = "#CB181D",   # BMJ标准红色
       col.square = "#238B45",           # 森林绿
       col.square.lines = "#238B45",
       col.study = "black",
       col.inside = "white",             
       # 坐标轴设置
       xlim = c(0.1, 10),
       xlab = "Odds Ratio (log scale)", 
       # 标签设置
       leftlabs = c("Study", "Events", "Total"),
       rightlabs = c("OR", "95% CI"),
       # 字体设置
       fontsize = 10,
       fs.heading = 11,
       ff.heading = "bold",
       # 图形参数
       squaresize = 0.8,
       lty.fixed = 0,                   # 隐藏固定效应线
       lty.random = 0,                  # 隐藏随机效应线
       # 统计显示
       print.I2 = TRUE,
       print.tau2 = TRUE,
       print.p = TRUE,
       print.Q = TRUE,
       print.subgroup.labels = TRUE,
       test.subgroup = TRUE,
       # 布局优化
       colgap.forest.left = "10mm",
       spacing = 1.2,
       calcwidth.subgroup = TRUE)

# 输出高清TIFF
tiff("BMJ_OR_Forest.tiff",
     width = 410, height = 297, units = "mm", 
     res = 600, compression = "lzw")
forest(meta_res)
dev.off()