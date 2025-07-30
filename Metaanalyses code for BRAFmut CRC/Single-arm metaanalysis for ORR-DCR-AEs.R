setwd('/Users/Desktop/BRAF-CRC/data/second-ORR')
data=read.table('ORR.txt',sep='\t',header=T)
#单个率的分析

# 添加连续性校正（当event=0或event=total时）
data <- data %>%
  mutate(
    event_adj = ifelse(event == 0 | event == total, event + 0.5, event),
    total_adj = ifelse(event == 0 | event == total, total + 1, total),
    p_adj = event_adj/total_adj
  )

data <- data %>%
  mutate(
    p = event_adj/total_adj,
    logit_p = log(p/(1-p)),
    se_logit = sqrt(1/(event_adj*(1-p)) + 1/((total_adj-event_adj)*p))
  )


# 进行meta分析
meta_res <- metagen(TE = logit_p,
                    seTE = se_logit,
                    studlab = study,
                    data = data,
                    sm = "PLOGIT",  # logit转换后的比例
                    subgroup = subgroup,
                    common = TRUE,
                    random = TRUE,
                    method.tau = "REML",
                    method.random.ci = "HK",  # Hartung-Knapp调整
                    title = "Single Proportion Meta-Analysis")

# 生成BMJ格式森林图
forest(meta_res,
       layout = "BMJ",
       # 核心修改：设置小数位数
       digits = 3,                      # 主效应量显示3位小数
       #digits.se = 3,                   # 标准误小数位数
       #digits.stat = 3,                 # 统计检验值小数位数
       
       col.diamond.fixed = "#045a8d",    # 固定效应菱形
       col.diamond.random = "#bd0026",   # 随机效应菱形
       col.square = "#2c7fb8",           # 研究方块
       col.study = "black",              
       col.inside = "black",             
       xlim = c(0, 1),                  # 比例范围0-1
       xlab = "Proportion (95% CI)",    
       leftlabs = c("Study", "Events", "Total"), 
       rightlabs = c("Proportion", "95% CI"),
       fontsize = 10,                   
       squaresize = 0.8,                
       lty.fixed = 2,                   
       lty.random = 3,                  
       print.I2 = TRUE,                 
       print.tau2 = TRUE,               
       print.p = TRUE,                  
       print.Q = TRUE,                  
       test.subgroup = TRUE,            # 显示亚组差异检验
       test.overall.fixed = TRUE,       
       test.overall.random = TRUE,
       
       # 添加比例转换
       backtransf = TRUE,               # 反向转换为原始比例
       irscale = 1,                     # 保持比例范围0-1
       irunit = "Proportion",
       
       # 亚组显示优化
       subgroup.name = "Category",      # 亚组标签名称
       print.subgroup.labels = TRUE,    # 显示亚组标签
       col.subgroup = "#666666",        # 亚组标题颜色
       fs.subgroup = 12)                # 亚组字体大小

# 输出高清TIFF
tiff("Single_Proportion_Forest.tiff", 
     width = 210, height = 297, units = "mm", res = 600, compression = "lzw")
forest(meta_res)  # 保持相同参数设置
dev.off()