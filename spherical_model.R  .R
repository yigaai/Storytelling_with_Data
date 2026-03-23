
## Import Libraries: --------------------------------- 


library(dplyr)  
library(ggplot2)
library(afex)
library(tidyverse) 
library(emmeans)  
library(nlme) 
library(lme4)
library(lmerTest) 
library(car) 
library(rstatix)
library(gridExtra)
library(performance)
library(interactions) 
library(mgcv) 
library(plotly) 
library(MixRF) 
library(geometry)  



### Read Data: --------------------------------------------------- 


## Data of Experiment:  
data_exp <- read.csv("C:/Users/GXlab_4/Desktop/data_completed/data_exp.csv",  
                     fileEncoding = "GBK")  

summary(data_exp)       


## Data of Questionnaire: 
scale <- read.csv("C:/Users/GXlab_4/Desktop/问卷数据/questionnaire.csv")   
summary(scale)
sum(is.na(scale))          


## check: 
key_vars <- c("participant_id", "v_row", "v_col", "angle", "distance")  


scale_dupli <- scale %>%  
  group_by(across(all_of(key_vars))) %>%      
  summarise(n = n(), .groups = "drop") %>% 
  filter(n>1)  

print(scale_dupli)  



## delete rows which n >1 :   

scale_clean <- scale%>% 
  anti_join(scale_dupli, by = key_vars)   



## Merge expriemnt data & subjetive data:  

data_total <- left_join(data_exp, scale_clean,     
                        by = c("participant_id", "v_row", "v_col", "angle", "distance")) 
summary(data_total)       



## check:  

nrow(data_exp) 
nrow(data_total) ## rows == 62720       


## Convert the unit of rt from s to ms: 
#dataset$rt_for_calc <- round(dataset$rt_for_calc*1000, 2)     



## Filling the NA of rt: -------------------------------------   


dataset <- data_total %>%       
  group_by(angle, distance, v_row, v_col) %>% 
  mutate(
    group_mean = mean(rt_for_calc[correct == 1], na.rm = TRUE),   
    rt_for_calc = case_when(
      correct == 1 & is.na(rt_for_calc) & !is.nan(group_mean) ~ group_mean,
      TRUE ~ rt_for_calc 
    )
  ) %>%
  select(-group_mean) %>%
  ungroup()  


sum(is.na(dataset$rt_for_calc[dataset$correct==1]))      ## NA； 0            

summary(dataset)     




## Conversion of Target：--------------------------------------------------   


idx <- dataset$v_col %in% 1:5     
dataset[idx, c("target_x_corr", "target_y_corr")] <- dataset[idx, c("target_y", "target_x")]
dataset[idx, "target_x_corr2"] <- dataset[idx, "target_x_corr"] - 0.4  
dataset[idx, "target_y_corr2"] <- dataset[idx, "target_y_corr"] * (-1)  


idx <- dataset$v_col %in% 6:7 
dataset[idx, c("target_x_corr", "target_y_corr")] <- dataset[idx, c("target_y", "target_x")]

dataset[idx, "target_x_corr2"] <- dataset[idx, "target_x_corr"] + 1.6    
dataset[idx, "target_y_corr2"] <- dataset[idx, "target_y_corr"] * (-1)      



## Conversion of Click : ----------------------- 

idx <- dataset$v_col %in% 1:5     
dataset[idx, c("click_x_corr", "click_y_corr")] <- dataset[idx, c("click_y", "click_x")] 
dataset[idx, "click_x_corr2"] <- dataset[idx, "click_x_corr"] - 0.4  
dataset[idx, "click_y_corr2"] <- dataset[idx, "click_y_corr"] *(-1)  


idx <- dataset$v_col %in% 6:7 
dataset[idx, c("click_x_corr", "click_y_corr")] <- dataset[idx, c("click_y", "click_x")]
dataset[idx, "click_x_corr2"] <- dataset[idx, "click_x_corr"] + 1.6    
dataset[idx, "click_y_corr2"] <- dataset[idx, "click_y_corr"] *(-1)     



## Conversion of the edge of grid ---------------------

idx <- dataset$v_col %in% 1:5      
dataset[idx, c("red_top_corr", "red_bottom_corr", "red_right_corr","red_left_corr" )] <-
  dataset[idx, c("red_left", "red_right", "red_top", "red_bottom")]   
dataset[idx, "red_bottom_corr2"] <- dataset[idx, "red_bottom_corr"] * (-0.1)   
dataset[idx, "red_top_corr2"] <- dataset[idx, "red_top_corr"] * (-0.1)   
dataset[idx, "red_left_corr2"] <- dataset[idx, "red_left_corr"] - 0.4 
dataset[idx, "red_right_corr2"] <- dataset[idx, "red_right_corr"] - 0.4   


idx <- dataset$v_col %in% 6:7     
dataset[idx, c("red_top_corr", "red_bottom_corr", "red_right_corr","red_left_corr" )] <- 
dataset[idx, c("red_left", "red_right", "red_top", "red_bottom")]   
dataset[idx, "red_bottom_corr2"] <- dataset[idx, "red_bottom_corr"] * (-0.1)      
dataset[idx, "red_top_corr2"] <- dataset[idx, "red_top_corr"] * (-0.1)   
dataset[idx, "red_left_corr2"] <- dataset[idx, "red_left_corr"] - 1.4 
dataset[idx, "red_right_corr2"] <- dataset[idx, "red_right_corr"] - 1.4  





## 将归一化坐标转换为厘米坐标： 


## 屏幕尺寸：-------------------------------------- 

screen_width  <- 112   # cm
screen_height <- 142   # cm
unit_max <- 1.4        # 2D坐标范围 [-1.4, +1.4]   


scale_x <- (screen_width  / 2) / unit_max   # 56/1.4 = 40 
scale_y <- (screen_height / 2) / unit_max    



dataset_3D <- within(dataset, {    
  
  
  
  # 1) 中心点（目标格子中心）2D: unit -> cm
  target_x_cm <- target_x_corr2 * scale_x
  target_y_cm <- target_y_corr2 * scale_y 
  
  
  # 2) 点击点 2D: unit -> cm
  click_x_cm  <- click_x_corr2  * scale_x
  click_y_cm  <- click_y_corr2  * scale_y  
  
  
  # 3) 四边中点 2D：用 red_* 补齐另一个维度（见顶部假设）
  # 左/右边中点：给了 x，用中心 y
  left_x_cm   <- red_left_corr2  * scale_x 
  left_y_cm   <- target_y_corr2  * scale_y  
  
  
  right_x_cm  <- red_right_corr2 * scale_x
  right_y_cm  <- target_y_corr2  * scale_y 
  
  
  # 上/下边中点：给了 y，用中心 x
  top_x_cm    <- target_x_corr2  * scale_x
  top_y_cm    <- red_top_corr2   * scale_y
  
  
  bottom_x_cm <- target_x_corr2  * scale_x
  bottom_y_cm <- red_bottom_corr2 * scale_y 
})


summary(dataset_3D)     



# 计算3D 点击坐标： -----------------------------------------------  


## click_X:  

dataset_3D$click_X2 <- ifelse(
  !is.na(dataset_3D$click_x_cm),
  ifelse(
    dataset_3D$angle == 0,
    dataset_3D$click_x_cm,
    ifelse(
      dataset_3D$angle == 45,
      dataset_3D$distance * cos(45 * pi / 180) + dataset_3D$click_x_cm * sin(45 * pi / 180),
      ifelse(
        dataset_3D$angle == 90,
        dataset_3D$distance,
        dataset_3D$distance * sin(45 * pi / 180) + dataset_3D$click_x_cm * cos(45 * pi / 180)
      )
    )
  ),
  NA
) 



## click_Z；   

dataset_3D$click_x_cm

   dataset_3D$click_Z2 <- 
  ifelse(
    !is.na(dataset_3D$click_x_cm),
    ifelse(
      dataset_3D$angle == 0,
      dataset_3D$distance,
      ifelse(
        dataset_3D$angle == 45,
        dataset_3D$distance * sin(45 * pi / 180) - dataset_3D$click_x_cm * cos(45 * pi / 180),
        ifelse(
          dataset_3D$angle == 90,
          -(dataset_3D$click_x_cm), 
          -(dataset_3D$distance * sin(45*pi/180) + dataset_3D$click_x_cm*cos(45*pi/180))
        )
      )
    ),
    NA
  )


## click_Y:  

dataset_3D$click_Y2 <- dataset_3D$click_y_cm   





## check 3D 坐标转换是否正确： 

summary(dataset_3D$click_Z2)   
summary(dataset_3D$click_X2)          






## 镜像数据对称： ------------------------------------- 

## 1. 先筛选需要的数据列并标注为"original"：    


data_original <- dataset_3D %>% 
  mutate(side = "original") %>%            
  select(participant_id, angle, distance,v_row, v_col,
                                  correct, rt_for_calc, physical_demand, 
                                  effort, satisfaction, arm_length, 
                                  click_X2, click_Y2, click_Z2, side)       



## 2. 数据对称：  
 
dataSym <- dataset_3D %>% mutate(
  side = "mirror",  
  click_X2 = -click_X2,         ## x-axis 取反  
  click_Y2 = click_Y2,  
  click_Z2 = click_Z2) %>% 
  select(participant_id, angle, distance,v_row, v_col,
         correct, rt_for_calc, physical_demand, 
         effort, satisfaction, arm_length, 
         click_X2, click_Y2, click_Z2, side)     


## 3.Merge：  

full_data <- bind_rows(data_original, dataSym)       


## check: 

View(full_data)
nrow(data_original)  ## 62720 
nrow(full_data)     ## 125,440   



## 筛选正确反应时数据： --------------------- 

data_corr <- full_data[full_data$correct==1,]          
summary(data_corr)




## Mean Dataset: 

meanRT <- data_corr %>%   
  group_by(angle, distance, v_row, v_col,side) %>%          
  summarise(mean_rt = mean(rt_for_calc, na.rm = TRUE),  
            click_X2=mean(click_X2, na.rm = TRUE),
            click_Y2=mean(click_Y2, na.rm = TRUE),
            click_Z2=mean(click_Z2, na.rm = TRUE),  
            physical= mean(physical_demand, na.rm = TRUE),  
            effort = mean(effort, na.rm = TRUE),  
            satisfaction = mean(satisfaction, na.rm = TRUE),    
            mean_arm=mean(arm_length,na.rm=T))                       

summary(meanRT)      




## Calculate Radius: 

meanRT$r2=sqrt(meanRT$click_X2^2+meanRT$click_Y2^2+meanRT$click_Z2^2) 



## Z-score normalization of subjective data:  

meanRT <- meanRT %>% 
  ungroup() %>%
  mutate(physical_z = as.numeric(scale(1/physical)), 
         effort_z = as.numeric(scale(1/effort)),  
         satisfaction_z = as.numeric(scale(satisfaction))  
         ) 



## check 区间范围：  
summary(meanRT[, c("physical_z", "effort_z", "satisfaction_z")])      



## check变量间相关性：  

cor(meanRT[, c("physical", "effort", "satisfaction")])         ## 原始数据相关性 
cor(meanRT[, c("physical_z", "effort_z", "satisfaction_z")])   ## 标准化后数据相关性   




meanRT$subjective_mean <- rowSums(meanRT[, c("physical_z", "effort_z", "satisfaction_z")]) 
meanRT$subj_obj <- rowSums(meanRT[, c("acc_rt", "subjective_mean")])  




## calculate theta & phi:    

meanRT$theta2 <- acos(meanRT$click_Z2 / meanRT$r2)      

meanRT$phi2 <- atan2(meanRT$click_Y2,
                     meanRT$click_X2)  



## Dataset for Avg. Accuracy:  

agg <- full_data %>%       
  group_by(angle, distance, v_row, v_col, side) %>%        
  summarise(mean_acc = mean(correct, na.rm = TRUE))             


 
summary(agg)   

#agg$r2=sqrt(agg$click_X2^2+agg$click_Y2^2+agg$click_Z2^2) 

#agg$theta2 <- acos(agg$click_Z2 / agg$r2)

#agg$phi2 <- atan2(agg$click_Y2,agg$click_X2)  





meanRT <- meanRT %>% 
  left_join(agg,  
            by = c("angle","distance","v_row","v_col", "side"))   


meanRT$acc_rt=meanRT$mean_acc/meanRT$mean_rt      




### plot: theta vs. acc_rt 
 
ggplot(meanRT, aes(x = theta2, y = acc_rt)) +  
  geom_point(color = "gray", size = 1) +      # 原始数据点（可选）
  geom_smooth(method = "loess", span = 0.75,       # loess 平滑，span 控制平滑程度
              color = "black", se =FALSE,lwd=0.5) +
  geom_smooth(method = "lm",formula =y~ x+log(x), 
              span = 0.6, n=120, color = "darkred", se =FALSE) +
  labs(x = "theta", y = "ACC/RT", 
       title = " ") +    
  theme_minimal()      




## phi vs. acc_rt: 

ggplot(meanRT, aes(x = phi2, y = acc_rt)) + 
  geom_point(color = "gray", size = 1) +          
  geom_smooth(method = "loess", span = 0.3,       
              color = "black", se =FALSE,lwd=0.5) +
  geom_smooth(method = "lm",formula =y ~ cos(x)+sin(x)+cos(2*x)+sin(2*x),  
              span = 0.3, color = "darkred", se =FALSE) +
  labs(x = "phi", y = "ACC/RT", 
       title = " ") +    
  theme_minimal()   



## plot: r vs. acc_rt  
ggplot(meanRT, aes(x = r2/mean_arm, y = acc_rt)) + 
  geom_point(color = "gray", size = 1) +          # 原始数据点（可选）
  geom_smooth(method = "loess", span = 0.3,       # loess 平滑，span 控制平滑程度
              color = "black", se =FALSE,lwd=0.5) +
  geom_smooth(method = "lm",formula =y~ x+I(x^2), span = 0.3,       # loess 平滑，span 控制平滑程度
              color = "darkred", se =FALSE) +
  labs(x = "R/arm", y = "ACC/RT", 
       title = " ") +    
  theme_minimal()  





###  check theta under different phi： 、

meanRT$phi_bin <- cut(meanRT$phi2,
                      breaks = seq(-pi, pi, length.out = 3),
                      include.lowest = TRUE) 


ggplot(meanRT, aes(x = theta2, y = mean_acc)) + 
  geom_point(color = "gray", size = 1) +
  geom_smooth(method = "loess",
              span = 0.3,
              color = "black",
              se = FALSE,
              lwd = 0.5) +
  geom_smooth(method = "lm",
              formula = y ~ x + I(x^2),
              color = "darkred",
              se = FALSE) +
  facet_wrap(~phi_bin) +
  labs(x = "theta", y = "RT",
       title = "RT vs Theta under different Phi") +
  theme_minimal()



###check phi under different theta

meanRT$theta_bin <- cut(meanRT$theta2,
                        breaks = seq(min(meanRT$theta2, na.rm = TRUE),
                                     max(meanRT$theta2, na.rm = TRUE),
                                     length.out = 3),
                        include.lowest = TRUE)

meanRT$r_bin <- cut(meanRT$r2/meanRT$mean_arm,
                        breaks = seq(min(meanRT$r2/meanRT$mean_arm, na.rm = TRUE),
                                     max(meanRT$r2/meanRT$mean_arm, na.rm = TRUE),
                                     length.out = 4),
                        include.lowest = TRUE)

ggplot(meanRT, 
       aes(x = phi2, y = mean_acc)) + 
  geom_point(color = "gray", size = 1, alpha = 0.5) +
  geom_smooth(method = "loess",
              span = 0.5,
              color = "black",
              se = FALSE,
              lwd = 0.5) +
  geom_smooth(method = "lm",
              formula = y ~ x+cos(x)+cos(2*x),
              color = "darkred",
              se = FALSE) +
  facet_wrap(~ theta_bin + r_bin,nrow = 2) +
  labs(x = "Phi", y = "RT",
       title = "RT vs Phi under different Theta and r bins") +   
  theme_minimal()  



## 3D 散点热力图(for acc/rt)：   

fig <- plot_ly(meanRT,  
               x = ~ click_X2,
               y = ~ click_Y2,
               z = ~ click_Z2,  
               color = ~mean_rt,   
               colors = c("#228B22", "#FFF68F", "#8B0000"),     
               type = "scatter3d", 
               mode = "markers", 
               size=0.5,
               marker = list(opacity = 0.8)) %>% 
  layout(scene = list(
    xaxis = list(title = "click_X2", range = c(100,-100)), 
    yaxis = list(title = "click_Y2"),  
    zaxis = list(title = "click_Z2")     
  )) 

fig  




## 三维散点包络面(for acc/rt): --------------------  

summary(meanRT$acc_rt) 

meanRT <- meanRT %>% 
  mutate(rt_group = case_when(
    acc_rt <= 0.570453 ~ "Low",   
    acc_rt > 0.570453 & acc_rt <= 0.627020 ~ "Medium",   
    acc_rt > 0.627020  ~ "High"  
  ))  


table(meanRT$rt_group)  



color_map <- c(
  "Low" = "#8B0000", 
  "Medium" = "#FFF68F", 
  "High" = "#228B40"    
)




## 初始化plotly对象：   


fig <- plot_ly()

for(group in c("Low", "Medium", "High")){
  group_data <- meanRT %>% filter(rt_group == group) 
  
  # --- 包络面 ---
  if(nrow(group_data) >= 4){
    hull_idx <- convhulln(cbind(group_data$click_X2,
                                group_data$click_Y2,
                                group_data$click_Z2), output.options = TRUE) 
    n_faces <- nrow(hull_idx$hull)
    face_colors <- rep(color_map[group], n_faces)
    
    fig <- fig %>%
      add_mesh(
        x = group_data$click_X2,
        y = group_data$click_Y2,
        z = group_data$click_Z2,  
        i = hull_idx$hull[,1] - 1,
        j = hull_idx$hull[,2] - 1,
        k = hull_idx$hull[,3] - 1,
        facecolor = face_colors,
        opacity = 0.5,
        showscale = FALSE,
        name = group,           # 图例名称
        showlegend = FALSE,     # 包络面不单独占图例
        lighting = list(ambient = 0.5, diffuse = 0.8, specular = 0),
        inherit = FALSE,
        hoverinfo = "skip" 
      )
  }
  
  # --- 散点：每组单独add_trace，直接指定颜色 ---
  
  hover_text <- paste0( "x:", round(group_data$click_X2, 4), "<br>", 
                        "y:", round(group_data$click_Y2, 4), "<br>",
                        "z:", round(group_data$click_Z2, 4), "<br>",
                        "acc/rt:", round(group_data$acc_rt, 4))     
  
  fig <- fig %>%
    add_trace(
      data = group_data,
      x = ~click_X2,
      y = ~click_Y2,
      z = ~click_Z2,  
      type = "scatter3d",
      mode = "markers",
      name = group,                          # 图例
      marker = list(
        color = unname(color_map[group]),    # 直接指定hex颜色，不单独映射
        opacity = 0.8,
        size = 2
      ), 
      text = hover_text, 
      hoverinfo = "text", 
      inherit = FALSE 
    )
}


fig <- fig %>% add_trace(x = 0, y = 0, z = 0,
                         type = "scatter3d",
                         mode = "markers+text",
                         name = "Origin", 
                         marker = list(color = "black", size = 6, symbol = "cross",  opacity = 0.8),
                         text = "(0,0,0)",
                         hoverinfo = "text",
                         showlegend = FALSE, 
                         inherit = FALSE) %>% 
  layout(
    hovermode = "closest", 
    hoverdistance = 20
  )


fig 


fig <- fig %>% layout(title = list(text = "", x=0.5, xanchor="center")) 




### Model for RT:---------------------------------------------- 
# meanRT$sin_theta <- sin(meanRT$theta2)
# meanRT$cos_theta <- cos(meanRT$theta2)
# 
# meanRT$sin_phi <- sin(meanRT$phi2)
# meanRT$cos_phi <- cos(meanRT$phi2)




rt_lmm0 <- lm(subjective_mean ~ sin(theta2) + cos(theta2) + sin(phi2) + cos(phi2) + r2+    
                mean_arm ,                                 
              data = meanRT)             

summary(rt_lmm0)      
r2(rt_lmm0)      



rt_lmm1 <- lm(mean_rt ~    cos(theta2) + sin(phi2) + cos(phi2)+  r2+    
                mean_arm ,                                 
              data = meanRT)             

summary(rt_lmm1)      
r2(rt_lmm1)      



rt_lmm2 <- lm(mean_rt ~    sin(phi2) + cos(phi2) +theta2 + I(theta2^2) + r2+ I(r2^2)+    
                mean_arm ,                                 
              data = meanRT)             

summary(rt_lmm2)      
r2(rt_lmm2)       



rt_lmm3 <- lm(mean_rt ~    sin(phi2) + cos(phi2) +theta2 + I(theta2^2) + r2+ I(r2^2)+    
                r2:theta2+
                mean_arm ,                                 
              data = meanRT)             

summary(rt_lmm3)      
r2(rt_lmm3)    



rt_lmm4 <-  lm(acc_rt ~     
                 cos(2*phi2)  +
                 log(theta2)+
                 r2+ I(r2^2) +mean_arm,                                  
               data = meanRT)                


summary(rt_lmm4)      
r2(rt_lmm4)        



### 主观融合指标的最优模型： 
rt_lmm4_1 <- lm( subjective_mean~    
                   phi2+cos(2*phi2)+cos(4*phi2) +
                   theta2 + log(theta2)+
                   I(r2^2) +mean_arm,                                  
                 data = meanRT)   

summary(rt_lmm4_1) 
r2(rt_lmm4_1)



### 主/客观指标融合的最优模型：  

rt_lmm4_2 <-  lm(subj_obj~    
                    phi2+ cos(2*phi2)+cos(4*phi2) +
                    theta2 + log(theta2)+
                    I(r2^2) +mean_arm,                                  
                  data = meanRT)  

summary(rt_lmm4_2) 
r2(rt_lmm4_2) 



 
# rt_lmm4 <- lm(acc_rt ~     
#                  cos(phi2) + cos(2*phi2) +
#                 theta2 + I(theta2^2)+
#                 r2/mean_arm+ I((r2/mean_arm)^2) ,                                 
#               data = meanRT)             
# 
# summary(rt_lmm4)      
# r2(rt_lmm4)      


# rt_lmm41 <- lm(acc_rt ~    
#                 cos(phi2) + cos(2*phi2) +
#                 theta2 + log(theta2)+
#                 ratio+ I(ratio^2) ,                                 
#               data = meanRT)             
# 
# summary(rt_lmm41)      
# r2(rt_lmm41)      
# 
# anova(rt_lmm4,rt_lmm41)
# 
# BIC(rt_lmm4,rt_lmm41)

#rt_range <- range(predict(rt_lmm4), na.rm = TRUE)
#rt_range

# rt_lmm5 <-lm(mean_rt ~ 
#                sin(phi2) + cos(phi2) +cos(2*phi2) +
#                theta2 + 
#                r2/mean_arm+ I(( r2/mean_arm)^2)+ 
#                cos(phi2):(r2/mean_arm),
#              data = meanRT)    ##best
# 
# summary(rt_lmm5)      
# 
# 
# rt_lmm6 <-lm(mean_rt ~ 
#                sin(phi2) +cos(phi2)+ cos(2*phi2) +
#                theta2 + 
#                r2/mean_arm+ I( (r2/mean_arm)^2)+ 
#                cos(2*phi2):theta2:r2+
#                mean_arm,
#              data = meanRT) 
# 
# 
# 
# summary(rt_lmm6)      
# 
# anova(rt_lmm4,rt_lmm6)

# library(kernlab)

# gpr_model <- gausspr(
#   mean_rt ~ r2 + sin_theta + cos_theta + sin_phi + cos_phi,
#   data = meanRT,
#   kernel = "rbfdot"
# )
# 
# pred <- predict(gpr_model, meanRT)
# R2 <- cor(meanRT$mean_rt, pred)^2
# R2   




### 球形可视化 ---------------------------------------------


r_fixed <- 120                   


theta_seq <- seq(0.01, pi, length.out = 80) 
phi_seq   <- seq(-pi,pi, length.out = 80) 


grid <- expand.grid(
  theta2 = theta_seq,
  phi2   = phi_seq 
)


grid$r2 <- r_fixed
grid$mean_arm <- mean(meanRT$mean_arm, na.rm = TRUE)
grid$sin_phi <- sin(grid$phi2)
grid$cos_phi <- cos(grid$phi2)
grid$cos2_phi <- cos(2 * grid$phi2) 



## predict model:   

grid$sub_pred <- predict(rt_lmm4, newdata = grid)     ## change here  

summary(grid$sub_pred)          



## 球坐标转为笛卡尔坐标： 
grid$x <- grid$r2 * sin(grid$theta2) * cos(grid$phi2)  
grid$y <- grid$r2 * sin(grid$theta2) * sin(grid$phi2)    
grid$z <- grid$r2 * cos(grid$theta2)    


n_theta <- length(theta_seq)   
n_phi <- length(phi_seq)     



## 形成矩阵： 
x_mat  <- matrix(grid$x, n_theta, n_phi) 
y_mat  <- matrix(grid$y, n_theta, n_phi)
z_mat  <- matrix(grid$z, n_theta, n_phi)
rt_mat <- matrix(grid$sub_pred,n_theta, n_phi)      ## change here  
 



## 球形可视化： ----------------------------------  

fig=plot_ly(
  x = x_mat,
  y = y_mat, 
  z = z_mat, 
  surfacecolor = rt_mat,    
  type = "surface", 
  colorscale = list( 
    c(0.00, "#8B0000"),
    c(0.50, "#FEE08B"),
    c(1.00, "#66BD63")   
  ),
    cmin = 0.2,    # change     
    cmax = 1.0,     ## change  
  colorbar = list(title = "Predicted Acc/RT"),      ##change     
  hovertemplate = paste(
    "<b>Point</b><br>", 
    "x: %{x}<br>",
    "y: %{y}<br>",
    "z: %{z}<br>",
    "pred_Acc/RT: %{surfacecolor:.2f}<br>",   ## change    
    "<extra></extra>"   
  )) %>% 
  layout(
    scene = list(
      xaxis = list(range=c(150, -150), autorange=FALSE),  
      yaxis = list(range=c(-150, 150), autorange = FALSE),    
      zaxis = list(range=c(-150, 150), autorange = FALSE)),            
    title = list(
      text = sprintf("Spherical Visualization with a radius of %g cm", r_fixed),  
                     x = 0.5, 
                     xanchor= "center") 
      ) %>%   
  add_trace(type = "scatter3d",
            mode = "markers+text",
            x = c(130, -130, 0, 0, 0, 0),
            y = c(0, 0, 130, -130, 0, 0),
            z = c(0, 0, 0, 0, 130, -130),                          
            text = c("Right(+X)", "Left(-X)", "Up(+Y)", "Down(-Y)", "Front(+Z)", "Back(-Z)"),
            textposition = "top center",
            marker = list(size = 4, color = "black"), 
            showlegend = FALSE)      


fig            






## 画2D 剖面图： -----------------------------------------------    


## When Y-axis == 0 (俯视图)：       

theta_seq <- seq(0.01, pi, length.out = 300)    ## theta： 控制Z轴 
r2_seq <- seq(10, 130, length.out = 100 )    



## 右半部分(phi=0, x为正)  
grid_right <- expand.grid( 
  theta2 = theta_seq, 
  r2 = r2_seq
)
grid_right$phi2 <- 0   



## 左半部分 ： 
grid_left <- expand.grid( 
  theta2 = theta_seq, 
  r2 = r2_seq
)
grid_left$phi2 <- pi  



## 合并左右两部分： 
print(ncol(grid_right)) 
print(ncol(grid_left))

grid <- rbind(grid_left, grid_right)   



## 用与模型变量名一致： 
grid$sin_phi <- sin(grid$phi2)
grid$cos_phi <- cos(grid$phi2)
grid$cos2_phi <- cos(2 * grid$phi2)      
grid$mean_arm <- mean(meanRT$mean_arm, na.rm = TRUE)  



## When X-axis == 0 ：---------------------------------------------------------  

theta_seq <- seq(0.01, pi, length.out = 300)     ## theta： 控制Z轴 
r2_seq <- seq(10, 130, length.out = 100 )     


## 前半部分： 
grid_front <- expand.grid(theta2 = theta_seq,
                          r2= r2_seq
)
grid_front$phi2 <- pi / 2    # cos(pi/2) = 0 → x = 0 


# 后半部分（phi = -pi/2）
grid_back <- expand.grid(theta2 = theta_seq,
                         r2     = r2_seq)
grid_back$phi2 <- -pi / 2    # cos(-pi/2) = 0 → x = 0   


# 合并
grid <- rbind(grid_front, grid_back)   



## when z-axis = 0（正视图）: ----------------------------------------  

phi_seq   <- seq(0, 2*pi, length.out = 300)    
r2_seq <- seq(10,130, length.out = 100)     


grid <- expand.grid(phi2= phi_seq,
                    r2= r2_seq
)

grid$theta2 <- pi / 2    


## 用于与模型变量名一致： 
grid$sin_phi <- sin(grid$phi2)
grid$cos_phi <- cos(grid$phi2)
grid$cos2_phi <- cos(2 * grid$phi2)      
grid$mean_arm <- mean(meanRT$mean_arm, na.rm = TRUE)  



## Predict Model（所有距离下）: ----------------------------------------------- 
grid$RT_pred <- predict(rt_lmm4, newdata = grid)     


## 球坐标转为笛卡尔坐标： 
grid$x <- grid$r2 * sin(grid$theta2) * cos(grid$phi2)  
grid$y <- grid$r2 * sin(grid$theta2) * sin(grid$phi2)    
grid$z <- grid$r2 * cos(grid$theta2)    



## 计算中位数： 
summary(grid$RT_pred) 

q50 <-  0.5529  
q75 <-  0.6050


# 给每个点分配圈层
grid$ring <- case_when(
  grid$RT_pred <= q50 ~ "Difficult(<50%)",
  grid$RT_pred <= q75 ~ "Moderate(50-75%)",
  TRUE                ~ "Easy(>75%)"    
)


grid$ring <- factor(grid$ring,
                    levels = c("Difficult(<50%)",
                               "Moderate(50-75%)",
                               "Easy(>75%)"))   



# 颜色映射

color_map <- c(
  "Easy(>75%)" = "#66BD63",
  "Moderate(50-75%)"= "#FEE08B",
  "Difficult(<50%)" = "#8B0000"    
) 



# 按圈层分组画散点（X-Z平面 = 俯视）  


fig_2d <- plot_ly()  

for (ring_name in levels(grid$ring)) { 
  d <- grid[grid$ring == ring_name, ] 
  fig_2d <- fig_2d %>%
    add_trace(
      type = "scatter", 
      mode = "markers",
      x = d$x, 
      y = d$y,               
      name = ring_name,   
      marker = list(
        color = color_map[ring_name],    
        size  = 4,
        opacity = 0.7
      ),
      hovertemplate = paste0( 
        "<b>", ring_name, "</b><br>",
        "r: ", round(d$r2, 2), "<br>",
        "phi: ", round(d$phi2, 2), "<br>",             ## change here 
        "pred_Acc/RT: ", round(d$RT_pred, 3), "<br>",  ## change here  
        "<extra></extra>"       
      )  
    )
}



## "Top Down View: Spherical Visualization (Y=0)" 
## "Right Side View: Spherical Visualization (X=0)"
## "Front View: Spherical Visualization (Z=0)"  
## "X axis (Left - / Right +)" 
##  "Y axis (Down - / Up +)" 
## "Z axis (Back - / Front +)"



fig_2d <- fig_2d %>%   
  layout(
    title = list(
      text = sprintf("Front View: Spherical Visualization (Z=0)"),    ##   
      x = 0.5,
      xanchor = "center"
    ),
    xaxis = list(
      title = "X axis (Left - / Right +)",    ##  
      scaleanchor = "y",   # 锁定X/Y比例 → 正圆
      scaleratio  = 1      # 改为 > 1 则X拉伸变椭圆 
    ),
    yaxis = list(
      title =  "Y axis (Down - / Up +)"   ##     
    ),
    legend = list(title = list(text = "Pred Acc/RT")),     ## change  
    plot_bgcolor  = "white",
    paper_bgcolor = "white"  
  )


fig_2d               

















######  展示切面 

## r/arm ratio
meanRT$ratio <- meanRT$r2 / meanRT$mean_arm
range(meanRT$ratio)
summary(meanRT$ratio)

ratio_levels <- c(0.6,0.8,1,1.2,1.4,1.6)


make_theta_phi_plot <- function(ratio_fixed){
  
  grid <- expand.grid(
    theta2 = theta_seq,
    phi2   = phi_seq
  )
  
  arm_fixed <- mean(meanRT$mean_arm, na.rm = TRUE)
  
  grid$mean_arm <- arm_fixed
  grid$ratio <- ratio_fixed
  grid$r2 <- ratio_fixed * arm_fixed
  
  grid$RT_pred <- predict(rt_lmm41, newdata = grid, re.form = NA)
  
  rt_mat <- matrix(
    grid$RT_pred,
    nrow = length(theta_seq),
    ncol = length(phi_seq),
    byrow = FALSE
  )
  
  plot_ly(
    x = phi_seq,
    y = theta_seq,
    z = rt_mat,
    type = "heatmap",
    coloraxis = "coloraxis"
  )
}
p1 <- make_theta_phi_plot(0.4)
p2 <- make_theta_phi_plot(0.6)
p3 <- make_theta_phi_plot(0.8)
p4 <- make_theta_phi_plot(1)
p5 <- make_theta_phi_plot(1.2)
p6 <- make_theta_phi_plot(1.4)
p7 <- make_theta_phi_plot(1.6)

fig <- subplot(
  p1, p2, p3, p4,p5,p6,p7,
  nrows = 2,
  shareX = TRUE,
  shareY = TRUE,
  titleX = TRUE,
  titleY = TRUE
) %>%
  layout(
    xaxis  = list(title = "Phi (φ)"),
    xaxis2 = list(title = "Phi (φ)"),
    xaxis3 = list(title = "Phi (φ)"),
    xaxis4 = list(title = "Phi (φ)"),
    xaxis5  = list(title = "Phi (φ)"),
    xaxis6 = list(title = "Phi (φ)"),
    xaxis7  = list(title = "Phi (φ)"),
    yaxis  = list(title = "Theta (θ)"),
    yaxis2 = list(title = "Theta (θ)"),
    yaxis3 = list(title = "Theta (θ)"),
    yaxis4 = list(title = "Theta (θ)"),
    yaxis5  = list(title = "Theta (θ)"),
    yaxis6 = list(title = "Theta (θ)"),
    yaxis7 = list(title = "Theta (θ)"),
    coloraxis = list(
      cmin = min(predict(rt_lmm41, re.form = NA), na.rm = TRUE),
      cmax = max(predict(rt_lmm41, re.form = NA), na.rm = TRUE),
      colorscale = list(
        c(0.00, "#8B0000"),
        c(0.50, "#FEE08B"),
        c(1.00, "#66BD63")
      ),
      colorbar = list(title = "Pred RT")
    ),
    annotations = list(
      list(x = 0.06, y = 1.05, text = "Ratio = 0.4", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.38, y = 1.05, text = "Ratio = 0.6", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.62, y = 1.05, text = "Ratio = 0.8", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.94, y = 1.05, text = "Ratio = 1",   showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.06, y = 0.50, text = "Ratio = 1.2", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.38, y = 0.50, text = "Ratio = 1.4", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.62, y = 0.50, text = "Ratio = 1.6", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12))
    )
  )

fig



## distance
make_theta_phi_plot <- function(r_fixed) {
  
  grid <- expand.grid(
    theta2 = theta_seq,
    phi2   = phi_seq
  )
  
  grid$r2 <- r_fixed
  grid$mean_arm <- mean(meanRT$mean_arm, na.rm = TRUE)
  grid$r_arm_ratio <- grid$r2 / grid$mean_arm
  
  grid$RT_pred <- predict(rt_lmm4, newdata = grid, re.form = NA)
  
  rt_mat <- matrix(grid$RT_pred,
                   nrow = length(theta_seq),
                   ncol = length(phi_seq),
                   byrow = FALSE)
  
  plot_ly(
    x = phi_seq,
    y = theta_seq,
    z = rt_mat,
    type = "heatmap",
    coloraxis = "coloraxis",
    hovertemplate = paste0(
      "θ: %{y:.2f}<br>",
      "φ: %{x:.2f}<br>",
      "Pred RT: %{z:.3f}<extra></extra>"
    )
  ) %>%
    layout(
      xaxis = list(title = "Phi (φ)"),
      yaxis = list(title = "Theta (θ)")
    )
}




p_theta_phi_40  <- make_theta_phi_plot(40)
p_theta_phi_50  <- make_theta_phi_plot(50)
p_theta_phi_60  <- make_theta_phi_plot(60)
p_theta_phi_70  <- make_theta_phi_plot(70)
p_theta_phi_80  <- make_theta_phi_plot(80)
p_theta_phi_90  <- make_theta_phi_plot(90)
p_theta_phi_100  <- make_theta_phi_plot(100)
p_theta_phi_110  <- make_theta_phi_plot(110)
p_theta_phi_120  <- make_theta_phi_plot(120)
p_theta_phi_130  <- make_theta_phi_plot(130)


subplot(
  p_theta_phi_40,
  p_theta_phi_50,
  p_theta_phi_60,
  p_theta_phi_70,
  p_theta_phi_80,
  p_theta_phi_90,
  p_theta_phi_100,
  p_theta_phi_110,
  p_theta_phi_120,
  p_theta_phi_130,
  nrows = 2,
  shareX = TRUE,
  shareY = TRUE,
  titleX = TRUE,
  titleY = TRUE
) %>%
  layout(
    coloraxis = list(
      cmin = rt_range[1],
      cmax = rt_range[2],
      colorscale = list(
        c(0.00,"#8B0000"),
        c(0.50, "#FEE08B"),
        c(1.00, "#66BD63")
      ),      
      colorbar = list(title = "ACC/RT")
    )
  ) %>%
  layout(
    annotations = list(
      list(x = 0.06, y = 1.06, text = "r = 40cm", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.26, y = 1.06, text = "r = 50cm", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.52, y = 1.06, text = "r = 60cm", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.74, y = 1.06, text = "r = 70cm",   showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.94, y = 1.06, text = "r = 80cm", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      
      list(x = 0.06, y = 0.50, text = "r = 90cm", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.26, y = 0.50, text = "r = 100cm", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.52, y = 0.50, text = "r = 110cm", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.74, y = 0.50, text = "r = 120cm", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12)),
      list(x = 0.94, y = 0.50, text = "r = 130cm", showarrow = FALSE, xref = "paper", yref = "paper", font = list(size = 12))
      
      )
  )
