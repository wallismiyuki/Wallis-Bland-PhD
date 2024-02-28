library(ggplot2)
library(generics)
library(tidyverse)
library(reshape2)
library(BioGeoBEARS)
library(ape)

wd = "/Users/wbla447/Desktop/Files/Tester"
setwd(wd)

df = read.csv("Julia_speedtest.csv")

names(df)
head(df)


# By model (No DECJ for 11)

ggplot(df, aes(x=N_areas, y=Total_Time_s, color=Model, shape=Model))+
    geom_point(size=5)

# By Language
plot2 <- ggplot(df, aes(x=N_areas, y=Total_Time_s, color=Language))+
    geom_point(size=5) +
    geom_smooth(se=F) +
    labs(x = "Number of Areas", y = "Total Time (s)")

ggsave("/Users/wbla447/Desktop/Files/Tester/Plots/language.png", plot2)

# By Language and Model
plot3 <- ggplot(df, aes(x=N_areas, y=Total_Time_s, color=Language, shape=Model))+
    geom_point(size=5) +
    geom_smooth(se=F) +
    labs(x = "Number of Areas", y = "Total Time (s)")

plot3

ggsave("/Users/wbla447/Desktop/Files/Tester/Plots/language_and_model.png", plot3)


#### Biogeobears lnl compare ####

df2 = read.csv("Julia_speedtest_bgb.csv")

names(df2)

N_areas = append(df2$N_ranges, df2$N_ranges)
lnls = append(df2$Total_lnl, df2$Bgb_lnl)
neg_lnl = -lnls
program = append(rep("PhyBEARS", 18), rep("BioGeoBEARS", 18))

df3 = data.frame(N_areas, lnls, neg_lnl, program)
names = c("PhyBEARS", "BioGeoBEARS")

plot4 = ggplot(df3, aes(x = N_areas, y = neg_lnl, color = program)) +
    geom_point(size=5)+
    scale_color_manual(values = c("PhyBEARS" = "#F8766D", "BioGeoBEARS" = "#00BFC4"), labels = names) +
    labs(x = "Number of Areas", y = " - Log Likelihood (-LnL)", color = "Inference")
plot4

ggsave("/Users/wbla447/Desktop/Files/Tester/Plots/bgb_v_total_lnl.png", plot4)

plot5 = ggplot(df2, aes(x = Bgb_lnl, y = Total_lnl))+
    geom_point(size = 5) +
    xlim(-1100, 0) + ylim(-1100,0)
plot5

plot6 = ggplot(df2, aes(x = -Bgb_lnl, y = -Total_lnl))+
    geom_point(size = 5) +
    xlim(0, 1100) + ylim(0, 1100) +
    labs(x = "- (BioGeoBEARS Log Likelihood)", y = " - (PhyBEARS Log Likelihood)")
plot6
ggsave("/Users/wbla447/Desktop/Files/Tester/Plots/bgb_v_total_lnl2.png", plot6)



###### Optimizers #####

df4 = read.csv("Julia_speedtest_opt2.csv")

names(df4)
head(df4)
df4


plot7 = ggplot(df4, aes(x=N_ranges, y=Total_time, color=Opt))+
    geom_point(size=5) +
    geom_smooth(se=F) +
    labs(x = "Number of Areas", y = "Total Time (s)", color = "Solver")
plot7

ggsave("/Users/wbla447/Desktop/Files/Tester/Plots/Optimizer_compare.png", plot7)
