library(ggplot2)
library(generics)
library(tidyverse)
library(reshape2)
library(BioGeoBEARS)
library(ape)
library(binsreg)



model = "spread_only"
version = "V1"
wd1 <- (paste0("/Users/wbla447/Desktop/School/Wallis_Castor/September2023/", model, "/", version))
setwd(wd1)


# plottable_cont <- NULL

# i = 1
# wd3 <- paste0("/Users/wbla447/Desktop/School/WallisCastor/September2023V2/ss", group, "/" , sim_folder[i])
# setwd(wd3)

# branchtops_con <- read.csv("con_branchtops.csv")
# branchtops <- branchtops_con


# columns = c("Null", "A", "B", "C", "AB", "BC", "AC", "ABC", "R_nodes", "brlen", "node_age", "node_name")
# colnames(branchtops) = columns
# branchtops
# simstates <- read.table("simstates_all.txt", header = TRUE)

# CSP = NULL    # csp stands for 'Correct State Probability'


# for (i in 1:nrow(simstates)){ # nolint
#     if (simstates$x[i] == 1){CSP <- append(CSP, branchtops$Null[i])}
#     else if (simstates$x[i] == 2){CSP <- append(CSP, branchtops$A[i])}
#     else if (simstates$x[i] == 3){CSP <- append(CSP, branchtops$B[i])}
#     else if (simstates$x[i] == 4){CSP <- append(CSP, branchtops$C[i])}
#     else if (simstates$x[i] == 5){CSP <- append(CSP, branchtops$AB[i])}
#     else if (simstates$x[i] == 6){CSP <- append(CSP, branchtops$BC[i])}
#     else if (simstates$x[i] == 7){CSP <- append(CSP, branchtops$AC[i])}
#     else if (simstates$x[i] == 8){CSP <- append(CSP, branchtops$ABC[i])}
# }


# control_csp <- CSP
# brlen <- branchtops$brlen
# node_age <- branchtops$node_age

# conttable = data.frame(node_age, CSP)

# plottable_cont <- rbind(plottable_cont,conttable)

# plottable_cont
# no_outliers <- plottable_cont
# no_outliers[no_outliers < 0] <- NA
# plottable_cont_na <- na.omit(no_outliers)

# ggplot(plottable_cont_na, aes(node_age, CSP)) +
#             geom_point()





####### FULL ########



for (i in 2:8){

    group = i
    print(group)
    wd2 <- paste0("/Users/wbla447/Desktop/School/Wallis_Castor/September2023/", model , "/" , version, "/ss", group)
    setwd(wd2)

    sim_folder <- list.files()
    sim_folder

    plottable_cont <- NULL

    for(i in 1:length(sim_folder)){ # nolint
        wd3 <- paste0("/Users/wbla447/Desktop/School/Wallis_Castor/September2023/", model, "/", version, "/ss", group, "/" , sim_folder[i])
        setwd(wd3)

        branchtops_con <- read.csv("con_branchtops_7000.csv")
        branchtops_spread <- read.csv("branchtops_inf.csv")

        branchbottoms_con <- read.csv("con_branchbottoms_7000.csv")
        branchbottoms_spread <- read.csv("branchbottoms_inf.csv")

        inftable_con <- read.csv("con_inferredtable_7000.csv")
        inftable_spread <- read.csv("inferredtable_inf.csv")

        branchtops_con_7000 <- read.csv("con_branchtops_7000.csv")
        branchtops_spread_7000 <- read.csv("branchtops_inf_7000.csv")

        branchbottoms_con_7000 <- read.csv("con_branchbottoms_7000.csv")
        branchbottoms_spread_7000 <- read.csv("branchbottoms_inf_7000.csv")

        inftable_con_7000 <- read.csv("con_inferredtable_7000.csv")
        inftable_spread_7000 <- read.csv("inferredtable_inf_7000.csv")


        ## We will mostly be using branch tops - as that is the actual nodes and end points

        simstates <- read.table("simstates_fixed.txt", header = TRUE)
        colnames(simstates) = "x"


        ### CHANGE THIS BASED ON CHICH ONE YOU ARE MAKING! ###

        ###########################
        ######## CONTROL ##########
        ###########################

        branchtops <- branchtops_con_7000

        columns = c("Null", "A", "B", "C", "AB", "BC", "AC", "ABC", "R_nodes", "brlen", "node_age", "node_name")
        colnames(branchtops) = columns
        branchtops

        CSP = NULL    # csp stands for 'Correct State Probability'


        for (i in 1:nrow(simstates)){ # nolint
            if (simstates$x[i] == 1){CSP <- append(CSP, branchtops$Null[i])}
            else if (simstates$x[i] == 2){CSP <- append(CSP, branchtops$A[i])}
            else if (simstates$x[i] == 3){CSP <- append(CSP, branchtops$B[i])}
            else if (simstates$x[i] == 4){CSP <- append(CSP, branchtops$C[i])}
            else if (simstates$x[i] == 5){CSP <- append(CSP, branchtops$AB[i])}
            else if (simstates$x[i] == 6){CSP <- append(CSP, branchtops$BC[i])}
            else if (simstates$x[i] == 7){CSP <- append(CSP, branchtops$AC[i])}
            else if (simstates$x[i] == 8){CSP <- append(CSP, branchtops$ABC[i])}
        }


        control_csp <- CSP
        brlen <- branchtops$brlen
        node_age <- branchtops$node_age

        conttable = data.frame(node_age, CSP)


        ###########################
        ########## spread ##########
        ###########################

        branchtops <- branchtops_spread_7000

        colnames(branchtops) = columns
        branchtops

        CSP = NULL    # csp stands for 'Correct State Probability'


        for (i in 1:nrow(simstates)){ # nolint
            if (simstates$x[i] == 1){CSP <- append(CSP, branchtops$Null[i])}
            else if (simstates$x[i] == 2){CSP <- append(CSP, branchtops$A[i])}
            else if (simstates$x[i] == 3){CSP <- append(CSP, branchtops$B[i])}
            else if (simstates$x[i] == 4){CSP <- append(CSP, branchtops$C[i])}
            else if (simstates$x[i] == 5){CSP <- append(CSP, branchtops$AB[i])}
            else if (simstates$x[i] == 6){CSP <- append(CSP, branchtops$BC[i])}
            else if (simstates$x[i] == 7){CSP <- append(CSP, branchtops$AC[i])}
            else if (simstates$x[i] == 8){CSP <- append(CSP, branchtops$ABC[i])}
        }


        spread_csp <- CSP
        brlen <- branchtops$brlen
        node_age <- branchtops$node_age

        spreadtable = data.frame(node_age, CSP)

        ##########################
        ### NOW MAKE BIG TABLE ###
        ##########################


        ## Differentiate between models
        conttable$Model = "cont"
        spreadtable$Model = "spread"


        bigtable <- rbind(conttable, spreadtable)
        no_outliers <- bigtable
        no_outliers[no_outliers < 0] <- NA
        bigtable_long <- na.omit(no_outliers)


        # ggplot(bigtable_long, aes(brlen, CSP, color = Model)) +
        #     geom_point()

        write.csv(bigtable_long, "bigtable_long_7000_fixed.csv")
        write.csv(conttable, "conttable_7000_fixed.csv")
        write.csv(spreadtable, "spreadtable_7000_fixed.csv")

        bigtable_wide <- data.frame(node_age, control_csp, spread_csp)
        write.csv(bigtable_wide, "bigtable_wide_7000_fixed.csv")

        ###############################
        ### NOW FOR THE OLD VERSION ###
        ###############################

        branchtops <- branchtops_con

        columns = c("Null", "A", "B", "C", "AB", "BC", "AC", "ABC", "R_nodes", "brlen", "node_age", "node_name")
        colnames(branchtops) = columns
        branchtops

        CSP = NULL    # csp stands for 'Correct State Probability'


        for (i in 1:nrow(simstates)){ # nolint
            if (simstates$x[i] == 1){CSP <- append(CSP, branchtops$Null[i])}
            else if (simstates$x[i] == 2){CSP <- append(CSP, branchtops$A[i])}
            else if (simstates$x[i] == 3){CSP <- append(CSP, branchtops$B[i])}
            else if (simstates$x[i] == 4){CSP <- append(CSP, branchtops$C[i])}
            else if (simstates$x[i] == 5){CSP <- append(CSP, branchtops$AB[i])}
            else if (simstates$x[i] == 6){CSP <- append(CSP, branchtops$BC[i])}
            else if (simstates$x[i] == 7){CSP <- append(CSP, branchtops$AC[i])}
            else if (simstates$x[i] == 8){CSP <- append(CSP, branchtops$ABC[i])}
        }


        control_csp <- CSP
        brlen <- branchtops$brlen
        node_age <- branchtops$node_age

        conttable = data.frame(node_age, CSP)


        # ###########################
        # ########## spread ##########
        # ###########################

        branchtops <- branchtops_spread

        colnames(branchtops) = columns
        branchtops

        CSP = NULL    # csp stands for 'Correct State Probability'


        for (i in 1:nrow(simstates)){ # nolint
            if (simstates$x[i] == 1){CSP <- append(CSP, branchtops$Null[i])}
            else if (simstates$x[i] == 2){CSP <- append(CSP, branchtops$A[i])}
            else if (simstates$x[i] == 3){CSP <- append(CSP, branchtops$B[i])}
            else if (simstates$x[i] == 4){CSP <- append(CSP, branchtops$C[i])}
            else if (simstates$x[i] == 5){CSP <- append(CSP, branchtops$AB[i])}
            else if (simstates$x[i] == 6){CSP <- append(CSP, branchtops$BC[i])}
            else if (simstates$x[i] == 7){CSP <- append(CSP, branchtops$AC[i])}
            else if (simstates$x[i] == 8){CSP <- append(CSP, branchtops$ABC[i])}
        }


        spread_csp <- CSP
        brlen <- branchtops$brlen
        node_age <- branchtops$node_age

        spreadtable = data.frame(node_age, CSP)

        # ##########################
        # ### NOW MAKE BIG TABLE ###
        # ##########################


        # ## Differentiate between models
        conttable$Model = "cont"
        spreadtable$Model = "spread"


        bigtable <- rbind(conttable, spreadtable)
        no_outliers <- bigtable
        no_outliers[no_outliers < 0] <- NA
        bigtable_long <- na.omit(no_outliers)


        # # ggplot(bigtable_long, aes(brlen, CSP, color = Model)) +
        # #     geom_point()

        write.csv(bigtable_long, "bigtable_long_fixed.csv")
        write.csv(conttable, "conttable_fixed.csv")
        write.csv(spreadtable, "spreadtable_fixed.csv")

        bigtable_wide <- data.frame(node_age, control_csp, spread_csp)
        write.csv(bigtable_wide, "bigtable_wide_fixed.csv")


    }
}

#######################################
# NOW BRING THEM ALL BACK IN AND PLOT #
#######################################

group = 8
version = "V1"
wd2 <- paste0("/Users/wbla447/Desktop/School/Wallis_Castor/September2023/", model , "/" , version, "/ss", group)
setwd(wd2)

sim_folder <- list.files()

plot_table <- NULL
plot_table_7000 <- NULL
models_table <- NULL
models_table_7000 <- NULL


for(i in 1:length(sim_folder)){ # nolint
    wd3 <- paste0("/Users/wbla447/Desktop/School/Wallis_Castor/September2023/", model, "/", version, "/ss", group, "/" , sim_folder[i])
    setwd(wd3)
    table <- read.csv("bigtable_long_7000_fixed.csv")
    table_model <- read.csv("bigtable_wide_7000_fixed.csv")


    plot_table_7000 <- rbind(plot_table_7000, table)
    models_table_7000 <- rbind(models_table_7000, table_model)
}


###### 9000 version ######
# nolint
for(i in 1:length(sim_folder)){ # nolint
    wd3 <- paste0("/Users/wbla447/Desktop/School/Wallis_Castor/September2023/", model, "/", version, "/ss", group, "/" , sim_folder[i])
    setwd(wd3) 
    table <- read.csv("bigtable_long_fixed.csv")
    table_model <- read.csv("bigtable_wide_fixed.csv")

    plot_table <- rbind(plot_table, table)
    models_table <- rbind(models_table, table_model)
}




#######
# Make sure to change this for every start state!
#######


#plot_table8_7000 <- plot_table_7000
#plot_table8 <- plot_table
#models_table8_7000 <- models_table_7000
#models_table8 <- models_table


wd1 <- (paste0("/Users/wbla447/Desktop/School/Wallis_Castor/September2023/", model, "/", version))
setwd(wd1)

write.csv(plot_table_7000, paste0("fixed_plot_table", group, "_7000.csv"))
write.csv(plot_table, paste0("fixed_plot_table", group, ".csv"))
write.csv(models_table_7000, paste0("fixed_models_table", group, "_7000.csv"))
write.csv(models_table, paste0("fixed_models_table", group, ".csv"))


################
# BIG BOY PLOT #
################
### REMEMBER TO CHANGE NUMBER IN THE PRINTOUT!


#plot_table = read.csv("/Users/wbla447/Desktop/School/Wallis_Castor/September2023/spread_only/V1/plot_table8_7000.csv")

model = "spread_only"
#list = c("V1", "V2", "V3", "V4", "V5")
list = c("V1")

for (version in list){

    for (i in 2:8){
        group = i
        plot_table = read.csv(paste0("/Users/wbla447/Desktop/WallisCastor_Git/", model, "/", version, "/plot_table", group, "_7000.csv"))
        testing <- as.data.frame(plot_table)[,-1]
        colnames(testing)
        colnames(plot_table)

        plot1 <- ggplot(plot_table, aes(node_age, CSP, color=Model)) +
                    geom_point() +
                    labs(x = "Node Age (mya)", y = "CSP", title = LETTERS[group - 1]) 

        plot1

        plot2 <- ggplot(plot_table, aes(node_age, CSP)) +
                    geom_point() + facet_wrap( ~ Model) +
                    labs(x = "Node Age (mya)", y = "CSP", title = LETTERS[group - 1])

        plot2

        pdf(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/points/", version, "ss", group, ".pdf"))
        plot1
        plot2
        dev.off()

        plot4 <- ggplot(plot_table, aes(node_age, CSP)) +
                geom_point() + facet_wrap( ~ Model) + geom_smooth(method = 'lm') +
                labs(x = "Node Age (mya)", y = "CSP", title = LETTERS[group - 1])

        plot4

        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lm/", version, "ss", group, "_gglm.png"), plot4)
        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/points/", version, "ss", group, ".png"), plot1)
        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/points/", version, "ss", group, "_sep.png"),plot2)
    }
}



for (version in list){

    for (i in 2:8){
        group = i
        plot_table = read.csv(paste0("/Users/wbla447/Desktop/WallisCastor_Git/", model, "/", version, "/plot_table", group, "_7000.csv"))
        testing <- as.data.frame(plot_table)[,-1]
        colnames(testing)
        colnames(plot_table)
        
        plot3 <- ggplot(plot_table, aes(node_age, CSP, color=Model)) +
            geom_point() + geom_smooth(method = 'lm') + ylim(0,1) +
            labs(x = "Node Age (mya)", y = "CSP", title = LETTERS[group - 1]) +
            scale_color_manual(values = c(cont = "#F8766D", spread = "#00BFC4"), labels = c("Control", "Spread"), name = "Model")

        plot3

        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lm/", version, "ss", group, "_gglm_color.png"), plot3)
    }
}

plot3_5 <- ggplot(plot_table2_7000, aes(node_age, CSP, color=Model)) +
            geom_point() + geom_smooth(method = 'lm') + ylim(0,1) +

plot3_5


for (version in list){

    for (i in 2:8){
        group = i
        plot_table = read.csv(paste0("/Users/wbla447/Desktop/WallisCastor_Git/", model, "/", version, "/plot_table", group, "_7000.csv"))
        testing <- as.data.frame(plot_table)[,-1]
        colnames(testing)
        colnames(plot_table)
        
    plot4 <- ggplot(plot_table, aes(node_age, CSP)) +
            geom_point() + facet_wrap( ~ Model) + geom_smooth(method = 'lm') + 
            ylim(0,1) +
            labs(x = "Node Age (mya)", y = "CSP", title = LETTERS[group - 1])

    plot4

    ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lm/", version, "ss", group, "_gglm.png"), plot4)
    }
}

pdf(file = paste0("/Users/wbla447/Desktop/School/Wallis_Castor/September2023/", model, "/", version, "/Plots/", version, "lm", 2, ".pdf"))
plot3
plot3_5
dev.off()

plot4 <- ggplot(plot_table, aes(node_age, CSP)) +
            geom_point() + facet_wrap( ~ Model) + geom_smooth(method = 'lm')

plot4



##################
# FOR EACH MODEL #
##################

models_table
models_table <- as.data.frame(models_table)[,-1]
colnames(models_table)

control_plot <- ggplot(models_table, aes(node_age, control_csp)) + geom_point()
spread_plot <- ggplot(models_table, aes(node_age, spread_csp)) + geom_point()







#########################
# TESTING OLD PLOT TREE #
#########################

version = "V3"
group = 8
simnum = "001"


setwd(paste0("/Users/wbla447/Desktop/School/WallisCastor/September2023", version, "/ss", group, "/sim_", simnum))
tree = read.tree(file="tree_wFossils.newick")
tree

states = read.table("simstates_living.txt")

# labels and colors
tmplabels = c("_", "A", "B", "C", "AB", "BC", "AC", "ABC")
tmpcols = c("white", "red", "orange", "yellow", "green", "blue", "purple", "black")

# internal nodes
intnums = (length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode)

# Tip node numbers
tipnums =  1:length(tree$tip.label)

plot(tree, label.offset=1); axisPhylo(); mtext(text="Millions of years ago (Ma)", side=1, line=2)
nodelabels(text=tmplabels[c(states[intnums,])], tip=tipnums, bg=tmpcols[c(states[intnums,])])
tiplabels(text=tmplabels[c(states[tipnums,])], tip=tipnums, bg=tmpcols[c(states[tipnums,])])
tiplabels()



#####################
# LINEAR REGRESSION #
#####################

version = "V2"
model = "spread_only"
group = 2

models_table <- read.csv(paste0("/Users/wbla447/Desktop/WallisCastor_Git/", model, "/", version, "/models_table", group, "_7000.csv"))
#models_table <- as.data.frame(models_table)[,-c(1,2)]
#colnames(models_table)

# control_plot <- ggplot(models_table, aes(node_age, control_csp)) + geom_point()
# ss_plot <- ggplot(models_table, aes(node_age, ss_csp)) + geom_point()
# split_plot <- ggplot(models_table, aes(node_age, split_csp)) + geom_point()
# spread_plot <- ggplot(models_table, aes(node_age, spread_csp)) + geom_point()

control_lm <- linear_regression_plot(models_table$node_age,models_table$control_csp)
spread_lm <- linear_regression_plot(models_table$node_age, models_table$spread_csp)




list = c("V1", "V2", "V3", "V4", "V5")
model = "spread_only"

for (version in list){
    for (group in 2:8){
        models_table_7000 <- read.csv(paste0("/Users/wbla447/Desktop/WallisCastor_Git/", model, "/", version, "/models_table", group, "_7000.csv"))
        models_table_7000 <- as.data.frame(models_table_7000)[,-c(1,2)]

        pdf(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lm/", version, "ss", group, ".pdf"))

        control_lm_7000 <- linear_regression_plot(models_table_7000$node_age, models_table_7000$control_csp)
        title(paste0("Control Model, start state", group))

        split_lm_7000 <- linear_regression_plot(models_table_7000$node_age, models_table_7000$spread_csp)
        title(paste0("Spread Model, start state", group))

        dev.off()
    }
}



# bin
# 0 - 1 mil
# take the average/95%
# then just do little box plots


### does this tell us something about how messy connectivity is when things are close?


#########################
#  LINEAR REGRESSION WITHOUT 1S AND 0S
#########################

coef <- NULL
points <- NULL
confident <- NULL
unconfident <- NULL


list = c("V1")
model = "spread_only"

for(version in list){
    for (group in 2:8){

        models_table <- read.csv(paste0("/Users/wbla447/Desktop/WallisCastor_Git/", model, "/", version, "/fixed_models_table", group, "_7000.csv"))
        models_table <- as.data.frame(models_table)[,-c(1,2)]

        TF1 = models_table$control_csp < 0.99999999
        TF2 = models_table$control_csp > 0.00000001

        control_csp_fractions = models_table$control_csp[(TF1 + TF2)==2]
        control_csp_fractions_node_ages = models_table$node_age[(TF1 + TF2)==2]

        lm_df1_control_CSP = linear_regression_plot(x=control_csp_fractions_node_ages, y=control_csp_fractions)
        title(paste0("start state", group, ", control_csp, no 0.0s and 1.0s"))

        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/No1_0/LM_", version, "ss", group, "control_9000.png"))
        lm_df1_control_CSP = linear_regression_plot(x=control_csp_fractions_node_ages, y=control_csp_fractions)
        title(paste0("start state", group, ", control_csp, no 0.0s and 1.0s") )
        dev.off()

        coef  <-rbind(coef, lm_df1_control_CSP$coefficients)

        TF7 = models_table$spread_csp < 0.99999999
        TF8 = models_table$spread_csp > 0.00000001

        spread_csp_fractions = models_table$spread_csp[(TF7 + TF8)==2]
        spread_csp_fractions_node_ages = models_table$node_age[(TF7 + TF8)==2]

        lm_df1_spread_CSP = linear_regression_plot(x=spread_csp_fractions_node_ages, y=spread_csp_fractions)
        title(paste0("start state", group, ", spread_csp, no 0.0s and 1.0s"))

        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/No1_0/LM_", version, "ss", group, "spread_9000.png"))
        lm_df1_spread_CSP = linear_regression_plot(x=spread_csp_fractions_node_ages, y=spread_csp_fractions)
        title(paste0("start state", group, ", spread_csp, no 0.0s and 1.0s"))
        dev.off()

        coef <- rbind(coef, lm_df1_spread_CSP$coefficients)

        points <- rbind(points, length(control_csp_fractions))
        points <- rbind(points, length(spread_csp_fractions))

        # 1.00?

        confident <- rbind(confident, table(TF1)["FALSE"])
        confident <- rbind(confident, table(TF7)["FALSE"])


        # 0.00?

        unconfident <- rbind(unconfident, table(TF2)["FALSE"])
        unconfident <- rbind(unconfident, table(TF8)["FALSE"])
    }
}


coef
points
confident
unconfident




### 
# Print to pdf
###

pdf(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/no1_0/", version, "ss", group, ".pdf"))

lm_df1_control_CSP = linear_regression_plot(x=control_csp_fractions_node_ages, y=control_csp_fractions)
title(paste0("start state", group, ", control_csp, no 0.0s and 1.0s"))

lm_df1_spread_CSP = linear_regression_plot(x=spread_csp_fractions_node_ages, y=spread_csp_fractions)
title(paste0("start state", group, ", spread_csp, no 0.0s and 1.0s"))

dev.off()




###################
# Root State Plot #
###################

model = "spread_only"
#version = "V1"


list = c("V1", "V2", "V3", "V4", "V5")

for(version in list){
    for(group in 2:8){
        wd = paste0("/Users/wbla447/Desktop/WallisCastor_Git/", model, "/", version, "/")
        setwd(wd)
        num = c("A", "B", "C", "AB", "AC", "BC", "ABC")
        num_title = c("A", "B", "C", "AB", "BC", "AC", "ABC") #because it turns out ss 6 and 7 should be reversed, theyre being fed in in the wrong order, so just make sure theyre labelled correctly after!

        conttable = read.csv(paste0("bigrootstateprobs", group, "_control_7000.csv"))
        spreadtable = read.csv(paste0("bigrootstateprobs", group, "inf_7000.csv"))

        conttable
        spreadtable


        ## Differentiate between models

        conttable$Model = "cont"
        spreadtable$Model = "spread"



        bigtable = rbind(conttable, spreadtable)
        colnames(bigtable)

        bigtable_NNSL  = bigtable[-c(1:2, 10:11)] # remove the liklihoods, Null AND SIM NUMBER
        summary(bigtable_NNSL)

        testdf = melt(bigtable_NNSL, id.vars="Model")
        summary(testdf)
        names = c("Control", "Spread")

        RootStatePlot = ggplot(testdf, aes(x=variable, y=value, color = Model)) +
            geom_boxplot() +
            labs(x = "Root State", y = "Root Probability", title = num[group - 1]) +
            scale_color_manual(values = c(cont = "#F8766D", spread = "#00BFC4"), labels = names)
        #RootStatePlot


        # pdf(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/rootstate/", version, "ss", group, ".pdf"))
        # RootStatePlot
        # dev.off()

        ###########################
        # Take only the highest one?
        ############################

        bigtable_count = bigtable
        bigtable_count$biggest = ""

        bigtable_count

        for (i in 1:length(bigtable_count$biggest)){
            bigtable_count$biggest[i] = "A"
            if (bigtable_count$A[i] < bigtable_count$B[i]){bigtable_count$biggest[i] = "B"}
            if (bigtable_count$A[i] < bigtable_count$C[i] && bigtable_count$B[i] < bigtable_count$C[i]) {bigtable_count$biggest[i] = "C"}
            if (bigtable_count$A[i] < bigtable_count$AB[i] && bigtable_count$B[i] < bigtable_count$AB[i] && bigtable_count$C[i] < bigtable_count$AB[i]) {bigtable_count$biggest[i] = "AB"}
            if (bigtable_count$A[i] < bigtable_count$BC[i] && bigtable_count$B[i] < bigtable_count$BC[i] && bigtable_count$C[i] < bigtable_count$BC[i] && bigtable_count$AB[i] < bigtable_count$BC[i]) {bigtable_count$biggest[i] = "AC"} # we needed to switch the AC and BC counts! So if old BC is the biggest, that's actually AC
            if (bigtable_count$A[i] < bigtable_count$AC[i] && bigtable_count$B[i] < bigtable_count$AC[i] && bigtable_count$C[i] < bigtable_count$AC[i] && bigtable_count$AB[i] < bigtable_count$AC[i] && bigtable_count$BC[i] < bigtable_count$AC[i]) {bigtable_count$biggest[i] = "BC"}
            if (bigtable_count$A[i] < bigtable_count$ABC[i] && bigtable_count$B[i] < bigtable_count$ABC[i] && bigtable_count$C[i] < bigtable_count$ABC[i] && bigtable_count$AB[i] < bigtable_count$ABC[i] && bigtable_count$BC[i] < bigtable_count$ABC[i] && bigtable_count$AC[i] < bigtable_count$ABC[i]) {bigtable_count$biggest[i] = "ABC"}

        }

        summary(bigtable_count)

        tester = ggplot(bigtable_count, aes(x = biggest)) +
            geom_bar()
        tester
        

        RootStatePlot_count = ggplot(bigtable_count, aes(x = biggest, fill = Model)) +
            geom_bar(position = "dodge") +
            scale_x_discrete(limits = num) +
            labs(x = "Highest Probability Rootstate", y = "Count", title = num_title[group - 1]) +
            scale_fill_manual(values = c(cont = "#F8766D", spread = "#00BFC4"), labels = names) +
            ylim(0,350)

        RootStatePlot_count


        pdf(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/rootstate/", version, "ss", group, ".pdf"))
        RootStatePlot
        RootStatePlot_count
        dev.off()

        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/rootstate/", version, "ss", group, ".png"), RootStatePlot)
        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/rootstate/", version, "ss", group, "_count.png"), RootStatePlot_count)
    }
}


###################
## BINNING PLOTS ##
###################


model = "spread_only"
#list = c("V1", "V2", "V3", "V4", "V5")
list = c("V1")

for(version in list){
    for(group in 2:8){

        wd = (paste0("/Users/wbla447/Desktop/WallisCastor_Git/", model, "/", version))
        setwd(wd)

        fn1 = paste0("fixed_models_table", group, "_7000.csv")
        fn2 = paste0("fixed_plot_table", group, "_7000.csv")

        models_table = read.csv(fn1)
        models_table <- as.data.frame(models_table)[,-c(1,2)]
        plot_table = read.csv(fn2)

        head(models_table)
        head(plot_table)


        #####
        ## Removing 0 / 1
        #####

        TF1 = models_table$control_csp < 0.99999999
        TF2 = models_table$control_csp > 0.00000001

        control_csp_fractions = models_table$control_csp[(TF1 + TF2)==2]
        control_csp_fractions_node_ages = models_table$node_age[(TF1 + TF2)==2]

        lm_df1_control_CSP = linear_regression_plot(x=control_csp_fractions_node_ages, y=control_csp_fractions)
        title(paste0("start state", group, ", control_csp, no 0.0s and 1.0s"))

        TF7 = models_table$spread_csp < 0.99999999
        TF8 = models_table$spread_csp > 0.00000001

        spread_csp_fractions = models_table$spread_csp[(TF7 + TF8)==2]
        spread_csp_fractions_node_ages = models_table$node_age[(TF7 + TF8)==2]

        lm_df1_spread_CSP = linear_regression_plot(x=spread_csp_fractions_node_ages, y=spread_csp_fractions)
        title(paste0("start state", group, ", spread_csp, no 0.0s and 1.0s"))

        # Binned plots
        # control bins

        # Bins plot *including* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)


        ###### Rename things! #######

        Control = NULL
        Spread = NULL

        Control$CSP <- models_table$control_csp
        Control$Node.Age <- models_table$node_age
        Spread$CSP <- models_table$spread_csp
        Spread$Node.Age <- models_table$node_age

        Control.CSP_NO_1_0 = control_csp_fractions 
        Control.NODE_AGE = control_csp_fractions_node_ages
        Spread.CSP_NO_1_0 = spread_csp_fractions
        Spread.NODE_AGE = spread_csp_fractions_node_ages 


        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/binned/", version, "ss", group, "control.png"))
        binsreg(y=Control$CSP, x=Control$Node.Age, 
                line=NULL, 
                plotyrange=c(0,1))
        dev.off()

        # with 95% CI
        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/binned/", version, "ss", group, "control_90CI.png"))
        binsreg(y=Control$CSP, x=Control$Node.Age, line=TRUE, ci=TRUE, plotyrange=c(0,1))
        dev.off()

        # Bins plot *excluding* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/no1_0/", version, "ss", group, "control_binned.png"))
        binsreg(y=Control.CSP_NO_1_0, x=Control.NODE_AGE, line=NULL, plotyrange=c(0,1))
        dev.off()

        # with 95% CI
        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/no1_0/", version, "ss", group, "control_binned90.png"))
        binsreg(y=Control.CSP_NO_1_0, x=Control.NODE_AGE, line=TRUE, ci=TRUE, plotyrange=c(0,1))
        dev.off()

        # 3. spread_csp

        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/binned/", version, "ss", group, "spread.png"))
        # Bins plot *including* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
        binsreg(y=Spread$CSP, x=Spread$Node.Age, 
            line=NULL, 
            noplot=FALSE, 
            plotyrange=c(0,1))

        dev.off()

        # with 95% CI
        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/binned/", version, "ss", group, "spread_90.png"))
        binsreg(y=Spread$CSP, x=Spread$Node.Age, line=TRUE, ci=TRUE, noplot=FALSE, plotyrange=c(0,1))
        dev.off()

        # Bins plot *excluding* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/no1_0/", version, "ss", group, "spread_binned.png"))
        binsreg(y=Spread.CSP_NO_1_0, x=Spread.NODE_AGE, line=NULL, plotyrange=c(0,1))
        dev.off()

        # with 95% CI
        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/no1_0/", version, "ss", group, "spread_binned90.png"))
        binsreg(y=Spread.CSP_NO_1_0, x=Spread.NODE_AGE, line=TRUE, ci=TRUE, plotyrange=c(0,1))
        dev.off()




        # 3. plot_table$CSP


        TF1 = plot_table$CSP < 0.99999999
        TF2 = plot_table$CSP > 0.00000001

        CSP_fractions = plot_table$CSP[(TF1 + TF2)==2]
        CSP_fractions_node_ages = plot_table$node_age[(TF1 + TF2)==2]
        CSP_fractions_model = plot_table$Model[(TF1 + TF2)==2]

        df_big = data.frame(CSP_fractions, CSP_fractions_node_ages, CSP_fractions_model) 

        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/no1_0/", version, "ss", group, "_Biglm.png"))
        lm_plot_table_CSP = linear_regression_plot(x=CSP_fractions_node_ages, y=CSP_fractions, xlab="Node Age", ylab="CSP")
        title("All points, no 0.0s and 1.0s")
        dev.off()

        png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/no1_0/", version, "ss", group, "_Control.png"))
        lm_df1_control_CSP = linear_regression_plot(x=control_csp_fractions_node_ages, y=control_csp_fractions, xlab="Node Age", ylab="CSP")
        title(paste0("start state", group, ", control_csp, no 0.0s and 1.0s"))
        dev.off()


        png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/no1_0/", version, "ss", group, "_Spread.png"))
        lm_df1_spread_CSP = linear_regression_plot(x=spread_csp_fractions_node_ages, y=spread_csp_fractions, xlab="Node Age", ylab="CSP")
        title(paste0("start state", group, ", spread_csp, no 0.0s and 1.0s"))
        dev.off()

        names = c("Control", "Spread")
        plot <-ggplot(df_big, aes(CSP_fractions_node_ages, CSP_fractions, color = CSP_fractions_model)) +
                geom_point() + geom_smooth(method = 'lm') +
                labs(x = "Node Age (mya)", y = "Correc CSP", title = LETTERS[group - 1]) +
                scale_color_manual(values = c(cont = "#F8766D", spread = "#00BFC4"), labels = names, name = "Model")


        plot

        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/no1_0/", version, "ss", group, "_gglm.png"), plot)
    }
}

model = "spread_only"
list = c("V1")

for(version in list){
    for(group in 2:8){

        wd = (paste0("/Users/wbla447/Desktop/WallisCastor_Git/", model, "/", version))
        setwd(wd)

        fn1 = paste0("fixed_plot_table", group, "_7000.csv")
        fn2 = paste0("fixed_plot_table", group, ".csv")


        head(read.csv(fn1))
        head(read.csv(fn2))
        hist(plot_table$CSP)


        #####
        ## Removing 0 / 1
        #####

        # Binned plots
        # control bins

        # Bins plot *including* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)

        # 3. plot_table$CSP

        plot_table = read.csv(fn1) #start with original 7000
        TF1 = plot_table$CSP < 0.99999999
        TF2 = plot_table$CSP > 0.00000001

        CSP_fractions = plot_table$CSP[(TF1 + TF2)==2]
        CSP_fractions_node_ages = plot_table$node_age[(TF1 + TF2)==2]
        CSP_fractions_model = plot_table$Model[(TF1 + TF2)==2]

        df_big_7000 = data.frame(CSP_fractions, CSP_fractions_node_ages, CSP_fractions_model) 
        df_big_7000$geo_iso = 7000
        head(df_big_7000)

        # then do for 9000

        plot_table = read.csv(fn2)
        TF1 = plot_table$CSP < 0.99999999
        TF2 = plot_table$CSP > 0.00000001

        CSP_fractions = plot_table$CSP[(TF1 + TF2)==2]
        CSP_fractions_node_ages = plot_table$node_age[(TF1 + TF2)==2]
        CSP_fractions_model = plot_table$Model[(TF1 + TF2)==2]

        df_big_9000 = data.frame(CSP_fractions, CSP_fractions_node_ages, CSP_fractions_model) 
        df_big_9000$geo_iso = 9000
        head(df_big_9000)


        df_big = rbind(df_big_7000, df_big_9000)


        names = c("Control", "Spread")

        plot <-ggplot(df_big, aes(CSP_fractions_node_ages, CSP_fractions, color = CSP_fractions_model)) +
                geom_point(size = 1) + geom_smooth(method = 'lm') +
                labs(x = "Node Age (mya)", y = "CSP", title = LETTERS[group - 1]) +
                xlim(0,25) +
                facet_wrap(~geo_iso)+
                scale_color_manual(values = c(cont = "#F8766D", spread = "#00BFC4"), labels = names, name = "Model")


        plot

        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/9000v7000/", version, "ss", group, "_gglm.png"), plot)


    }
}

###########
# NO TIPS
###########

model = "spread_only"
coef <- NULL
list = c("V1")
version = "V1"

#for(version in list){
   for(group in 2:8){
        print(group)      
        wd = (paste0("/Users/wbla447/Desktop/WallisCastor_Git/", model, "/", version))
        setwd(wd)

        fn1 = paste0("fixed_models_table", group, "_7000.csv")
        fn2 = paste0("fixed_plot_table", group, "_7000.csv")

        models_table = read.csv(fn1)
        models_table <- as.data.frame(models_table)[,-c(1,2)]
        plot_table = read.csv(fn2)

        #head(models_table)
        #head(plot_table)

        #hist(models_table$control_csp)

        #hist(plot_table$CSP)


        #####
        ## Removing tips
        #####

        TF1 = models_table$node_age > 0.00000001
        TF1
        control_csp_fractions = models_table$control_csp[TF1 == TRUE]
        control_csp_fractions_node_ages = models_table$node_age[TF1 == TRUE]

        lm_df1_control_CSP = linear_regression_plot(x=control_csp_fractions_node_ages, y=control_csp_fractions)
        title(paste0("start state", group, ", control model, no tips"))

        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/LM_", version, "ss", group, "control.png"))
        lm_df1_control_CSP = linear_regression_plot(x=control_csp_fractions_node_ages, y=control_csp_fractions)
        title(paste0("start state", group, ", control model, no tips"))
        dev.off()

        coef  <-rbind(coef, lm_df1_control_CSP$coefficients)

        spread_csp_fractions = models_table$spread_csp[TF1 == TRUE]
        spread_csp_fractions_node_ages = models_table$node_age[TF1 == TRUE]

        lm_df1_spread_CSP = linear_regression_plot(x=spread_csp_fractions_node_ages, y=spread_csp_fractions)
        title(paste0("start state", group, ", spread_csp, no tips"))

        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/LM_", version, "ss", group, "spread.png"))
        lm_df1_spread_CSP = linear_regression_plot(x=spread_csp_fractions_node_ages, y=spread_csp_fractions)
        title(paste0("start state", group, ", spread_csp, no tips"))
        dev.off()

        coef <- rbind(coef, lm_df1_spread_CSP$coefficients)

        # Binned plots
        # control bins

        # Bins plot *including* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
        ###### Rename things! #######

        ######
        #CONTROL binned_scale
        #######   
        CSP <- models_table$control_csp[TF1 == TRUE]
        node_age <- models_table$node_age[TF1 == TRUE]
        Control<- data.frame(CSP, node_age) 

        Control.CSP_NO_Tip = control_csp_fractions 
        Control.NODE_AGE = control_csp_fractions_node_ages
        Spread.CSP_NO_Tip = spread_csp_fractions
        Spread.NODE_AGE = spread_csp_fractions_node_ages 

        breaks <- c(0:20,30)
        tags <- c("[0-1)","[1-2)","[2-3)","[3-4)","[4-5)","[5-6)","[6-7)","[7-8)","[8-9)","[9-10)","[10-11)","[11-12)","[12-13)","[13-14)","[14-15)","[15-16)","[16-17)","[17-18)","[18-19)","[19-20)","[20+)" )
        #tags <- c("[0-1)","[1-2)","[2-3)","[3-4)","[4-5)","[5-6)","[6-7)","[7-8)","[8-9)","[9-10)","[10-11)","[11-12)","[12-13)","[13-14)","[14-15)","[15-16)","[16-17)","[17-18)","[18-19)","[19-20)","[20-21)","[21-22)","[22-23)","[23-24)","[24-25)","[25-26)","[26-27)","[27-28)","[28-29)","[29-30)" )

        group_tags <- cut(Control$node_age, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
        summary(group_tags)


        v <- Control %>% select(node_age, CSP) #pick the variable 
        vgroup <- as_tibble(v) %>% 
            mutate(tag = case_when(
                node_age < 1 ~ tags[1],
                node_age >= 1 & node_age < 2 ~ tags[2],
                node_age >= 2 & node_age < 3 ~ tags[3],
                node_age >= 3 & node_age < 4 ~ tags[4],
                node_age >= 4 & node_age < 5 ~ tags[5],
                node_age >= 5 & node_age < 6 ~ tags[6],
                node_age >= 6 & node_age < 7 ~ tags[7],
                node_age >= 7 & node_age < 8 ~ tags[8],
                node_age >= 8 & node_age < 9 ~ tags[9],
                node_age >= 9 & node_age < 10 ~ tags[10],
                node_age >= 10 & node_age < 11 ~ tags[11],
                node_age >= 11 & node_age < 12 ~ tags[12],
                node_age >= 12 & node_age < 13 ~ tags[13],
                node_age >= 13 & node_age < 14 ~ tags[14],
                node_age >= 14 & node_age < 15 ~ tags[15],
                node_age >= 15 & node_age < 16 ~ tags[16],
                node_age >= 16 & node_age < 17 ~ tags[17],
                node_age >= 17 & node_age < 18 ~ tags[18],
                node_age >= 18 & node_age < 19 ~ tags[19],
                node_age >= 19 & node_age < 20 ~ tags[20],
                node_age > 20 ~ tags[21]
                )
            )
        summary(vgroup)

        vgroup$tag <- factor(vgroup$tag,
                            levels = tags,
                            ordered = FALSE)
        summary(vgroup$tag)

        control_plot <- ggplot() + 
            geom_point(data = vgroup %>% 
            # Group the data by brand then get means
            group_by(tag) %>% 
            summarise(mean_csp = mean(CSP)), mapping = aes(y = mean_csp, x = tag),  size = 5, color = 'red') +
            labs(x='Node Age (mya)', y = 'Correst State Probability (CSP)',  title = paste0("Start State ", group, ', Control Model, No Tips')) +
            ylim(0,1) +
            theme_minimal()

        control_plot

        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/", version, "SS", group, "_controlbinned.png"), control_plot)
        vgroup1 <- vgroup

        ######
        # SPREAD BINED NO TIPS
        #####

        CSP <- models_table$spread_csp[TF1 == TRUE]
        node_age <- models_table$node_age[TF1 == TRUE]

        Spread<- data.frame(CSP, node_age) 

        breaks <- c(0:20,30)
        tags <- c("[0-1)","[1-2)","[2-3)","[3-4)","[4-5)","[5-6)","[6-7)","[7-8)","[8-9)","[9-10)","[10-11)","[11-12)","[12-13)","[13-14)","[14-15)","[15-16)","[16-17)","[17-18)","[18-19)","[19-20)","[20+)" )
        #tags <- c("[0-1)","[1-2)","[2-3)","[3-4)","[4-5)","[5-6)","[6-7)","[7-8)","[8-9)","[9-10)","[10-11)","[11-12)","[12-13)","[13-14)","[14-15)","[15-16)","[16-17)","[17-18)","[18-19)","[19-20)","[20-21)","[21-22)","[22-23)","[23-24)","[24-25)","[25-26)","[26-27)","[27-28)","[28-29)","[29-30)" )

        group_tags <- cut(Spread$node_age, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)
        summary(group_tags)


        v <- Spread %>% select(node_age, CSP) #pick the variable 
        vgroup2 <- as_tibble(v) %>% 
            mutate(tag = case_when(
                node_age < 1 ~ tags[1],
                node_age >= 1 & node_age < 2 ~ tags[2],
                node_age >= 2 & node_age < 3 ~ tags[3],
                node_age >= 3 & node_age < 4 ~ tags[4],
                node_age >= 4 & node_age < 5 ~ tags[5],
                node_age >= 5 & node_age < 6 ~ tags[6],
                node_age >= 6 & node_age < 7 ~ tags[7],
                node_age >= 7 & node_age < 8 ~ tags[8],
                node_age >= 8 & node_age < 9 ~ tags[9],
                node_age >= 9 & node_age < 10 ~ tags[10],
                node_age >= 10 & node_age < 11 ~ tags[11],
                node_age >= 11 & node_age < 12 ~ tags[12],
                node_age >= 12 & node_age < 13 ~ tags[13],
                node_age >= 13 & node_age < 14 ~ tags[14],
                node_age >= 14 & node_age < 15 ~ tags[15],
                node_age >= 15 & node_age < 16 ~ tags[16],
                node_age >= 16 & node_age < 17 ~ tags[17],
                node_age >= 17 & node_age < 18 ~ tags[18],
                node_age >= 18 & node_age < 19 ~ tags[19],
                node_age >= 19 & node_age < 20 ~ tags[20],
                node_age > 20 ~ tags[21]
                )
            )
        summary(vgroup2)

        vgroup2$tag <- factor(vgroup2$tag,
                            levels = tags,
                            ordered = FALSE)
        summary(vgroup2$tag)

        spread_plot <- ggplot() + 
            geom_point(data = vgroup %>% 
            # Group the data by brand then get means
            group_by(tag) %>% 
            summarise(mean_csp = mean(CSP)), mapping = aes(y = mean_csp, x = tag),  size = 5, color = "#00BFC4") +
            labs(x='Node Age (mya)', y = 'Correst State Probability (CSP)',  title = paste0("Start State ", group, ', Spread Model, No Tips')) +
            guides(color=FALSE) + ylim(0,1) +
            theme_minimal()  

        spread_plot

        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/", version, "ss", group, "_Spreadbinned.png"), spread_plot)


        # ALL TOGETHER NOW


        vgroup1$Model <- "Control"
        vgroup2$Model <- "Spread"
        
        combo_table <- rbind(vgroup1, vgroup2)
        final_table <- combo_table %>% 
            # Group the data by brand then get means
            group_by(Model,tag) %>% 
            summarise(mean_csp = mean(CSP))
        
        final_table$Model = as.factor(final_table$Model)
        summary(final_table)

        combo_plot <- ggplot(final_table, aes(y=mean_csp, x=tag, color=Model)) + 
            geom_point(size = 5) +
            labs(x='Node Age (mya)', y = 'Correst State Probability (CSP)',  title = paste0("Start State ", group, ' No Tips')) +
            ylim(0,1) +
            scale_color_manual(values = c(Control = "#F8766D", Spread = "#00BFC4")) +
            theme_minimal()  

        combo_plot

        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/", version, "ss", group, "_Allbinned2.png"), combo_plot)



        # OLD BINS

        # Bins plot *excluding* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/", version, "ss", group, "control_binned_old.png"))
        binsreg(y=Control.CSP_NO_Tip, x=Control.NODE_AGE, line=NULL, plotyrange=c(0,1))
        dev.off()

        # with 95% CI
        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/", version, "ss", group, "control_binned90_old.png"))
        binsreg(y=Control.CSP_NO_Tip, x=Control.NODE_AGE, line=TRUE, ci=TRUE, plotyrange=c(0,1))
        dev.off()



        # 4. spread_csp
        # Bins plot *excluding* 0s and 1s ("plotyrange" doesn't appear to work to set y-axis, boo!)
        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/", version, "ss", group, "spread_binned_old.png"))
        binsreg(y=Spread.CSP_NO_Tip, x=Spread.NODE_AGE, line=NULL, plotyrange=c(0,1))
        dev.off()

        # with 95% CI
        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/", version, "ss", group, "spread_90_old.png"))
        binsreg(y=Spread.CSP_NO_Tip, x=Spread.NODE_AGE, line=TRUE, ci=TRUE, plotyrange=c(0,1))
        dev.off()




        # 3. plot_table$CSP

        TF2 = plot_table$node_age > 0.00000001

        CSP_fractions = plot_table$CSP[TF2 == TRUE]
        CSP_fractions_node_ages = plot_table$node_age[TF2 == TRUE]
        CSP_fractions_model = plot_table$Model[TF2 == TRUE]

        df_big = data.frame(CSP_fractions, CSP_fractions_node_ages, CSP_fractions_model) 

        png(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/", version, "ss", group, "_Biglm.png"))
        lm_plot_table_CSP = linear_regression_plot(x=CSP_fractions_node_ages, y=CSP_fractions, xlab="Node Age", ylab="Correct State Probability (CSP)")
        title(paste0("Start State", group, "All points, No Tips"))
        dev.off()

        png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/", version, "ss", group, "_Control.png"))
        lm_df1_control_CSP = linear_regression_plot(x=control_csp_fractions_node_ages, y=control_csp_fractions, xlab="Node Age", ylab="Correct State Probability (CSP)")
        title(paste0("Start State", group, ", Control Model, No Tips"))
        dev.off()

        png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/", version, "ss", group, "_Spread.png"))
        lm_df1_spread_CSP = linear_regression_plot(x=spread_csp_fractions_node_ages, y=spread_csp_fractions, xlab="Node Age", ylab="Correct State Probability (CSP)")
        title(paste0("Start State", group, ", Spread Model, No Tips"))
        dev.off()

        names = c("Control", "Spread")

        plot <-ggplot(df_big, aes(CSP_fractions_node_ages, CSP_fractions, color = CSP_fractions_model)) +
                geom_point() + geom_smooth(method = 'lm') +
                labs(x = "Node Age (mya)", y = "Correct State Probability (CSP)", title = paste0("Start State ", group)) +
                xlim(0,25) +
                scale_color_manual(values = c(cont = "#F8766D", spread = "#00BFC4"), labels = names, name = "Model")


        plot

        ggsave(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/NoTips/", version, "ss", group, "_gglm.png"), plot)
    }
#}

coef

#########
# lnl histogram!
#########


model = "spread_only"
version = "V1"
big = NULL
big2 = NULL


for (group in 2:8){
    wd2 <- paste0("/Users/wbla447/Desktop/School/Wallis_Castor/September2023/", model , "/" , version, "/")
    setwd(wd2)

    sim_folder <- list.files()
    sim_folder

    plottable_cont <- NULL

    control_table = read.csv(paste0("bigrootstateprobs", group, "_control_7000.csv"))
    spread_table = read.csv(paste0("bigrootstateprobs", group, "inf_7000.csv"))

    head(control_table)

    control_table$totallnl
    spread_table$totallnl

    test = spread_table$totallnl - control_table$totallnl
    test

    png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lnL/", version, "ss", group, "_TotalLnL_title.png"))
    hist(test, main=paste0("Histogram of Spread lnL - Control lnL \n Across 600 tests for Start State", group), xlab = "Net LogLikelihood")        
    dev.off()

    png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lnL/", version, "ss", group, "_TotalLnL_notitle.png"))
    hist(test, main=paste0("Start State", group), xlab = "Net LogLikelihood")        
    dev.off()

    big = append(big,test)

    control_table$rootlnl
    spread_table$rootlnl

    test2 = spread_table$rootlnl - control_table$rootlnl
    test2

    png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lnL/", version, "ss", group, "_RootLnL_title.png"))
    hist(test2, main=paste0("Histogram of Spread root lnL - Control root lnL \n Across 600 tests for Start State", group), xlab = "Net LogLikelihood")        
    dev.off()

    png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lnL/", version, "ss", group, "_RootLnL_notitle.png"))
    hist(test2, main=paste0("Start State", group), xlab = "Net Root LogLikelihood")        
    dev.off()

    big2 = append(big2,test2)


}

summary(big)

above_sig = sum(big > 1.92)
above_sig / length(big)

neg_sig = sum(big < -1.92)
neg_sig / length(big)

post_only = subset(big, big > 0)
post_only

png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lnL/OneBigHistogram_title.png"))
hist(big, main = "Histogram of Spread lnL - Control lnL \n Across all Start States", xlab = "Net LogLikelihood")
dev.off()

png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lnL/OneBigHistogram_notitle.png"))
hist(big, xlab = "Net LogLikelihood")
dev.off()

summary(big2)

png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lnL/OneBigHistogram_Root_title.png"))
hist(big, main = "Histogram of Spread Root lnL - Control Root lnL \n Across all Start States", xlab = "Net Root LogLikelihood")
dev.off()

png(file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/lnL/OneBigHistogram_Root_notitle.png"))
hist(big, xlab = "Net Root LogLikelihood")
dev.off()







#####
# TREE
#####

model = "spread_only"
version = "V1"
startstate = 5
simulation = "050"

dir.create(file.path(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/trees"), paste0("sim", simulation)), showWarnings = FALSE)
#setwd(file.path(paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/trees"), paste0("sim", simulation)))

for(i in 1){
    wd2 <- paste0("/Users/wbla447/Desktop/School/Wallis_Castor/September2023/", model , "/" , version, "/ss", startstate, "/sim_", simulation)
    setwd(wd2)

    # bring tree in
    tr = read.tree("living_tree.newick")
    tr

    # choose colors
    colors = c(1 == "grey", 2 == "red", 3 == "blue", 4 == "yellow", 5 == "purple", 6 == "orange", 7 == "green", 8 == "black")
    colors2 = c("grey", "red", "blue", "yellow", "purple", "orange", "green", "black")
    colors3 = c("black", "#f65356", "#5daad8", "#feea80", "#c56bba", "#fea671", "#a5de94", "grey")
    statenums = c(1:8)




    # bring original states in
        simstates <- read.table("simstates_fixed.txt", header = TRUE)
        colnames(simstates) = "x"

    # get node labels
    node_text = simstates$x[51:99]
    node_node = c(51:99)

    tip_text = simstates$x[1:50]
    tip_tip = c(1:50)

    png(file = (file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/trees/sim", simulation ,"/original_states.png")), height = 850, width = 850)
    plot(tr, label.offset = 0.25)
    nodelabels(text = node_text, node = node_node, bg = colors3[node_text], cex = 2)
    tiplabels(text = tip_text, tip = tip_tip, bg = colors3[tip_text])
    axisPhylo()
    mtext("millions of years ago", side = 1, line = 2)
    mtext("Simulated States")
    dev.off()



    # bring control inferred states in
    test = read.csv("con_branchtops_7000.csv")
    colnames(test)

    # separate out the state proportions
    control_states = test[1:8]
    head(control_states)
    colnames(control_states) = c(1:8)

    #nodes for pies
    node_text = control_states[51:99,1:8]
    node_node = c(51:99)
    head(node_text)

    tip_text = simstates$x[1:50]
    tip_tip = c(1:50)

    png(file = (file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/trees/sim", simulation ,"/control_states.png")), height = 850, width = 850)
    plot(tr, label.offset = 1.25)
    nodelabels(pie = node_text, node = node_node, piecol = colors3)
    tiplabels(text = tip_text, tip = tip_tip, bg = colors3[tip_text])
    axisPhylo()
    mtext("millions of years ago", side = 1, line = 2)
    mtext("Control Model Inferred States")
    dev.off()


    # Test Inference States!
    test = read.csv("branchtops_inf_7000.csv")
    colnames(test)

    # separate out the state proportions
    inf_states = test[1:8]
    head(inf_states)
    colnames(inf_states) = c(1:8)

    #nodes for pies
    node_text = inf_states[51:99,1:8]
    node_node = c(51:99)
    head(node_text)

    tip_text = simstates$x[1:50]
    tip_tip = c(1:50)

    png(file = (file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/trees/sim", simulation ,"/inf_states.png")), height = 850, width = 850)
    plot(tr, label.offset = 1.25)
    nodelabels(pie = node_text, node = node_node, piecol = colors3)
    tiplabels(text = tip_text, tip = tip_tip, bg = colors3[tip_text])
    axisPhylo()
    mtext("millions of years ago", side = 1, line = 2)
    mtext("Spread Model Inferred States")
    dev.off()




    #####
    #Maybe do the whole tings with Letters?
    #####

    tr = read.tree("living_tree.newick")
    tr

    # choose colors
    colors = c(1 == "grey", 2 == "red", 3 == "blue", 4 == "yellow", 5 == "purple", 6 == "orange", 7 == "green", 8 == "black")
    colors2 = c("grey", "red", "blue", "yellow", "purple", "orange", "green", "black")

    colors_julia = c("black", "#f65356", "#5daad8", "#feea80", "#c56bba", "#fea671", "#a5de94", "grey")
    #colors_r = c("black", "#f65356", "#5daad8", "#feea80", "#c56bba", "#a5de94", "#fea671", "grey")

    state_areas_julia = c("Null", "A", "B", "C", "AB", "AC", "BC", "ABC") # for inferred version
    #state_areas_r = c("Null", "A", "B", "C", "AB", "BC", "AC", "ABC") # for simulated fucked version

    # bring original states in
        simstates <- read.table("simstates_fixed.txt", header = TRUE)
        colnames(simstates) = "x"

    # get node labels
    node_text = simstates$x[51:99]
    node_node = c(51:99)
    #node_text_r = state_areas_r[node_text]
    node_text_true = state_areas_julia[node_text]

    tip_text = simstates$x[1:50]
    tip_tip = c(1:50)
    #tip_text_r = state_areas_r[tip_text]
    tip_text_true = state_areas_julia[tip_text]

    png(file = (file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/trees/sim", simulation ,"/letters_original_states.png")), height = 850, width = 850)
    plot(tr, label.offset = 0.25)
    nodelabels(text = node_text_true, node = node_node, bg = colors_julia[node_text], cex = 1.5)
    tiplabels(text = tip_text_true, tip = tip_tip, bg = colors_julia[tip_text])
    axisPhylo()
    mtext("millions of years ago", side = 1, line = 2)
    mtext("Simulated States")
    dev.off()



    # bring control inferred states in
    test = read.csv("con_branchtops_7000.csv")
    colnames(test)

    # separate out the state proportions
    control_states = test[1:8]
    head(control_states)
    colnames(control_states) = state_areas_julia

    #nodes for pies
    node_text = control_states[51:99,1:8]
    node_node = c(51:99)
    head(node_text)

    tip_text = simstates$x[1:50]
    tip_tip = c(1:50)
    tip_text_true = state_areas_julia[tip_text]

    png(file = (file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/trees/sim", simulation ,"/letters_control_states.png")), height = 850, width = 850)
    plot(tr, label.offset = 1.25)
    nodelabels(pie = node_text, node = node_node, piecol = colors_julia)
    tiplabels(text = tip_text_true, tip = tip_tip, bg = colors_julia[tip_text])
    axisPhylo()
    mtext("millions of years ago", side = 1, line = 2)
    mtext("Control Model Inferred States")
    dev.off()


    # Test Inference States!
    test = read.csv("branchtops_inf_7000.csv")
    colnames(test)

    # separate out the state proportions
    inf_states = test[1:8]
    head(inf_states)
    colnames(inf_states) = state_areas_julia

    #nodes for pies
    node_text = inf_states[51:99,1:8]
    node_node = c(51:99)
    head(node_text)

    tip_text = simstates$x[1:50]
    tip_tip = c(1:50)

    png(file = (file = paste0("/Users/wbla447/Desktop/WallisCastor_Git/Plots/", model, "/new_nodes/trees/sim", simulation ,"/letters_inf_states.png")), height = 850, width = 850)
    plot(tr, label.offset = 1.25)
    nodelabels(pie = node_text, node = node_node, piecol = colors_julia)
    tiplabels(text = tip_text_true, tip = tip_tip, bg = colors_julia[tip_text])
    axisPhylo()
    mtext("millions of years ago", side = 1, line = 2)
    mtext("Spread Model Inferred States")
    dev.off()
    }
