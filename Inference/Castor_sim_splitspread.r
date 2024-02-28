

# Remember:
# 5 - AB
# 6 - AC
# 7 - BC
# 8 - ABC

# Note that D_rate is 0.2 in this!

library(castor)
library(ape)
library(BioGeoBEARS)
library(zeallot)
library(tools)
library(cladoRcpp)


Times_3 = rep(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),3)
Land1_3 = rep(c("test1","test1","test2"),each=11)
Land2_3 = rep(c("test2","test3","test3"),each=11)
Distances_km_3 = c(0, 1000, 2000, 3500, 4100, 4900, 6000, 5800, 7000, 7500, 8000,0, 2000, 3000, 4000, 5100, 5900, 7000, 6800, 8000, 8500, 9000,0, 500, 1000, 2500, 3100, 3900, 5000, 4800, 6000, 6500, 7000) 
testdf = data.frame(Times_3, Distances_km_3, Land1_3, Land2_3)
geographic_iso = 7000

getwd()
setwd("/Users/wbla447/Desktop/School/WallisCastor/April_2023")
#    Times_3 Distances_km_3 Land1_3 Land2_3
# 1        0              0   test1   test2
# 2        1           1000   test1   test2
# 3        2           2000   test1   test2
# 4        3           3500   test1   test2
# 5        4           4100   test1   test2
# 6        5           4900   test1   test2
# 7        6           6000   test1   test2
# 8        7           5800   test1   test2
# 9        8           7000   test1   test2
# 10       9           7500   test1   test2
# 11      10           8000   test1   test2
# 12       0              0   test1   test3
# 13       1           2000   test1   test3
# 14       2           3000   test1   test3
# 15       3           4000   test1   test3
# 16       4           5100   test1   test3
# 17       5           5900   test1   test3
# 18       6           7000   test1   test3
# 19       7           6800   test1   test3
# 20       8           8000   test1   test3
# 21       9           8500   test1   test3
# 22      10           9000   test1   test3
# 23       0              0   test2   test3
# 24       1            500   test2   test3
# 25       2           1000   test2   test3
# 26       3           2500   test2   test3
# 27       4           3100   test2   test3
# 28       5           3900   test2   test3
# 29       6           5000   test2   test3
# 30       7           4800   test2   test3
# 31       8           6000   test2   test3
# 32       9           6500   test2   test3
# 33      10           7000   test2   test3


 
 
landnames = c("A" = "test1","B" = "test2","C" = "test3")
statenames = c("1" = "E", "2" = "A", "3" = "B", "4" = "C", "5" = "AB", "6" = "BC", "7" = "AC", "8" = "ABC")
Nlands = 3
Nstates = 8

time_grid = seq(0,10,0.1) 



state8table_julia = read.csv("/Users/wbla447/Desktop/School/WallisCastor/April_2023/state8.csv") # Read in Julia style table

#    event i j k pair wt      prob       rate        val
# 1      y 2 2 2    1  1 1.0000000 0.32881640 0.32881640
# 2      y 3 3 3    1  1 1.0000000 0.32881640 0.32881640
# 3      y 4 4 4    1  1 1.0000000 0.32881640 0.32881640
# 4      v 5 3 2    2  2 0.3333333 0.10960547 0.10960547
# 5      s 5 5 2    2  2 0.3333333 0.10960547 0.10960547
# 6      s 5 5 3    2  2 0.3333333 0.10960547 0.10960547
# 7      v 6 4 2    2  2 0.3333333 0.10960547 0.10960547
# 8      s 6 6 2    2  2 0.3333333 0.10960547 0.10960547
# 9      s 6 6 4    2  2 0.3333333 0.10960547 0.10960547
# 10     v 7 4 3    2  2 0.3333333 0.10960547 0.10960547
# 11     s 7 7 3    2  2 0.3333333 0.10960547 0.10960547
# 12     s 7 7 4    2  2 0.3333333 0.10960547 0.10960547
# 13     v 8 5 4    2  2 0.1666667 0.05480273 0.05480273
# 14     v 8 6 3    2  2 0.1666667 0.05480273 0.05480273
# 15     v 8 7 2    2  2 0.1666667 0.05480273 0.05480273
# 16     s 8 8 2    2  2 0.1666667 0.05480273 0.05480273
# 17     s 8 8 3    2  2 0.1666667 0.05480273 0.05480273
# 18     s 8 8 4    2  2 0.1666667 0.05480273 0.05480273

# Below is the version we make by hand? (7 state version)

#        [,1] [,2] [,3]      [,4]  Column
# tmprow    1    1    1 1.0000000  1
# tmprow    2    2    2 1.0000000  2
# tmprow    3    3    3 1.0000000  3
# tmprow    4    4    4 1.0000000  4
# tmprow    5    2    3 0.3333333  5
# tmprow    5    2    5 0.3333333  6
# tmprow    5    3    5 0.3333333  7
# tmprow    5    5    5 0.0000000  8
# tmprow    6    3    4 0.3333333  9
# tmprow    6    3    6 0.3333333  10
# tmprow    6    4    6 0.3333333  11
# tmprow    6    6    6 0.0000000  12
# tmprow    7    2    4 0.3333333  13
# tmprow    7    2    7 0.3333333  14
# tmprow    7    4    7 0.3333333  15
# tmprow    7    7    7 0.0000000  16

# In order to just use the julia ones we will need:

state8table = state8table_julia[,c(2:4,7)] #take only needed rows
state8table = rbind(c(1,1,1,1),state8table) # insert null states
state8table

#    i j k      prob
# 1  1 1 1 1.0000000
# 2  2 2 2 1.0000000
# 3  3 3 3 1.0000000
# 4  4 4 4 1.0000000
# 5  5 3 2 0.3333333
# 6  5 5 2 0.3333333
# 7  5 5 3 0.3333333
# 8  6 4 2 0.3333333
# 9  6 6 2 0.3333333
# 10 6 6 4 0.3333333
# 11 7 4 3 0.3333333
# 12 7 7 3 0.3333333
# 13 7 7 4 0.3333333
# 14 8 5 4 0.1666667
# 15 8 6 3 0.1666667
# 16 8 7 2 0.1666667
# 17 8 8 2 0.1666667
# 18 8 8 3 0.1666667
# 19 8 8 4 0.1666667

state8table_hand = NULL # compare to how we previously made the table
tmprow = c(1,1,1,1.0)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(2,2,2,1.0)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(3,3,3,1.0)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(4,4,4,1.0)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(5,2,3,1/3)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(5,2,5,1/3)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(5,3,5,1/3)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(5,5,5,0.0)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(6,3,4,1/3)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(6,3,6,1/3)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(6,4,6,1/3)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(6,6,6,0.0)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(7,2,4,1/3)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(7,2,7,1/3)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(7,4,7,1/3)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(7,7,7,0.0)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(8,5,4,1/6)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(8,6,2,1/6)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(8,7,3,1/6)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(8,8,5,1/6)
state8table_hand= rbind(state8table_hand,tmprow)
tmprow = c(8,8,6,1/6)
state8table_hand = rbind(state8table_hand,tmprow)
tmprow = c(8,8,7,1/6)
state8table_hand = rbind(state8table_hand,tmprow)
tmprow = c(8,8,8,0.0)
state8table_hand = rbind(state8table_hand,tmprow)
state8table_hand

#        [,1] [,2] [,3]      [,4]
# tmprow    1    1    1 1.0000000
# tmprow    2    2    2 1.0000000
# tmprow    3    3    3 1.0000000
# tmprow    4    4    4 1.0000000
# tmprow    5    2    3 0.3333333 V
# tmprow    5    2    5 0.3333333
# tmprow    5    3    5 0.3333333
# tmprow    5    5    5 0.0000000
# tmprow    6    3    4 0.3333333 V
# tmprow    6    3    6 0.3333333
# tmprow    6    4    6 0.3333333
# tmprow    6    6    6 0.0000000
# tmprow    7    2    4 0.3333333 V
# tmprow    7    2    7 0.3333333
# tmprow    7    4    7 0.3333333
# tmprow    7    7    7 0.0000000
# tmprow    8    5    4 0.1666667 V
# tmprow    8    6    2 0.1666667 V
# tmprow    8    7    3 0.1666667 V
# tmprow    8    8    5 0.1666667
# tmprow    8    8    6 0.1666667
# tmprow    8    8    7 0.1666667
# tmprow    8    8    8 0.0000000

# NOTE: Differences between Julia Table and Handdone Table:
# Julia table does not have AB -> AB + AB
# We must inset Null -> Null + Null (Done above)
# The order of AB, BC, and AC have changed! (just keep this in mine. I think this matters with the way castor ran the simulation?)


# WE WILL USE HAND TABLE!

A = get_random_mk_transition_matrix(Nstates=8, rate_model="ER", max_rate=0.1)
A[,] = 0.0

# Make it more like a DEC model (anagenetic)
d_rate = 0.3 # range expansion
e_rate = 0.3 # range contraction

# 8 states are:
# null, A, B, C, AB, BC, AC, ABC
A[2,5] = d_rate # A->AB
A[3,5] = d_rate # B->AB
A[3,6] = d_rate # B->BC
A[4,6] = d_rate # C->BC
A[2,7] = d_rate # A->AC
A[4,7] = d_rate # C->AC
A[5,8] = d_rate # AB->ABC
A[6,8] = d_rate # BC->ABC
A[7,8] = d_rate # AC->ABC

A[2,1] = e_rate # A->null
A[3,1] = e_rate # B->null
A[4,1] = e_rate # C->null

A[5,2] = e_rate # AB->A
A[5,3] = e_rate # AB->B
A[6,3] = e_rate # BC->B
A[6,4] = e_rate # BC->C
A[7,2] = e_rate # AC->A
A[7,4] = e_rate # AC->C
A[8,5] = e_rate # ABC->AB
A[8,6] = e_rate # ABC->BC
A[8,7] = e_rate # ABC->AC
A

# In Q transition matrices, the diagonals = -sum(off-diagonal for that row)
diag(A) = 0.0
A
diag(A) = -rowSums(A)
A

Amatrix_array_const = array(data=0.0, dim=c(nrow(A), ncol(A), length(time_grid)))

#       [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]
# [1,] 0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000
# [2,] 0.028 -0.096  0.000  0.000  0.034  0.000  0.034  0.000
# [3,] 0.028  0.000 -0.096  0.000  0.034  0.034  0.000  0.000
# [4,] 0.028  0.000  0.000 -0.096  0.000  0.034  0.034  0.000
# [5,] 0.000  0.028  0.028  0.000 -0.090  0.000  0.000  0.034
# [6,] 0.000  0.000  0.028  0.028  0.000 -0.090  0.000  0.034
# [7,] 0.000  0.028  0.000  0.028  0.000  0.000 -0.090  0.034
# [8,] 0.000  0.000  0.000  0.000  0.028  0.028  0.028 -0.084

state8_matrix = get_random_mk_transition_matrix(Nstates=8, rate_model="ER", max_rate=0.1)

state8_matrix[,]=0.0 # make entire matrix 0.0

# fill matrix
state8_matrix[1,1] = 1.0  # null -> null
state8_matrix[2,2] = 1.0  # A -> A
state8_matrix[3,3] = 1.0  # B -> B
state8_matrix[4,4] = 1.0  # C -> C

# fill collumn 5 AB

state8_matrix[5,1] = 0.0
state8_matrix[5,2] = 2/6 # AB -> A,AB or AB,A
state8_matrix[5,3] = 2/6 # AB -> B,AB or AB,B
state8_matrix[5,5] = 2/6 # AB -> A,B or B,A??? I think that's what this one is? check in again.

# fill collumn 6 BC

state8_matrix[6,1] = 0.0
state8_matrix[6,3] = 2/6 # BC -> B,BC or BC,B
state8_matrix[6,4] = 2/6 # BC -> C,BC or BC,C
state8_matrix[6,6] = 2/6 # BC -> C,B or B,C??? I think that's what this one is? check in again.

# fill collumn 7 AC

state8_matrix[7,1] = 0.0
state8_matrix[7,2] = 2/6 # AC -> A,AC or A,C
state8_matrix[7,4] = 2/6 # AC -> C,AC or C,A
state8_matrix[7,7] = 2/6 # AC -> AC,C or AC,A

# fill collum 8 ABC

state8_matrix[8,1] = 0.0
state8_matrix[8,2] = 1/12 # ABC -> A,BC
state8_matrix[8,3] = 1/12 # ABC -> B,AC
state8_matrix[8,4] = 1/12 # ABC -> C,AB
state8_matrix[8,5] = 2/12 # ABC -> AB,C or AB,ABC
state8_matrix[8,6] = 2/12 # ABC -> BC,A or BC,ABC
state8_matrix[8,7] = 2/12 # ABC -> AC,B or AC,ABC
state8_matrix[8,8] = 3/12 # ABC -> ABC,AB or ABC,BC or ABC,AC 

state8_matrix

#      [,1]       [,2]       [,3]       [,4]      [,5]      [,6]      [,7] [,8]
# [1,]    1 0.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.00
# [2,]    0 1.00000000 0.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.00
# [3,]    0 0.00000000 1.00000000 0.00000000 0.0000000 0.0000000 0.0000000 0.00
# [4,]    0 0.00000000 0.00000000 1.00000000 0.0000000 0.0000000 0.0000000 0.00
# [5,]    0 0.33333333 0.33333333 0.00000000 0.3333333 0.0000000 0.0000000 0.00
# [6,]    0 0.00000000 0.33333333 0.33333333 0.0000000 0.3333333 0.0000000 0.00
# [7,]    0 0.33333333 0.00000000 0.33333333 0.0000000 0.0000000 0.3333333 0.00
# [8,]    0 0.08333333 0.08333333 0.08333333 0.1666667 0.1666667 0.1666667 0.25

state8_table_indc =  state8table_hand[,1:3] - 1
state8_table_probs = state8table_hand[,4]
state8_table_indc_R = state8table_hand[,1:3]

state8_trans_matrix_array_const = array(data=0.0, dim=c(nrow(state8_matrix), ncol(state8_matrix), length(time_grid)))




for(i in 1:length(time_grid)){
    state8_trans_matrix_array_const[,,i] = state8_matrix
}

distance_interp_big <- function(df, time, land1, land2){
    correct_lands = NULL

    for (i in 1:nrow(df)){
        if (df$Distances_km_3[i] != "NA" && (df$Land1_3[i] == land1 && df$Land2_3[i] == land2 || df$Land1_3[i] == land2 && df$Land2_3[i] == land1)){
            a = as.numeric(df$Times_3[i])
            b = as.numeric(df$Distances_km_3[i])
            correct_lands = rbind(correct_lands, c(a,b))
        }
    }

    x = 1
    for (timecheck in correct_lands[,1]){
        if (time - timecheck <= 0){
            if(time == 0){
                time_low = correct_lands[x,1]
                time_high = correct_lands[x,1]
                i = x
                break
            }else{
            time_low = correct_lands[x-1,1]
            time_high = correct_lands[x,1]
            i = x
            break
            }
        }
        x = x + 1
    }
    

    if(time==0){
        dist_high = correct_lands[i,2]
        dist_low = correct_lands[i,2]
    }else{
    dist_high = correct_lands[i,2]
    dist_low = correct_lands[i - 1,2]
    }
    
    t2b_high = abs(time - time_high)
    t2b_low = abs(time - time_low)
    t2bsum = t2b_high + t2b_low
    
    weight_high = (t2bsum - t2b_high)/t2bsum
    weight_low = (t2bsum - t2b_low)/t2bsum
    
    weight_dist_high = weight_high * dist_high
    weight_dist_low = weight_low * dist_low
    
    distance =  weight_dist_high + weight_dist_low

    if(distance == "NaN"){
        distance = 0
    }
    
    distance
}

split_interp_mult = function(df, time_list, land1, land2, geographic_iso, time = NA){ # Time = NA is for if you are asking for a specific time
    
    dist_list = NULL
    for (i in 1:length(time_list)){
        dist_list = append(dist_list, distance_interp_big(df,time_list[i],land1,land2)) # Next Steps: In order to check all lands within one function, may need a matrix?
    }

    scale_high = geographic_iso
    scale_low = 0
    all_split_rates = NULL

    for (i in 1:length(dist_list)){
        if(dist_list[i] > geographic_iso){all_split_rates = rbind(all_split_rates,1) }
        else{
        all_split_rates = rbind(all_split_rates,(dist_list[i])/scale_high)}
    }

    if (is.na(time)){
        all_split_rates
    } else {
        t = which(time_list == time)
        single_split_rate = all_split_rates[t]
        single_split_rate
    }

    #split_rates = NULL
    #for(i in 1:length(base_rates)){
    #    split_rates = split_rates
    #}
    
}

spread_interp_mult = function(df,time_list, land1, land2, geographic_iso, time = NA){

    dist_list = NULL

    for (i in 1:length(time_list)){
        dist_list = append(dist_list, distance_interp_big(df,time_list[i],land1,land2))
    }

    scale_high = geographic_iso
    scale_low = 4000
    spread_rates = NULL

    for (i in 1:length(dist_list)){
        if (dist_list[i] > geographic_iso){spread_rates = rbind(spread_rates,0)}
        else if(dist_list[i] <= scale_low){spread_rates = rbind(spread_rates,d_rate)}
        else{
            spread_rates = rbind(spread_rates,((1-((dist_list[i]-scale_low)/(scale_high-scale_low))))*d_rate)}
    }

    if (is.na(time)){
        spread_rates
    } else {
        t = which(time_list == time)
        single_spread_rate = spread_rates[t]
        single_spread_rate
    }

}


getlands = function(landnames, statenames, daughter){

    daughter_let = statenames[daughter] #get the letter state of the daughters

    n=1
    daughter_lands = substring(daughter_let, seq(1, nchar(daughter_let), n), seq(n, nchar(daughter_let), n)) # split the letters in the land individually

    return(landnames[daughter_lands])
}


comparingsplits = function(distancedf, time_list, state_indc, landnames, statenames, geographic_iso, row, time = NA){

    all_split_rates = NULL 

 # if the daughters are the same as the parent, skip

    daughter1 = state_indc[row,2] 
    daughter2 = state_indc[row,3]

    d1landnames = getlands(landnames,statenames,daughter1)
    d2landnames = getlands(landnames,statenames,daughter2)


    for(i in time_list){
        split_rate = 1

        for(land1 in d1landnames){
            for(land2 in d2landnames){ # compare both sets to each other
                if (land1 == land2){
                    split_rate = 0
                    break
                }

                ini_rate = split_interp_mult(distancedf, time_list, land1, land2, geographic_iso, i)
                if (ini_rate <= split_rate){split_rate = ini_rate} #if the new rate is lower then it becomes the new split_rate
            }
        if (split_rate == 0) {break} 
        }

        all_split_rates = append(all_split_rates, split_rate) #create that list for each row?
    }
    
    if (is.na(time)){
        return(all_split_rates)
    } else {
        t = which(time_list == time)
        single_split_rate = all_split_rates[t]
        return(single_split_rate)
    }

} 

compare_spread = function(df,time_list, state1, state2, newstate, geographic_iso, time = NA){ # This is for hard coding

    spread_rates1 = spread_interp_mult(df, time_list, state1, newstate, geographic_iso) # rates for first letter to new letter
    spread_rates2 = spread_interp_mult(df, time_list, state2, newstate, geographic_iso) # rates for second land to new land
    spread_rates = NULL

    for (i in 1:length(spread_rates1)){
        if(spread_rates1[i] >= spread_rates2[i]){ bigger_rate = spread_rates1[i]} #if A to C is bigger or equal to B to C, take A to C (becasue that rate is bettter!)
        else{bigger_rate = spread_rates2[i]}
        spread_rates = append(spread_rates, bigger_rate)
    }  #take whichever rate is bigger, put them all in a list

    return(spread_rates)
}

# create counts for each type at the table

# numberofevent_table = NULL
# for (s in 1:Nstates){  #for each state type
#     n_ofrows = 0
#     n_ofs = 0
#     n_ofv = 0
#     for (r in 1:nrow(state8_table_indc_R)){
#         if (state8_table_indc_R[r,1] == s){ #if the row matches what you state we are on
#             if(state_indc[r,1] == state_indc[r,2] && state_indc[r,1] == state_indc[r,3]){
#                 next
#             }else if(state8_table_indc_R[r,1] == state8_table_indc_R[r,2] || state8_table_indc_R[r,1] == state8_table_indc_R[r,3]){
#                 n_ofs = n_ofs + 1 # +1 subset
#                 n_ofrows = n_ofrows + 1 #if the parent matches the state type we count how many we have
#             }else{
#                 n_ofv = n_ofv + 1
#                 n_ofrows = n_ofrows + 1
#             }
#         }
#     }

#     tmprow = c(s, n_ofrows, n_ofv, n_ofs)
#     numberofevent_table = rbind(numberofevent_table, tmprow)
# }
# numberofevent_table
# #        [,1] [,2] [,3] [,4]
# # tmprow    1    0    0    0
# # tmprow    2    0    0    0
# # tmprow    3    0    0    0
# # tmprow    4    0    0    0
# # tmprow    5    3    1    2
# # tmprow    6    3    1    2
# # tmprow    7    3    1    2
# # tmprow    8    6    3    3


fillsplitmatrix = function(distancedf, time_list, state_matrix, statetable, landnames, statenames, geographic_iso, Nstates){
    
    state_indc =  statetable[,1:3]
    state_probs = statetable[,4]

    timed_probarray = array(data=0.0, dim=c(length(state_probs), 2, length(time_list)))

    for(i in 1:length(time_list)){
        timed_probarray[,1,i] = state_probs
        timed_probarray[,2,i] = state_probs
    }

    numberofevent_table = NULL # we need the number of events to tell us what we need to divide by
    for (s in 1:Nstates){  #for each state type
        n_ofrows = 0
        n_ofs = 0
        n_ofv = 0
        for (r in 1:nrow(state_indc)){
            if (state_indc[r,1] == s){ #if the row matches what you state we are on
                if(state_indc[r,1] == state_indc[r,2] && state_indc[r,1] == state_indc[r,3]){
                    next
                }else if(state_indc[r,1] == state_indc[r,2] || state_indc[r,1] == state_indc[r,3]){
                    n_ofs = n_ofs + 1 # +1 subset
                    n_ofrows = n_ofrows + 1 #if the parent matches the state type we count how many we have
                }else{
                    n_ofv = n_ofv + 1
                    n_ofrows = n_ofrows + 1
                }
            }
        }

        tmprow = c(s, n_ofrows, n_ofv, n_ofs)
        numberofevent_table = rbind(numberofevent_table, tmprow)
    }

    for (i in 1:length(time_list)){
        for (s in 1:Nstates){  #for each state type
            sumrate = 0

            for (r in 1:nrow(state_indc)){ # run through it once to get the sum of all the rates?
                if (state_indc[r,1] == s){ #if the row matches what you state we are on

                    if(state_indc[r,1] == state_indc[r,2] && state_indc[r,1] == state_indc[r,3]){
                        next
                    }else if(state_indc[r,1] == state_indc[r,2] || state_indc[r,1] == state_indc[r,3]){ # this might be where we say if a land is too 'small' then it couldn't handle two species in this area?
                        rate = 1/(numberofevent_table[s,2]) # atm model will use 1/the number of events, instead of the inverse like before (inverse okay with only one vicariance event)
                    }else{
                        rate = comparingsplits(distancedf, time_list, state_indc, landnames, statenames, geographic_iso, r, time_list[i])
                        timed_probarray[r,1,i] = rate
                        timed_probarray[r,2,i] = rate
                    }

                    sumrate = sumrate + rate
                }
            }
            sumrate

            for (r in 1:nrow(state_indc)){ # run through it a second time to weight those rates
                if (state_indc[r,1] == s){ #if the row matches what you state we are on
                    if(state_indc[r,1] == state_indc[r,2] && state_indc[r,1] == state_indc[r,3]){
                        next
                    }else if(state_indc[r,1] == state_indc[r,2] || state_indc[r,1] == state_indc[r,3]){
                        timed_probarray[r,1,i] = (1/(numberofevent_table[s,2]))/sumrate
                        timed_probarray[r,2,i] = (1/(numberofevent_table[s,2]))/sumrate
                    }else{
                        rate = comparingsplits(distancedf, time_list, state_indc, landnames, statenames, geographic_iso, r, time_list[i])
                        timed_probarray[r,1,i] = timed_probarray[r,1,i]/sumrate
                        timed_probarray[r,2,i] = timed_probarray[r,2,i]/sumrate
                    }
                }            
            }        

        }
    }

    return(timed_probarray)

}



Amatrix_array_const
A
spread_ratesAB = spread_interp_mult(testdf, time_grid, "test1", "test2", geographic_iso) #anything to AB
spread_ratesBC = spread_interp_mult(testdf, time_grid, "test2", "test3", geographic_iso)
spread_ratesAC = spread_interp_mult(testdf, time_grid, "test1", "test3", geographic_iso)
spread_ratesAB_C = compare_spread(testdf, time_grid, "test1", "test2", "test3", geographic_iso) # AB to ABC
spread_ratesBC_A = compare_spread(testdf, time_grid, "test3", "test2", "test1", geographic_iso) # BC to ABC
spread_ratesAC_B = compare_spread(testdf, time_grid, "test1", "test3", "test2", geographic_iso) # AC to ABC 



for(i in 1:length(time_grid)){
    Amatrix_array_const[2,5,i] = spread_ratesAB[i]
    Amatrix_array_const[3,5,i] = spread_ratesAB[i] 
    Amatrix_array_const[3,6,i] = spread_ratesBC[i] 
    Amatrix_array_const[4,6,i] = spread_ratesBC[i] 
    Amatrix_array_const[2,7,i] = spread_ratesAC[i]
    Amatrix_array_const[4,7,i] = spread_ratesAC[i]
    Amatrix_array_const[5,8,i] = spread_ratesAB_C[i]
    Amatrix_array_const[6,8,i] = spread_ratesBC_A[i]
    Amatrix_array_const[7,8,i] = spread_ratesAC_B[i]

    diag(Amatrix_array_const[,,i]) = 0.0
    diag(Amatrix_array_const[,,i]) = -rowSums(Amatrix_array_const[,,i])
}

Amatrix_array_const

#state8_trans_table_probs_matrix_const = fillsplitmatrix(testdf, time_grid, state8_matrix, state8table_hand, landnames, statenames, geographic_iso, 8)
#save(state8_trans_table_probs_matrix_const, file = "/Users/wbla447/Desktop/School/WallisCastor/September2023/8statetimedistprobs.RData")
load("/Users/wbla447/Desktop/School/WallisCastor/September2023/8statetimedistprobs.RData")
state8_trans_table_probs_matrix_const

# THIS TAKES 1 HOUR AND 7 MINUTES
# for (i in 1:length(time_grid)){
#     for (s in 1:Nstates){  #for each state type
#         sumrate = 0

#         for (r in 1:nrow(state8_table_indc_R)){ # run through it once to get the sum of all the rates?
#             if (state8_table_indc_R[r,1] == s){ #if the row matches what you state we are on

#                 if(state8_table_indc_R[r,1] == state8_table_indc_R[r,2] && state8_table_indc_R[r,1] == state8_table_indc_R[r,3]){
#                     next
#                 }else if(state8_table_indc_R[r,1] == state8_table_indc_R[r,2] || state8_table_indc_R[r,1] == state8_table_indc_R[r,3]){ # this might be where we say if a land is too 'small' then it couldn't handle two species in this area?
#                     rate = 1/(numberofevent_table[s,2])
#                 }else{
#                     rate = comparingsplits(testdf, time_grid, state8_table_indc_R, landnames, statenames, geographic_iso, time_grid[i])
#                     state8_trans_table_probs_matrix_const[r,1,i] = rate
#                     state8_trans_table_probs_matrix_const[r,2,i] = rate
#                 }

#                 sumrate = sumrate + rate
#             }
#         }
#         sumrate

#         for (r in 1:nrow(state8_table_indc_R)){ # run through it a second time to weight those rates
#             if (state8_table_indc_R[r,1] == s){ #if the row matches what you state we are on
#                 if(state8_table_indc_R[r,1] == state8_table_indc_R[r,2] && state8_table_indc_R[r,1] == state8_table_indc_R[r,3]){
#                     next
#                 }else if(state8_table_indc_R[r,1] == state8_table_indc_R[r,2] || state8_table_indc_R[r,1] == state8_table_indc_R[r,3]){
#                     state8_trans_table_probs_matrix_const[r,1,i] = (1/(numberofevent_table[s,2]))/sumrate
#                     state8_trans_table_probs_matrix_const[r,2,i] = (1/(numberofevent_table[s,2]))/sumrate
#                 }else{
#                     rate = comparingsplits(testdf, time_grid, state8_table_indc_R, landnames, statenames, geographic_iso, time_grid[i])
#                     state8_trans_table_probs_matrix_const[r,1,i] = state8_trans_table_probs_matrix_const[r,1,i]/sumrate
#                     state8_trans_table_probs_matrix_const[r,2,i] = state8_trans_table_probs_matrix_const[r,2,i]/sumrate
#                 }
#             }            
#         }        

#     }
# }


parameters = list(birth_rates                = 0.3,
                  death_rates                = 0.05,
                  transition_matrix_A        = Amatrix_array_const,
                  transition_matrix_C        = state8_trans_matrix_array_const,
                  transition_table_indices_C = state8_table_indc,
                  transition_table_probs_C   = state8_trans_table_probs_matrix_const 
                  )
                  
# simulation = simulate_tdsse2( Nstates        = 8, 
#                              parameters      = parameters, 
#                              start_state	 = 2,
#                              max_tips        = 50, 
#                              time_grid       = time_grid
#                              )

# plot(simulation$tree); axisPhylo(); mtext(text="Millions of years ago (Ma)", side=1, line=2) # PLOT ISN'T WORKING ATM?
# simulation$Ntransitions_A
# simulation$Ntransitions_C

# oldtree = simulation$tree
# oldtable = prt(simulation$tree)
# oldstring = write.tree(simulation$tree, file="")
# oldstring
# oldtabletip = prt(simulation$tree, get_tipnames=TRUE)

# newtree = read.tree(file="", text=oldstring)
# newtree
# newtable = prt(newtree)
# newtabletip = prt(newtree, get_tipnames=TRUE)

# # NOTE THAT THIS NEW TABLE HAS A DIFFERENT NODE ORDER THAN SIMULATION TREE TABLE


# tablematches = match(x=newtabletip$tipnames, table=oldtabletip$tipnames)
# oldtabletip$tipnames[tablematches] == newtabletip$tipnames

# oldstates = c(as.numeric(simulation$tip_states), simulation$node_states)
# newstates = oldstates[tablematches]

# statenames = c("1" = "E", "2" = "A", "3" = "B", "4" = "C", "5" = "AB", "6" = "BC", "7" = "AC", "8" = "ABC")
# statelabels = NULL 

# for (i in 1:length(newstates)){
#     numericstate = newstates[i]
#         statelabels = rbind(statelabels, statenames[numericstate])
# }

# check =  data.frame(newstates,statelabels)
# check
# newtreenolabel = newtree

# tipnums = 1:length(newtree$tip.label)
# newtreenolabel$tip.label = statelabels[tipnums]

# plot(newtreenolabel); axisPhylo(); mtext(text="Millions of years ago (Ma)", side=1, line=2)

# nodenumstart = length(newtree$tip.label) + 1
# nodenumend = length(newtree$tip.label) + newtree$Nnode
# nodenums = nodenumstart:nodenumend
# nodelabels(text=statelabels[nodenums])


# setwd("/Users/wbla447/Desktop/School/WallisCastor/8StateSpread8_095")

source("/Users/wbla447/Desktop/School/WallisCastor/castor_helpers.R")


ss = 5 #start state

v = "3" #version

for (ss in 2:8){
    failures = 0
    for (i in 1:600){
        simulation = simulate_tdsse2(Nstates        = 8, 
                                     parameters      = parameters, 
                                     start_state	 = ss,
                                     max_tips        = 51, 
                                     time_grid       = time_grid
        )
        
        
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }
        
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }
        
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }
        
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }    
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }
        if(simulation$success == FALSE){
            failures = failures + 1
            simulation = simulate_tdsse2( Nstates        = 8, 
                                          parameters      = parameters, 
                                          start_state	 = ss,
                                          max_tips        = 51, 
                                          time_grid       = time_grid
            )
        }
    
        wd = paste0("/Users/wbla447/Desktop/School/WallisCastor/September2023/splitspread/V", v,"/ss", ss, "/")
        setwd(wd)
        dir.create(paste0("sim_", sprintf("%03s", i)))
    
        wd = paste0("/Users/wbla447/Desktop/School/WallisCastor/September2023/splitspread/V", v, "/ss", ss, "/sim_", sprintf("%03s", i),"/")
        setwd(wd)
    
        write_out_original_castor_simfiles(simulation, wd, fn="rawsim")
        simulation2 = remove_last_tip_from_simulation(simulation)
        plot(simulation2$tree)
    
        write_out_original_castor_simfiles(simulation2, wd, fn="rawsim_wo_tip51")
    
        simulation3 = reorder_castor_sim_to_default_ape_node_order(simulation2)
    
        # External inputs for making geography files etc.
        area_names = c("A","B","C")
        states_list = rcpp_areas_list_to_states_list(areas=area_names, maxareas=length(area_names), include_null_range=TRUE)
    
        # Reorder, and also write out the reordered files plus geography etc.
        write_out_reordered_castor_simfiles(simulation3, wd, area_names=c("A","B","C"), states_list=states_list)
    }
    
    failures
    print(failures)
}

# setwd("/Users/wbla447/Desktop/School/WallisCastor/April_2023/ss2")
# testing = list.dirs()
# testing

# for (i in 1:length(testing)){
#     file.rename(testing[i], paste0("sim_", sprintf("%03s", i -1)))
# }


# Option one is the rate of v is 0 when the distance is zero
# Option two is the rate of v is the base rate and the rate increases when the distance increases



