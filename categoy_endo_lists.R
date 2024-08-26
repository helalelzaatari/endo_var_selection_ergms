##################################################################################################
################### CATEGORIZING THE ENDOGENOUS VARIABLES INTO D1 D2 and D3 #####################################################




#### Loading Libraries #### 

library(Bergm)
library(ergm)

library(igraph)
library(statnet)
library(intergraph)
library(tidyverse)



######################## Load Lists of Endogenous variables obtained from upper bound input parameter #####################

load("C:\\Users\\helal\\Desktop\\final_endo_list.RData")







####################################################################################################



uni_ergm <- list()


relative_AIC_change <- c()




###### THIS FUNCTION IS VERY COMPUTATIONALL INTESNVIE #####

relative_aic <- function(network_data, all_endo_list) {
  
  ################## FIT NULL MODEL ###################
  
  
  null_model <- ergm(network_data~edges, estimate="CD")
  
  null_model <- logLik(null_model, add=TRUE)
  b_0 <- summary(null_model)
  print(b_0)
  
  ### report initial NULL AIC ###
  aic_0 <- b_0$aic
  
  
  
  
  for (i in 1:length(all_endo_list)){
    print(paste("This is univariate ERGM", i, "with endo variable ", gsub("~","", all_endo_list[i])))
    uni_ergm[[i]] <- ergm(as.formula(paste("network_data~","edges + ", all_endo_list[i]))  , estimate = "CD" )
    uni_ergm[[i]] <- logLik(uni_ergm[[i]], add=TRUE)
    
    aic_1 <- summary(uni_ergm[[i]])$aic
    
    if (is.null(aic_1)){
      
      relative_AIC_change <- cbind( relative_AIC_change, "NA")
      
      
    } else {
      b_i <- (as.numeric(aic_0) - as.numeric(aic_1)) / (as.numeric(aic_0))
      
      relative_AIC_change <- cbind(relative_AIC_change, b_i)
    }}
  
  
  return(relative_AIC_change)
  
  
  
  
}


#######################################################


###### DO THIS FOR EACH NETWORK #############
### OUTPUT MATRIX OF RELATIVE AIC CHANGE ########


## relative AIC change matrix for network 1 ##
relative_1  <- replicate(n=5, relative_aic(network_1, all_endo_list_1), simplify="TRUE")


## relative AIC change matrix for network 2 ##
relative_2  <- replicate(n=5, relative_aic(network_2, all_endo_list_2), simplify="TRUE")



## relative AIC change matrix for network 3 ##
relative_3 <- replicate(n=5, relative_aic(network_3, all_endo_list_3), simplify="TRUE")


## relative AIC change matrix for network 4 ##
relative_4  <- replicate(n=5, relative_aic(network_4, all_endo_list_4), simplify="TRUE")


## relative AIC change matrix for network 5 ##
relative_5  <- replicate(n=1, relative_aic(network_5, all_endo_list_5), simplify="TRUE")



## relative AIC change matrix for network 6 ##
relative_6 <- replicate(n=1, relative_aic(network_6, all_endo_list_6), simplify="TRUE")



## relative AIC change matrix for network 1 ##
relative_7  <- replicate(n=1, relative_aic(network_7, all_endo_list_7), simplify="TRUE")


## relative AIC change matrix for network 2 ##
relative_8  <- replicate(n=1, relative_aic(network_8, all_endo_list_8), simplify="TRUE")



## relative AIC change matrix for network 3 ##
relative_9 <- replicate(n=1, relative_aic(network_9, all_endo_list_9), simplify="TRUE")


## relative AIC change matrix for network 4 ##
relative_10  <- replicate(n=1, relative_aic(network_10, all_endo_list_10), simplify="TRUE")


## relative AIC change matrix for network 11 ##
relative_11  <- replicate(n=1, relative_aic(network_11, all_endo_list_11), simplify="TRUE")


















######################################################################################################
############# CLASSIFY ENDOGENOUS VARIABLES INTO DIFFERENT GROUPS BASED ON RELATIVE AIC CHANGE #######




## Create 4 sets of endogenous variables which are D1, D2, D3 and D4 (NAs) ##




D1 <- c()
D2 <- c()
D3 <- c()
D4 <- c()



#### First IDENTIFY THE endogenous variables which yield NAs ####


### input matrix of relative AIC change ###
### input endo list ###

endo_NA_category <- function(r1, all_endo_list){
  
  row_number <- c()
  
  for (i in 1:length(all_endo_list)){
    
    
    temp_1 <- quantile(as.numeric(r1[i,]), probs=c(0.01,0.1,0.2,0.5,0.8), na.rm=TRUE) 
    
    
    
    
    if(is.na(temp_1[5])){
      row_number <- cbind(row_number, i)
    } 
  }
  
  ### Remove the endogenous variables that do not have an AIC ###
  na_removed_endo_list <-  cbind(row_number)
  
  
  return(na_removed_endo_list)
  
  ############### Returns a vector of numbers which is the position of the NAs ########
}




####### STORES VECTOR OF POSITION OF THE NAS FOR EACH NETWORK ######
h1 <- as.numeric(endo_NA_category(relative_1, all_endo_list_1))
h2 <- as.numeric(endo_NA_category(relative_2, all_endo_list_2))
h3 <- as.numeric(endo_NA_category(relative_3, all_endo_list_3))
h4 <- as.numeric(endo_NA_category(relative_4, all_endo_list_4))
h5 <- as.numeric(endo_NA_category(relative_5, all_endo_list_5))
h6 <- as.numeric(endo_NA_category(relative_6, all_endo_list_6))
h7 <- as.numeric(endo_NA_category(relative_7, all_endo_list_7))
h8 <- as.numeric(endo_NA_category(relative_8, all_endo_list_8))
h9 <- as.numeric(endo_NA_category(relative_9, all_endo_list_9))
h10 <- as.numeric(endo_NA_category(relative_10, all_endo_list_10))
h11 <- as.numeric(endo_NA_category(relative_11, all_endo_list_11))









### CREATES SET OF ENDOGENOUS VARIABLES CALLED D4 for each network ###
D4_1 <- all_endo_list_1[h1]
D4_2 <- all_endo_list_2[h2]
D4_3 <- all_endo_list_3[h3]
D4_4 <- all_endo_list_4[h4]
D4_5 <- all_endo_list_1[h5]
D4_6 <- all_endo_list_2[h6]
D4_7 <- all_endo_list_3[h7]
D4_8 <- all_endo_list_4[h8]
D4_9 <- all_endo_list_2[h9]
D4_10 <- all_endo_list_3[h10]
D4_11 <- all_endo_list_4[h11]










#### RETURNS  D1 CATEGORY ####


endo_category_D1 <- function(r1, all_endo_list){
  D1 <- c()
  
  ## identify the NAs ##
  h1 <- as.numeric(endo_NA_category(r1, all_endo_list))
  
  ## remove row of NAs ##
  r1 <- r1[-h1,]
  all_endo_list <- all_endo_list[-h1]
  
  ### number of rows in relative AIC change ###
  
  n <- dim(r1)[1]
  
  for (i in 1:n){
    temp_1 <- quantile(as.numeric(r1[i,]), probs=c(0.01,0.1,0.2,0.5,0.8), na.rm=TRUE)
    if (temp_1[1] > 0){
      D1<- cbind(D1, all_endo_list[i])
    }
  }
  
  return(D1)
}










#### RETURNS  D2 CATEGORY ####


endo_category_D2 <- function(r1, all_endo_list){
  D2 <- c()
  
  
  ## identify the NAs ##
  h1 <- as.numeric(endo_NA_category(r1, all_endo_list))
  
  ## remove row of NAs ##
  r1 <- r1[-h1,]
  all_endo_list <- all_endo_list[-h1]
  
  ### number of rows in relative AIC change ###
  
  n <- dim(r1)[1]
  
  for (i in 1:n){
    temp_1 <- quantile(as.numeric(r1[i,]), probs=c(0.01,0.1,0.2,0.5,0.8), na.rm=TRUE)
    if (temp_1[3] < 0){
      D2<- cbind(D2, all_endo_list[i])
    }
  }
  
  return(D2)
}






####### returns list of D3 ###


endo_category_D3 <- function(r1, all_endo_list){
  D3 <- c()
  
  
  ## identify the NAs ##
  h1 <- as.numeric(endo_NA_category(r1, all_endo_list))
  
  ## remove row of NAs ##
  r1 <- r1[-h1,]
  all_endo_list <- all_endo_list[-h1]
  
  ### number of rows in relative AIC change ###
  
  n <- dim(r1)[1]
  
  for (i in 1:n){
    temp_1 <- quantile(as.numeric(r1[i,]), probs=c(0.01,0.1,0.2,0.5,0.8), na.rm=TRUE)
    if (temp_1[2] > 0 && temp[1] < 0){
      D3<- cbind(D3, all_endo_list[i])
    }
  }
  
  return(D3)
}


###################### ENDO VARIABLES IN D1 for each network ###########################

D1_1 <- endo_category_D1(relative_1, all_endo_list_1)

D1_2 <- endo_category_D1(relative_2, all_endo_list_2)


D1_3 <- endo_category_D1(relative_3, all_endo_list_3)

D1_4 <- endo_category_D1(relative_4, all_endo_list_4)


D1_5 <- endo_category_D1(relative_5, all_endo_list_5)

D1_6 <- endo_category_D1(relative_6, all_endo_list_6)


D1_7 <- endo_category_D1(relative_7, all_endo_list_7)

D1_8 <- endo_category_D1(relative_8, all_endo_list_8)

D1_9 <- endo_category_D1(relative_9, all_endo_list_9)


D1_10 <- endo_category_D1(relative_10, all_endo_list_10)

D1_11 <- endo_category_D1(relative_11, all_endo_list_11)

###################### ENDO VARIABLES IN D2 for each network ###########################


D2_1 <- endo_category_D2(relative_1, all_endo_list_1)

D2_2 <- endo_category_D2(relative_2, all_endo_list_2)


D2_3 <- endo_category_D2(relative_3, all_endo_list_3)

D2_4 <- endo_category_D2(relative_4, all_endo_list_4)


D2_5 <- endo_category_D2(relative_5, all_endo_list_5)

D2_6 <- endo_category_D2(relative_6, all_endo_list_6)


D2_7 <- endo_category_D2(relative_7, all_endo_list_7)

D2_8 <- endo_category_D2(relative_8, all_endo_list_8)

D2_9 <- endo_category_D2(relative_9, all_endo_list_9)


D2_10 <- endo_category_D2(relative_10, all_endo_list_10)

D2_11 <- endo_category_D2(relative_11, all_endo_list_11)



###################### ENDO VARIABLES IN D3 for each network ###########################


D3_1 <- endo_category_D3(relative_1, all_endo_list_1)

D3_2 <- endo_category_D3(relative_2, all_endo_list_2)


D3_3 <- endo_category_D3(relative_3, all_endo_list_3)

D3_4 <- endo_category_D3(relative_4, all_endo_list_4)


D3_5 <- endo_category_D3(relative_5, all_endo_list_5)

D3_6 <- endo_category_D3(relative_6, all_endo_list_6)


D3_7 <- endo_category_D3(relative_7, all_endo_list_7)

D3_8 <- endo_category_D3(relative_8, all_endo_list_8)

D3_9 <- endo_category_D3(relative_9, all_endo_list_9)


D3_10 <- endo_category_D3(relative_10, all_endo_list_10)

D3_11 <- endo_category_D3(relative_11, all_endo_list_11)







save.image("endo_categories.RData")




















