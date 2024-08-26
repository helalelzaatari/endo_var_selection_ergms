# Endogenous Variable Selection Code for ERGMs 

# Algorithm 1 #

########################### BOUNDING INPUT PARAMETER for UNDIRECTED Networks #######################################################################

#### Loading Libraries ####

library(Bergm)
library(ergm)

library(igraph)
library(statnet)
library(intergraph)
library(tidyverse)



########################################################################################################
##################### LIST OF ENDOGENOUS VARIABLES #####################################################




##########  list of endogenous variables that do not take on a natural nubmer #############

endo_list_other <- c("degreepopularity", "degcrossprod", "degcor", "gwesp", "gwnsp", "gwdsp", "triangle",
                     "gwdegree" ,"isolates", "sociality" )


######## list of k-stars up until 10 ########

### KSTAR ###
kstar_list_ <- c("kstar(2)", "kstar(3)", "kstar(4)", "kstar(5)", "kstar(6)", "kstar(7)", "kstar(8)", "kstar(9)", "kstar(10)", "kstar(11)", "kstar(12)","kstar(13)","kstar(14)", "kstar(15)", "kstar(16)", "kstar(17)", "kstar(18)", "kstar(19)", "kstar(20)", "kstar(21)", "kstar(22)", "kstar(23)", "kstar(24)", "kstar(25)", "kstar(26)", "kstar(27)", "kstar(28)","kstar(29)", "kstar(30)", "kstar(31)", "kstar(32)", "kstar(33)", "kstar(34)", "kstar(35)", "kstar(36)", "kstar(37)", "kstar(38)", "kstar(39)", "kstar(40)", "kstar(41)", "kstar(42)", "kstar(43)", "kstar(44)", "kstar(45)", "kstar(46)","kstar(47)","kstar(48)", "kstar(49)", "kstar(50)", "kstar(51)", "kstar(52)", "kstar(53)", "kstar(54)", "kstar(55)", "kstar(56)", "kstar(57)", "kstar(58)", "kstar(59)", "kstar(60)", "kstar(61)", "kstar(62)","kstar(63)", "kstar(64)", "kstar(65)", "kstar(66)", "kstar(67)", "kstar(68)", "kstar(69)", "kstar(70)", "kstar(71)", "kstar(72)", "kstar(73)", "kstar(74)", "kstar(75)", "kstar(76)", "kstar(77)", "kstar(78)","kstar(79)", "kstar(80)", "kstar(81)", "kstar(82)", "kstar(83)", "kstar(84)", "kstar(85)", "kstar(86)", "kstar(87)", "kstar(88)","kstar(89)", "kstar(90)", "kstar(91)", "kstar(92)", "kstar(93)", "kstar(94)", "kstar(95)", "kstar(96)", "kstar(97)", "kstar(98)", "kstar(99)","kstar(100)", "kstar(101)", "kstar(102)", "kstar(103)", "kstar(104)", "kstar(105)", "kstar(106)","kstar(107)","kstar(108)", "kstar(109)", "kstar(110)","kstar(111)", "kstar(112)", "kstar(113)", "kstar(114)", "kstar(115)", "kstar(116)", "kstar(117)","kstar(118)","kstar(119)", "kstar(120)", "kstar(121)", "kstar(122)", "kstar(123)", "kstar(124)", "kstar(125)","kstar(126)", "kstar(127)", "kstar(128)", "kstar(129)", "kstar(130)", "kstar(131)", "kstar(132)","kstar(133)","kstar(134)", "kstar(135)", "kstar(136)","kstar(137)", "kstar(138)", "kstar(139)", "kstar(140)", "kstar(141)", "kstar(142)", "kstar(143)","kstar(144)","kstar(145)", "kstar(146)", "kstar(147)", "kstar(148)", "kstar(149)", "kstar(150)")



#### CYCLE ####
cycle_list <- c("cycle(3)", "cycle(4)","cycle(5)", "cycle(6)", "cycle(7)", "cycle(8)", "cycle(9)", "cycle(10)")



###### ESP ######
esp_list <- c("esp(2)", "esp(3)", "esp(4)", "esp(5)", "esp(6)", "esp(7)","esp(8)", "esp(9)", "esp(10)",
              "esp(11)", "esp(12)", "esp(13)", "esp(14)", "esp(15)", "esp(16)","esp(17)", "esp(18)", "esp(19)",
              "esp(20)", "esp(21)", "esp(22)", "esp(23)", "esp(24)", "esp(25)","esp(26)", "esp(27)", "esp(28)", "esp(29)", "esp(30)", "esp(31)","esp(32)", "esp(33)", "esp(34)",
              "esp(35)", "esp(36)", "esp(37)", "esp(38)", "esp(39)", "esp(40)","esp(41)", "esp(42)", "esp(43)",
              "esp(44)", "esp(45)", "esp(46)", "esp(47)", "esp(48)", "esp(49)", "esp(50)")

####### DSP #######
dsp_list <- c("dsp(2)", "dsp(3)", "dsp(4)", "dsp(5)", "dsp(6)", "dsp(7)","dsp(8)", "dsp(9)", "dsp(10)",
              "dsp(11)", "dsp(12)", "dsp(13)", "dsp(14)", "dsp(15)", "dsp(16)","dsp(17)", "dsp(18)", "dsp(19)",
              "dsp(20)", "dsp(21)", "dsp(22)", "dsp(23)", "dsp(24)", "dsp(25)","dsp(26)", "dsp(27)", "dsp(28)", "dsp(29)", "dsp(30)", "dsp(31)","dsp(32)", "dsp(33)", "dsp(34)",
              "dsp(35)", "dsp(36)", "dsp(37)", "dsp(38)", "dsp(39)", "dsp(40)","dsp(41)", "dsp(42)", "dsp(43)",
              "dsp(44)", "dsp(45)", "dsp(46)", "dsp(47)", "dsp(48)", "dsp(49)", "dsp(50)")


####### NSP #######
nsp_list <- c("nsp(2)", "nsp(3)", "nsp(4)", "nsp(5)", "nsp(6)", "nsp(7)","nsp(8)", "nsp(9)", "nsp(10)",
              "nsp(11)", "nsp(12)", "nsp(13)", "nsp(14)", "nsp(15)", "nsp(16)","nsp(17)", "nsp(18)", "nsp(19)",
              "nsp(20)", "nsp(21)", "nsp(22)", "nsp(23)", "nsp(24)", "nsp(25)", "nsp(26)", "nsp(27)", "nsp(28)", "nsp(29)", "nsp(30)", "nsp(31)","nsp(32)", "nsp(33)", "nsp(34)",
              "nsp(35)", "nsp(36)", "nsp(37)", "nsp(38)", "nsp(39)", "nsp(40)","nsp(41)", "nsp(42)", "nsp(43)",
              "nsp(44)", "nsp(45)", "nsp(46)", "nsp(47)", "nsp(48)", "nsp(49)", "nsp(50)" )



#######################################################################################################

############################ LIST OF NETWORKS ############









###############################################################################################################



############### INITIAL SCREENING FUNCTION FOR EACH ENDO LIST ##############################################






## FUNCTION THAT FINDS THE UPPER BOUND OF THE INPUT PARAMETER ##################################


### OUTPUTS ###  




upper_bound_input_parameter <- function(obs_network, endo_list){
  
  uni_ergm <- list()
  param_estimates <- c()
  
  for (i in 1:length(endo_list)){
    print(paste("This is univariate ERGM", i, "with endo variable ", gsub("~","", endo_list[i])))
    uni_ergm[[i]] <- ergm(as.formula(paste("obs_network~","edges + ", endo_list[i]))  , estimate = "CD" )
    print(uni_ergm[[i]])
    
    param_estimates <- cbind(param_estimates, uni_ergm[[i]]$coefficients[2])
    
  }
  
  
  return(param_estimates)
  
  
}



############## LOAD NETWORKS ################################

####### load data sets first ##########
data("lazega")
data("kapferer")
data("zach")
data("windsurfers")
data("molecule")
data("faux.mesa.high")
data("ecoli")
data("moodyContactSim")
data("florentine")


network_1 <- lazega
network_2 <- kapferer
network_3 <- kapferer2
network_4 <- zach
network_5 <- windsurfers
network_6 <- molecule
network_7 <- faux.mesa.high
network_8 <- ecoli2
network_9 <- moodyContactSim
network_10 <- flomarriage
network_11 <- flobusiness


#########################################################

### network 1 ###

esp_1 <- replicate(n=5, upper_bound_input_parameter(network_1, esp_list), simplify = TRUE )
dsp_1 <- replicate(n=5, upper_bound_input_parameter(network_1, dsp_list), simplify = TRUE )
nsp_1 <- replicate(n=5, upper_bound_input_parameter(network_1, nsp_list), simplify = TRUE )
#cycle_1 <- replicate(n=5, upper_bound_input_parameter(network_1, cycle_list), simplify = TRUE )
kstar_1 <- replicate(n=5, upper_bound_input_parameter(network_1, kstar_list_), simplify = TRUE )


### network 2 ###

esp_2 <- replicate(n=5, upper_bound_input_parameter(network_2, esp_list), simplify = TRUE )
dsp_2 <- replicate(n=5, upper_bound_input_parameter(network_2, dsp_list), simplify = TRUE )
nsp_2 <- replicate(n=5, upper_bound_input_parameter(network_2, nsp_list), simplify = TRUE )
#cycle_2 <- replicate(n=5, upper_bound_input_parameter(network_2, cycle_list), simplify = TRUE )
kstar_2 <- replicate(n=5, upper_bound_input_parameter(network_2, kstar_list_), simplify = TRUE )


### network 3 ###

esp_3 <- replicate(n=5, upper_bound_input_parameter(network_3, esp_list), simplify = TRUE )
dsp_3 <- replicate(n=5, upper_bound_input_parameter(network_3, dsp_list), simplify = TRUE )
nsp_3 <- replicate(n=5, upper_bound_input_parameter(network_3, nsp_list), simplify = TRUE )
#cycle_3 <- replicate(n=5, upper_bound_input_parameter(network_3, cycle_list), simplify = TRUE )
kstar_3 <- replicate(n=5, upper_bound_input_parameter(network_3, kstar_list_), simplify = TRUE )


### network 4 ###

esp_4 <- replicate(n=5, upper_bound_input_parameter(network_4, esp_list), simplify = TRUE )
dsp_4 <- replicate(n=5, upper_bound_input_parameter(network_4, dsp_list), simplify = TRUE )
nsp_4 <- replicate(n=5, upper_bound_input_parameter(network_4, nsp_list), simplify = TRUE )
#cycle_4 <- replicate(n=5, upper_bound_input_parameter(network_4, cycle_list), simplify = TRUE )
kstar_4 <- replicate(n=5, upper_bound_input_parameter(network_4, kstar_list_), simplify = TRUE )


### network 5 ###
esp_5 <- replicate(n=5, upper_bound_input_parameter(network_5, esp_list), simplify = TRUE )
dsp_5 <- replicate(n=5, upper_bound_input_parameter(network_5, dsp_list), simplify = TRUE )
nsp_5 <- replicate(n=5, upper_bound_input_parameter(network_5, nsp_list), simplify = TRUE )
#cycle_5 <- replicate(n=5, upper_bound_input_parameter(network_5, cycle_list), simplify = TRUE )
kstar_5 <- replicate(n=5, upper_bound_input_parameter(network_5, kstar_list_), simplify = TRUE )



### network 6 ###

esp_6 <- replicate(n=5, upper_bound_input_parameter(network_6, esp_list), simplify = TRUE )
dsp_6 <- replicate(n=5, upper_bound_input_parameter(network_6, dsp_list), simplify = TRUE )
nsp_6 <- replicate(n=5, upper_bound_input_parameter(network_6, nsp_list), simplify = TRUE )
#cycle_6 <- replicate(n=5, upper_bound_input_parameter(network_6, cycle_list), simplify = TRUE )
kstar_6 <- replicate(n=5, upper_bound_input_parameter(network_6, kstar_list_), simplify = TRUE )


### network 7 ###

esp_7 <- replicate(n=5, upper_bound_input_parameter(network_7, esp_list), simplify = TRUE )
dsp_7 <- replicate(n=5, upper_bound_input_parameter(network_7, dsp_list), simplify = TRUE )
nsp_7 <- replicate(n=5, upper_bound_input_parameter(network_7, nsp_list), simplify = TRUE )
#cycle_7 <- replicate(n=5, upper_bound_input_parameter(network_7, cycle_list), simplify = TRUE )
kstar_7 <- replicate(n=5, upper_bound_input_parameter(network_7, kstar_list_), simplify = TRUE )




### network 8 ###

esp_8 <- replicate(n=5, upper_bound_input_parameter(network_8, esp_list), simplify = TRUE )
dsp_8 <- replicate(n=5, upper_bound_input_parameter(network_8, dsp_list), simplify = TRUE )
nsp_8 <- replicate(n=5, upper_bound_input_parameter(network_8, nsp_list), simplify = TRUE )
#cycle_8 <- replicate(n=5, upper_bound_input_parameter(network_8, cycle_list), simplify = TRUE )
kstar_8 <- replicate(n=5, upper_bound_input_parameter(network_8, kstar_list_), simplify = TRUE )




### network 9 ###

esp_9 <- replicate(n=5, upper_bound_input_parameter(network_9, esp_list), simplify = TRUE )
dsp_9 <- replicate(n=5, upper_bound_input_parameter(network_9, dsp_list), simplify = TRUE )
nsp_9 <- replicate(n=5, upper_bound_input_parameter(network_9, nsp_list), simplify = TRUE )
#cycle_9 <- replicate(n=5, upper_bound_input_parameter(network_9, cycle_list), simplify = TRUE )
kstar_9 <- replicate(n=5, upper_bound_input_parameter(network_9, kstar_list_), simplify = TRUE )



## network 10 ##

esp_10 <- replicate(n=5, upper_bound_input_parameter(network_10, esp_list), simplify = TRUE )
dsp_10 <- replicate(n=5, upper_bound_input_parameter(network_10, dsp_list), simplify = TRUE )
nsp_10 <- replicate(n=5, upper_bound_input_parameter(network_10, nsp_list), simplify = TRUE )
#cycle_9 <- replicate(n=5, upper_bound_input_parameter(network_10, cycle_list), simplify = TRUE )
kstar_10 <- replicate(n=5, upper_bound_input_parameter(network_10, kstar_list_), simplify = TRUE )




## network 11 ##

esp_11 <- replicate(n=5, upper_bound_input_parameter(network_11, esp_list), simplify = TRUE )
dsp_11 <- replicate(n=5, upper_bound_input_parameter(network_11, dsp_list), simplify = TRUE )
nsp_11 <- replicate(n=5, upper_bound_input_parameter(network_11, nsp_list), simplify = TRUE )
#cycle_9 <- replicate(n=5, upper_bound_input_parameter(network_10, cycle_list), simplify = TRUE )
kstar_11 <- replicate(n=5, upper_bound_input_parameter(network_11, kstar_list_), simplify = TRUE )



##########################################################################################################
################### FINALIZE LIST OF ENDOGENOUS VARIABLES FOR EACH NETWORK FOR INITIAL SCREENING
###############################################################################################################

############ SELECT LIST OF ENDOGENOUS VARIABLES FOR EACH NETWORK #####################################



######################################## NETWORK 1 #######################################################


esp_list_1 <- c()

for (i in 1:length(esp_list)){
  if (esp_1[i] > -10000){
    esp_list_1 <- cbind(esp_list_1,esp_list[i])
  }
}



dsp_list_1 <- c()

for (i in 1:length(dsp_list)){
  if (dsp_1[i] > -10000){
    dsp_list_1 <- cbind(dsp_list_1,dsp_list[i])
  }
}


nsp_list_1 <- c()
for (i in 1:length(nsp_list)){
  if (nsp_1[i] > -10000){
    nsp_list_1 <- cbind(nsp_list_1,nsp_list[i])
  }
}



#cycle_list_1 <- c()
#for (i in 1:length(cycle_list)){
 # if (cycle_1[i] > -10000){
 #   cycle_list_1 <- cbind(cycle_list_1,cycle_list[i])
#  }
#}


kstar_list_1 <- c()
for (i in 1:length(kstar_list_)){
  if (kstar_1[i] > -10000){
    kstar_list_1 <- cbind(kstar_list_1,kstar_list_[i])
  }
}

all_endo_list_1 <- c(esp_list_1, dsp_list_1, nsp_list_1, kstar_list_1, endo_list_other)


###########                        NETWORK 2                  ######################


esp_list_2 <- c()

for (i in 1:length(esp_list)){
  if (esp_2[i] > -10000){
    esp_list_2 <- cbind(esp_list_2,esp_list[i])
  }
}



dsp_list_2 <- c()

for (i in 1:length(dsp_list)){
  if (dsp_2[i] > -10000){
    dsp_list_2 <- cbind(dsp_list_2,dsp_list[i])
  }
}


nsp_list_2 <- c()
for (i in 1:length(nsp_list)){
  if (nsp_2[i] > -10000){
    nsp_list_2 <- cbind(nsp_list_2, nsp_list[i])
  }
}



kstar_list_2 <- c()
for (i in 1:length(kstar_list_)){
  if (kstar_2[i] > -10000){
    kstar_list_2 <- cbind(kstar_list_2,kstar_list_[i])
  }
}

all_endo_list_2 <- c(esp_list_2, dsp_list_2, nsp_list_2, kstar_list_2, endo_list_other)




###########                        NETWORK 3                  ######################


esp_list_3 <- c()

for (i in 1:length(esp_list)){
  if (esp_3[i] > -10000){
    esp_list_3 <- cbind(esp_list_3,esp_list[i])
  }
}



dsp_list_3 <- c()

for (i in 1:length(dsp_list)){
  if (dsp_3[i] > -10000){
    dsp_list_3 <- cbind(dsp_list_3,dsp_list[i])
  }
}


nsp_list_3 <- c()
for (i in 1:length(nsp_list)){
  if (nsp_3[i] > -10000){
    nsp_list_3 <- cbind(nsp_list_3, nsp_list[i])
  }
}





kstar_list_3 <- c()
for (i in 1:length(kstar_list_)){
  if (kstar_3[i] > -10000){
    kstar_list_3 <- cbind(kstar_list_3,kstar_list_[i])
  }
}

all_endo_list_3 <- c(esp_list_3, dsp_list_3, nsp_list_3, kstar_list_3, endo_list_other)





###########                        NETWORK 4                  ######################


esp_list_4 <- c()

for (i in 1:length(esp_list)){
  if (esp_4[i] > -10000){
    esp_list_4 <- cbind(esp_list_4,esp_list[i])
  }
}



dsp_list_4 <- c()

for (i in 1:length(dsp_list)){
  if (dsp_4[i] > -10000){
    dsp_list_4 <- cbind(dsp_list_4,dsp_list[i])
  }
}


nsp_list_4 <- c()
for (i in 1:length(nsp_list)){
  if (nsp_4[i] > -10000){
    nsp_list_4 <- cbind(nsp_list_4, nsp_list[i])
  }
}





kstar_list_4 <- c()
for (i in 1:length(kstar_list_)){
  if (kstar_4[i] > -10000){
    kstar_list_4 <- cbind(kstar_list_4,kstar_list_[i])
  }
}

all_endo_list_4 <- c(esp_list_4, dsp_list_4, nsp_list_4, kstar_list_4, endo_list_other)







###########                        NETWORK 5                  ######################


esp_list_5 <- c()

for (i in 1:length(esp_list)){
  if (esp_5[i] > -10000){
    esp_list_5 <- cbind(esp_list_5,esp_list[i])
  }
}



dsp_list_5 <- c()

for (i in 1:length(dsp_list)){
  if (dsp_5[i] > -10000){
    dsp_list_5 <- cbind(dsp_list_5,dsp_list[i])
  }
}


nsp_list_5 <- c()
for (i in 1:length(nsp_list)){
  if (nsp_5[i] > -10000){
    nsp_list_5 <- cbind(nsp_list_5, nsp_list[i])
  }
}




kstar_list_5 <- c()
for (i in 1:length(kstar_list_)){
  if (kstar_5[i] > -10000){
    kstar_list_5 <- cbind(kstar_list_5,kstar_list_[i])
  }
}

all_endo_list_5 <- c(esp_list_5, dsp_list_5, nsp_list_5, kstar_list_5, endo_list_other)







###########                        NETWORK 6                  ######################


esp_list_6 <- c()

for (i in 1:length(esp_list)){
  if (esp_6[i] > -10000){
    esp_list_6 <- cbind(esp_list_6,esp_list[i])
  }
}



dsp_list_6 <- c()

for (i in 1:length(dsp_list)){
  if (dsp_6[i] > -10000){
    dsp_list_6 <- cbind(dsp_list_6,dsp_list[i])
  }
}


nsp_list_6 <- c()
for (i in 1:length(nsp_list)){
  if (nsp_6[i] > -10000){
    nsp_list_6 <- cbind(nsp_list_6, nsp_list[i])
  }
}





kstar_list_6 <- c()
for (i in 1:length(kstar_list_)){
  if (kstar_6[i] > -10000){
    kstar_list_6 <- cbind(kstar_list_6,kstar_list_[i])
  }
}

all_endo_list_6 <- c(esp_list_6, dsp_list_6, nsp_list_6, kstar_list_6, endo_list_other)






###########                        NETWORK 7                  ######################



esp_list_7 <- c()

for (i in 1:length(esp_list)){
  if (esp_7[i] > -10000){
    esp_list_7 <- cbind(esp_list_7,esp_list[i])
  }
}



dsp_list_7 <- c()

for (i in 1:length(dsp_list)){
  if (dsp_7[i] > -10000){
    dsp_list_7 <- cbind(dsp_list_7,dsp_list[i])
  }
}


nsp_list_7 <- c()
for (i in 1:length(nsp_list)){
  if (nsp_7[i] > -10000){
    nsp_list_7 <- cbind(nsp_list_7, nsp_list[i])
  }
}



kstar_list_7 <- c()
for (i in 1:length(kstar_list_)){
  if (kstar_7[i] > -10000){
    kstar_list_7 <- cbind(kstar_list_7,kstar_list_[i])
  }
}

all_endo_list_7 <- c(esp_list_7, dsp_list_7, nsp_list_7, kstar_list_7, endo_list_other)






###########                        NETWORK 8                  ######################

esp_list_8 <- c()

for (i in 1:length(esp_list)){
  if (esp_8[i] > -10000){
    esp_list_8 <- cbind(esp_list_8,esp_list[i])
  }
}



dsp_list_8 <- c()

for (i in 1:length(dsp_list)){
  if (dsp_8[i] > -10000){
    dsp_list_8 <- cbind(dsp_list_8,dsp_list[i])
  }
}


nsp_list_8 <- c()
for (i in 1:length(nsp_list)){
  if (nsp_8[i] > -10000){
    nsp_list_8 <- cbind(nsp_list_8, nsp_list[i])
  }
}





kstar_list_8 <- c()
for (i in 1:length(kstar_list_)){
  if (kstar_8[i] > -10000){
    kstar_list_8 <- cbind(kstar_list_8,kstar_list_[i])
  }
}

all_endo_list_8 <- c(esp_list_8, dsp_list_8, nsp_list_8, kstar_list_8, endo_list_other)





###########                        NETWORK 9                  ######################



esp_list_9 <- c()

for (i in 1:length(esp_list)){
  if (esp_9[i] > -10000){
    esp_list_9 <- cbind(esp_list_9,esp_list[i])
  }
}



dsp_list_9 <- c()

for (i in 1:length(dsp_list)){
  if (dsp_9[i] > -10000){
    dsp_list_9 <- cbind(dsp_list_9,dsp_list[i])
  }
}


nsp_list_9 <- c()
for (i in 1:length(nsp_list)){
  if (nsp_9[i] > -10000){
    nsp_list_9 <- cbind(nsp_list_9, nsp_list[i])
  }
}






kstar_list_9 <- c()
for (i in 1:length(kstar_list_)){
  if (kstar_9[i] > -10000){
    kstar_list_9 <- cbind(kstar_list_9,kstar_list_[i])
  }
}

all_endo_list_9 <- c(esp_list_9, dsp_list_9, nsp_list_9, kstar_list_9, endo_list_other)





###########                        NETWORK 10                  ######################



esp_list_10 <- c()

for (i in 1:length(esp_list)){
  if (esp_10[i] > -10000){
    esp_list_10 <- cbind(esp_list_10,esp_list[i])
  }
}



dsp_list_10 <- c()

for (i in 1:length(dsp_list)){
  if (dsp_10[i] > -10000){
    dsp_list_10 <- cbind(dsp_list_10,dsp_list[i])
  }
}


nsp_list_10 <- c()
for (i in 1:length(nsp_list)){
  if (nsp_10[i] > -10000){
    nsp_list_10 <- cbind(nsp_list_10, nsp_list[i])
  }
}






kstar_list_10 <- c()
for (i in 1:length(kstar_list_)){
  if (kstar_10[i] > -10000){
    kstar_list_10 <- cbind(kstar_list_10,kstar_list_[i])
  }
}

all_endo_list_10 <- c(esp_list_10, dsp_list_10, nsp_list_10, kstar_list_10, endo_list_other)




###########                        NETWORK 11                  ######################



esp_list_11 <- c()

for (i in 1:length(esp_list)){
  if (esp_11[i] > -10000){
    esp_list_11 <- cbind(esp_list_11,esp_list[i])
  }
}



dsp_list_11 <- c()

for (i in 1:length(dsp_list)){
  if (dsp_11[i] > -10000){
    dsp_list_11 <- cbind(dsp_list_11,dsp_list[i])
  }
}


nsp_list_11 <- c()
for (i in 1:length(nsp_list)){
  if (nsp_11[i] > -10000){
    nsp_list_11 <- cbind(nsp_list_11, nsp_list[i])
  }
}



kstar_list_11 <- c()
for (i in 1:length(kstar_list_)){
  if (kstar_11[i] > -10000){
    kstar_list_11 <- cbind(kstar_list_11,kstar_list_[i])
  }
}

all_endo_list_11 <- c(esp_list_11, dsp_list_11, nsp_list_11, kstar_list_11, endo_list_other)







save.image("final_endo_lists.RData")

