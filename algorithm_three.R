################################################################## MODEL SELECTION CODE ################################################################################




#### Loading Libraries ####

library(Bergm)
library(ergm)

library(igraph)
library(statnet)
library(intergraph)
library(tidyverse)




######################## Load categories of Endogenous variables obtained from categorizing endo variables #####################




load("//nas//longleaf//home//helal//endo_categories_lazega_1.RData")



#############################################################################################################################################################




#### Create Network Motifs ####

###  edge ###       ### order 1 ###

adjacency_edge <- matrix( c(0,1,1,0), nrow=2, ncol=2)
graph_edge <- graph.adjacency(adjacency_edge, mode="undirected")



##########  2-stars ##############  #### order 2 ###

adjacency_2star <- matrix( c(0,1,0,1,0,1,0,1,0), nrow=3 , ncol=3)

graph_2star <- graph.adjacency(adjacency_2star, mode="undirected")




#### triangles ###### order 3 ####

triangle <- graph.full(3, directed=FALSE)




################################## Function that counts the isomorphisms between an adjacency matrix and Network Motif ##########


isomorphism_function_undirected <- function(adjcaency_matrix, shape){
  ################## Turn network object into graph #################
  
  
  ## Transform into adjacency matrix ##
  
  
  ## Transform adjacency matrix into igraph object ##
  graph_1 <- graph.adjacency(adjcaency_matrix, mode="undirected")
  
  ## Count the number of isomorphisms between graph object and shape ##
  n1 <- graph.count.subisomorphisms.vf2(graph_1, shape)
  return(n1)
}


######################### UNI-VARIATE ERGM SELECTION based on NETWORK MOTIFS #########################################



#### UNIVARIATE MODEL SELECTION FOR endo variables in D1 or D3 ####




### note input is D1 but it can be D3 or any list of endo variables this function just counts the number of triangles, edges and 2-stars ###

uni_count <- function(obs_network, D1){
  
  results <- c()
  
  total_num_triangle_1 <- c()
  total_num_2stars_1 <- c()
  total_num_edges_1 <- c()
  
  for (i in 1:length(D1)){
    
    print(paste("This is univariate ERGM", i, "with endo variable ", gsub("~","", D1[i])))
    uni_ergm[[i]] <- ergm(as.formula(paste("obs_network~","edges + ", D1[i]))  , estimate = "CD" )
    print(uni_ergm[[i]])
    
    b <- summary(uni_ergm[[i]])
    ### save coefficient estimate for simulation ##
    coef_estimates <- b$coefficients[,1]
    
    
    
    ### STORE THE NUMBER OF NETWORK MOTIFS ###
    
    num_triangle <- c()
    num_2stars <- c()
    num_edges <- c()
    
    
    for(j in 1:10){
      
      simulate_lazega_full_model_1 <- simulate(as.formula(paste("obs_network~","edges + ", D1[i])) , nsim=10,
                                               coef= coef_estimates)
      
      
      A1 <- lapply(simulate_lazega_full_model_1, as.matrix)
      
      
      
      
      
      a1 <- A1[[1]]
      a2 <- A1[[2]]
      a3 <- A1[[3]]
      a4 <- A1[[4]]
      a5 <- A1[[5]]
      a6 <- A1[[6]]
      a7 <- A1[[7]]
      a8 <- A1[[8]]
      a9 <- A1[[9]]
      a10 <- A1[[10]]
      
      
      n_2star_1 <- isomorphism_function_undirected(a1, triangle)
      n_2star_2 <- isomorphism_function_undirected(a2, triangle)
      n_2star_3 <- isomorphism_function_undirected(a3, triangle)
      n_2star_4 <- isomorphism_function_undirected(a4, triangle)
      n_2star_5 <- isomorphism_function_undirected(a5, triangle)
      n_2star_6 <- isomorphism_function_undirected(a6, triangle)
      n_2star_7 <- isomorphism_function_undirected(a7, triangle)
      n_2star_8 <- isomorphism_function_undirected(a8, triangle)
      n_2star_9 <- isomorphism_function_undirected(a9, triangle)
      n_2star_10 <- isomorphism_function_undirected(a10, triangle)
      
      temp_2stars <- cbind(n_2star_1, n_2star_2, n_2star_3, n_2star_4, n_2star_5, n_2star_6, n_2star_7, n_2star_8, n_2star_9, n_2star_10
      )
      
      num_triangle <- cbind(num_triangle, temp_2stars)
      
      
      n_2star_1 <- isomorphism_function_undirected(a1, graph_2star)
      n_2star_2 <- isomorphism_function_undirected(a2, graph_2star)
      n_2star_3 <- isomorphism_function_undirected(a3, graph_2star)
      n_2star_4 <- isomorphism_function_undirected(a4, graph_2star)
      n_2star_5 <- isomorphism_function_undirected(a5, graph_2star)
      n_2star_6 <- isomorphism_function_undirected(a6, graph_2star)
      n_2star_7 <- isomorphism_function_undirected(a7, graph_2star)
      n_2star_8 <- isomorphism_function_undirected(a8, graph_2star)
      n_2star_9 <- isomorphism_function_undirected(a9, graph_2star)
      n_2star_10 <- isomorphism_function_undirected(a10, graph_2star)
      
      temp_2stars <- cbind(n_2star_1, n_2star_2, n_2star_3, n_2star_4, n_2star_5, n_2star_6, n_2star_7, n_2star_8, n_2star_9, n_2star_10
      )
      
      ######### STORE 2-stars ##########
      
      
      num_2stars <- cbind(num_2stars, temp_2stars)
      
      
      n_2star_1 <- isomorphism_function_undirected(a1, graph_edge)
      n_2star_2 <- isomorphism_function_undirected(a2, graph_edge)
      n_2star_3 <- isomorphism_function_undirected(a3, graph_edge)
      n_2star_4 <- isomorphism_function_undirected(a4, graph_edge)
      n_2star_5 <- isomorphism_function_undirected(a5, graph_edge)
      n_2star_6 <- isomorphism_function_undirected(a6, graph_edge)
      n_2star_7 <- isomorphism_function_undirected(a7, graph_edge)
      n_2star_8 <- isomorphism_function_undirected(a8, graph_edge)
      n_2star_9 <- isomorphism_function_undirected(a9, graph_edge)
      n_2star_10 <- isomorphism_function_undirected(a10, graph_edge)
      
      temp_2stars <- cbind(n_2star_1, n_2star_2, n_2star_3, n_2star_4, n_2star_5, n_2star_6, n_2star_7, n_2star_8, n_2star_9, n_2star_10
      )
      
      num_edges <- cbind(num_edges, temp_2stars)
    }
    
    total_num_triangle_1 <- cbind(total_num_triangle_1, num_triangle)
    total_num_2stars_1 <- cbind(total_num_2stars_1, num_2stars)
    total_num_edges_1 <- cbind(total_num_edges_1, num_edges)
    
    final_network_motif <- rbind(total_num_edges_1, total_num_2stars_1, total_num_triangle_1)
    
    
  }
  results <- rbind(results, final_network_motif)
  return(results)
}






#### creates matrix with 3 columns ###






####### Bi Variate ERGM network motif COUNT #############

### expand.grid gives us the combinations ###

bivariate_count <- function(obs_network, endo_list_1, endo_list_2){
  
  ### create list of endogenous variables PAIRS ##
  
  t1 <- as.data.frame(expand.grid(endo_list_1, endo_list_2))
  
  
  results <- c()
  
  total_num_triangle_1 <- c()
  total_num_2stars_1 <- c()
  total_num_edges_1 <- c()
  
  
  
  for (i in 1:dim(t1)[1]){
    
    print(paste("This is bivariate ERGM", i, "with endo variable ", gsub("~","", t1[i,])))
    uni_ergm[[i]] <- ergm(as.formula(paste("obs_network~","edges + ", t1[i,1], "+", t1[i,2]))  , estimate = "CD" )
    print(uni_ergm[[i]])
    
    b <- summary(uni_ergm[[i]])
    ### save coefficient estimate for simulation ##
    coef_estimates <- b$coefficients[,1]
    
    
    
    ### STORE THE NUMBER OF NETWORK MOTIFS ###
    
    num_triangle <- c()
    num_2stars <- c()
    num_edges <- c()
    
    
    for(j in 1:10){
      
      simulate_lazega_full_model_1 <- simulate(as.formula(paste("obs_network~","edges + ", t1[i,1], "+", t1[i,2])) , nsim=10,
                                               coef= coef_estimates)
      
      
      A1 <- lapply(simulate_lazega_full_model_1, as.matrix)
      
      
      
      
      
      a1 <- A1[[1]]
      a2 <- A1[[2]]
      a3 <- A1[[3]]
      a4 <- A1[[4]]
      a5 <- A1[[5]]
      a6 <- A1[[6]]
      a7 <- A1[[7]]
      a8 <- A1[[8]]
      a9 <- A1[[9]]
      a10 <- A1[[10]]
      
      
      n_2star_1 <- isomorphism_function_undirected(a1, triangle)
      n_2star_2 <- isomorphism_function_undirected(a2, triangle)
      n_2star_3 <- isomorphism_function_undirected(a3, triangle)
      n_2star_4 <- isomorphism_function_undirected(a4, triangle)
      n_2star_5 <- isomorphism_function_undirected(a5, triangle)
      n_2star_6 <- isomorphism_function_undirected(a6, triangle)
      n_2star_7 <- isomorphism_function_undirected(a7, triangle)
      n_2star_8 <- isomorphism_function_undirected(a8, triangle)
      n_2star_9 <- isomorphism_function_undirected(a9, triangle)
      n_2star_10 <- isomorphism_function_undirected(a10, triangle)
      
      temp_2stars <- cbind(n_2star_1, n_2star_2, n_2star_3, n_2star_4, n_2star_5, n_2star_6, n_2star_7, n_2star_8, n_2star_9, n_2star_10
      )
      
      num_triangle <- cbind(num_triangle, temp_2stars)
      
      
      n_2star_1 <- isomorphism_function_undirected(a1, graph_2star)
      n_2star_2 <- isomorphism_function_undirected(a2, graph_2star)
      n_2star_3 <- isomorphism_function_undirected(a3, graph_2star)
      n_2star_4 <- isomorphism_function_undirected(a4, graph_2star)
      n_2star_5 <- isomorphism_function_undirected(a5, graph_2star)
      n_2star_6 <- isomorphism_function_undirected(a6, graph_2star)
      n_2star_7 <- isomorphism_function_undirected(a7, graph_2star)
      n_2star_8 <- isomorphism_function_undirected(a8, graph_2star)
      n_2star_9 <- isomorphism_function_undirected(a9, graph_2star)
      n_2star_10 <- isomorphism_function_undirected(a10, graph_2star)
      
      temp_2stars <- cbind(n_2star_1, n_2star_2, n_2star_3, n_2star_4, n_2star_5, n_2star_6, n_2star_7, n_2star_8, n_2star_9, n_2star_10
      )
      
      ######### STORE 2-stars ##########
      
      
      num_2stars <- cbind(num_2stars, temp_2stars)
      
      
      n_2star_1 <- isomorphism_function_undirected(a1, graph_edge)
      n_2star_2 <- isomorphism_function_undirected(a2, graph_edge)
      n_2star_3 <- isomorphism_function_undirected(a3, graph_edge)
      n_2star_4 <- isomorphism_function_undirected(a4, graph_edge)
      n_2star_5 <- isomorphism_function_undirected(a5, graph_edge)
      n_2star_6 <- isomorphism_function_undirected(a6, graph_edge)
      n_2star_7 <- isomorphism_function_undirected(a7, graph_edge)
      n_2star_8 <- isomorphism_function_undirected(a8, graph_edge)
      n_2star_9 <- isomorphism_function_undirected(a9, graph_edge)
      n_2star_10 <- isomorphism_function_undirected(a10, graph_edge)
      
      temp_2stars <- cbind(n_2star_1, n_2star_2, n_2star_3, n_2star_4, n_2star_5, n_2star_6, n_2star_7, n_2star_8, n_2star_9, n_2star_10
      )
      
      num_edges <- cbind(num_edges, temp_2stars)
    }
    
    total_num_triangle_1 <- cbind(total_num_triangle_1, num_triangle)
    total_num_2stars_1 <- cbind(total_num_2stars_1, num_2stars)
    total_num_edges_1 <- cbind(total_num_edges_1, num_edges)
    
    final_network_motif <- rbind(total_num_edges_1, total_num_2stars_1, total_num_triangle_1)
    
    
  }
  results <- rbind(results, final_network_motif)
  return(results)
}


D1_1 <- cbind(D1_1, D3_1)


uni_results_1 <- uni_count(network_1, D1_1)  # 900 rows --> 100 rows for each ERGM
uni_results_2 <- uni_count(network_1, D1_1)
uni_results_3 <- uni_count(network_1, D1_1)
uni_results_4 <- uni_count(network_1, D1_1)
uni_results_5 <- uni_count(network_1, D1_1)
uni_results_6 <- uni_count(network_1, D1_1)
uni_results_7 <- uni_count(network_1, D1_1)
uni_results_8 <- uni_count(network_1, D1_1)
uni_results_9 <- uni_count(network_1, D1_1)
uni_results_10 <- uni_count(network_1, D1_1)



bi_results_1 <- bivariate_count(network_1, D1_1, D1_1) ## 8100 rows --> 100 rows for each ERGM
bi_results_2 <- bivariate_count(network_1, D1_1, D1_1)
bi_results_3 <- bivariate_count(network_1, D1_1, D1_1)
bi_results_4 <- bivariate_count(network_1, D1_1, D1_1)
bi_results_5 <- bivariate_count(network_1, D1_1, D1_1)
bi_results_6 <- bivariate_count(network_1, D1_1, D1_1)
bi_results_7 <- bivariate_count(network_1, D1_1, D1_1)
bi_results_8 <- bivariate_count(network_1, D1_1, D1_1)
bi_results_9 <- bivariate_count(network_1, D1_1, D1_1)
bi_results_10 <- bivariate_count(network_1, D1_1, D1_1)


### For a given model, we will have 10 draws of parameter estimates and for each parameter estimate 
# we will simulate the model 100 times. So in total every ERGM should have 1000 observations for number of edges,
# triangles and 2-stars ... 



save.image("lazega_model_selection.RData")

