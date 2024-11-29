library(gSeg)
library(NetworkDistance)
library(igraph)
library(ape)
library(ade4)

####### ####### #######
# LOAD DATA
####### ####### #######

# standard deviation of whole data set
data_sd <- sd(as.matrix(data[,3:93]))
data_stand <- data[,3:93]/data_sd

# standard deviation per month
data_sd <- apply(as.matrix(data[,3:93]),2,sd)
data_stand_col <- sweep(as.matrix(data[,3:93]),2,data_sd,FUN='/')

# standard deviation per edge
data_rm <- data[-c(350),]
data_sd <- apply(as.matrix(data_rm[,3:93]),1,sd)
data_stand_col <- sweep(as.matrix(data_rm[,3:93]),1,data_sd,FUN='/')
vertex_set <- union(unique(data_rm[,2]),unique(data_rm[,1]))

# stationary
data_edges=data_rm[,c(1,2)]
data_rm=data_rm[,-c(1,2)]
data_rm <- data.matrix(data_rm)
# difference and normalise
dtn <- stationary_ts(data_rm)
data_norm <- t(dtn[[1]])
data_sd <- dtn[[2]] 
vertex_set <- union(unique(data_edges[,2]),unique(data_edges[,1]))

# stl
library(stats)
stl_list <- apply(data[,-c(1,2)],1,function(x) stl(ts(as.vector(x),start = c(2015,1), end=c(2022,7),frequency = 12),s.window = "periodic"))
stl_mat <- t(sapply(stl_list, function(x) x$time.series[,3]))

vertex_set <- union(unique(data[,2]),unique(data[,1]))
# create list of adjacency matrices
adj_list <- list()
adj_list_binary <- list()
for(i in 1:91){
  #edgelist_loc <- cbind(data[,1:2],data_stand[,i]) # rescaled data with /sd of whole df
  #edgelist_loc <- cbind(data[,1:2],data[,2+i]) # data as recorded
  #edgelist_loc <- cbind(data[,1:2],data_stand_col[,i]) # rescaled data with /sd of cols
  #edgelist_loc <- cbind(data_rm[,1:2],data_stand_col[,i]) # rescaled data with /sd of rows
  #edgelist_loc <- cbind(data_edges,data_norm[,i]) # stationary
  edgelist_loc <- cbind(data[,1:2],stl_mat[,i]) # stl
  edgelist_loc <- edgelist_loc[which(edgelist_loc[,3]!=0),]
  colnames(edgelist_loc) <- c("From","To","Weight")
  graph_loc <- graph.data.frame(edgelist_loc,vertices = vertex_set)
  adj_list[[i]] <- get.adjacency(graph_loc,attr = 'Weight')
  
  # get binary adjacencies
  adj_list_binary[[i]] <- get.adjacency(graph_loc)
}


# Obtaine distance matrices for weighted graphs and calculate statistics
dist_mat_frob <- nd.edd(adj_list,out.dist = TRUE)
dist_mat_spec <- nd.wsd(adj_list,out.dist = TRUE) # Normalized Laplacian matrix contains topological information of a corresponding network via its spectrum. nd.wsd adopts weighted spectral distribution of eigenvalues and brings about a metric via binning strategy.
dist_mat_diff <- nd.gdd(adj_list,out.dist = TRUE) # Graph Diffusion Distance (nd.gdd) quantifies the difference between two weighted graphs of same size. It takes an idea from heat diffusion process on graphs via graph Laplacian exponential kernel matrices.

#mst_adj <- mst(dist_mat$D) # did not work for gseg
#mst_graph <- graph.adjacency(mst_adj,mode = "undirected")
#mst_edeglist <- get.edgelist(mst_graph)
mst_edeglist <- mstree(dist_mat_frob$D)
mst_edeglist <- mstree(dist_mat_diff$D)
gseg1(91,mst_edeglist,statistics = "all")

gseg2(91,mst_edeglist,statistics = "all")


# Binary case
dist_mat_hamm <- nd.hamming(adj_list_binary,out.dist = TRUE)
mst_edeglist <- mstree(dist_mat_hamm$D)
gseg1(91,mst_edeglist,statistics = "all")

################ ################ ################ ################ ################ ################ 
################ FOR COMPERABLE RESULTS DO SAME PREPROCESSING WITH LOGLIKELIHOOD RATIO ###################
################ ################ ################ ################ ################ ################ 

# check for possible change points, which time series are all 0s (before or after) to remove corresponding edges
# remove those edges for all change point cases for comparable results
ts_rem_bef <- ts_rem_aft <- c()
for (c in c(14:81)){ # seq(39,72,3)
  data_bef <- data[,3:c]
  data_aft <- data[,(c+1):93]
  ts_rem_bef <- c(ts_rem_bef,which(apply(data_bef,1,sum)==0))
  ts_rem_aft <- c(ts_rem_aft,which(apply(data_aft,1,sum)==0))
}

edge_ind_rem <- unique(c(ts_rem_bef,ts_rem_aft))

data <- data[-edge_ind_rem,]
data_edges=data[,c(1,2)]
data=data[,-c(1,2)]
data <- data.matrix(data)
vertex_set <- union(unique(data_edges[,2]),unique(data_edges[,1]))
###########################
# differencing before split
###########################

# difference and normalise
dtn <- stationary_ts(data)
data_norm <- t(dtn[[1]])
data_sd <- dtn[[2]] # keep sd of data_train


# create list of adjacency matrices
adj_list <- list()
adj_list_binary <- list()
ind_adj <- 1
for(i in 1:90){#1:90 , c((14-2):(81-2))
  #edgelist_loc <- cbind(data[,1:2],data_stand[,i]) # rescaled data with /sd of whole df
  #edgelist_loc <- cbind(data[,1:2],data[,2+i]) # data as recorded
  #edgelist_loc <- cbind(data[,1:2],data_stand_col[,i]) # rescaled data with /sd of cols
  #edgelist_loc <- cbind(data_rm[,1:2],data_stand_col[,i]) # rescaled data with /sd of rows
  edgelist_loc <- cbind(data_edges,data_norm[,i]) # stationary
  #edgelist_loc <- cbind(data[,1:2],stl_mat[,i]) # stl
  edgelist_loc <- edgelist_loc[which(edgelist_loc[,3]!=0),]
  colnames(edgelist_loc) <- c("From","To","Weight")
  graph_loc <- graph.data.frame(edgelist_loc,vertices = vertex_set)
  adj_list[[ind_adj]] <- get.adjacency(graph_loc,attr = 'Weight')
  
  # get binary adjacencies
  adj_list_binary[[ind_adj]] <- get.adjacency(graph_loc)
  ind_adj <- ind_adj+1
}

# Obtain distance matrices for weighted graphs and calculate statistics
dist_mat_frob <- nd.edd(adj_list,out.dist = TRUE)

mst_edeglist <- mstree(dist_mat_frob$D)
gseg1(90,mst_edeglist,statistics = "all")

###########################
# differencing after split
###########################

# create list of adjacency matrices
adj_list <- list()
adj_list_binary <- list()
ind_adj <- 1
for (c in c((14-2):(81-2))){
  data_bef <- data[,1:c]
  data_aft <- data[,(c+1):91]
  # difference and normalise
  dtn_bef <- stationary_ts(data_bef)
  data_norm_bef <- t(dtn_bef[[1]])
  data_sd_bef <- dtn_bef[[2]] 
  
  dtn_aft <- stationary_ts(data_aft)
  data_norm_aft <- t(dtn_aft[[1]])
  data_sd_aft <- dtn_aft[[2]] 
  
  data_norm <- cbind(data_norm_bef,data_norm_aft)
  edgelist_loc <- cbind(data_edges,data_norm[,c]) # stationary
  
  edgelist_loc <- edgelist_loc[which(edgelist_loc[,3]!=0),]
  colnames(edgelist_loc) <- c("From","To","Weight")
  graph_loc <- graph.data.frame(edgelist_loc,vertices = vertex_set)
  adj_list[[ind_adj]] <- get.adjacency(graph_loc,attr = 'Weight')
  
  # get binary adjacencies
  adj_list_binary[[ind_adj]] <- get.adjacency(graph_loc)
  ind_adj <- ind_adj+1
}

# Obtain distance matrices for weighted graphs and calculate statistics
dist_mat_frob <- nd.edd(adj_list,out.dist = TRUE)

mst_edeglist <- mstree(dist_mat_frob$D)
gseg1(68,mst_edeglist,statistics = "all")
