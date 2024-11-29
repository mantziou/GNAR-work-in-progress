library(igraph)
library(VARDetect) # library for VAR model change point detection
library(NetworkDistance)
library(ade4)
# simulation experiments for change point detection using log-likel ratio and GNAR-edge model

# ER, 20-node, density=.2
# GNAR-edge(2,[2,2])
# T1=100, T2=100

# generate seeds
set.seed(44)
seeds <- sample(1:500,50,replace=FALSE)

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 
# Scenario with change in autoregressive parameters and residual variance:
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

# set parameters
alpha_par_bef <- c(-.1,.2)
beta_par_bef <- list(c(-0.3,-.2),c(0.1,.1))

alpha_par_aft <- c(.01,-.05)
beta_par_aft <- list(c(0.05,.01),c(-0.03,-0.02))

sigma_bef <- 1
sigma_aft <- .5

###### ###### ###### ###### ###### ###### ###### ###### 
# Scenario with change only in autoregressive parameters:
###### ###### ###### ###### ###### ###### ###### ###### 

# set parameters
alpha_par_bef <- c(-.1,.2)
beta_par_bef <- list(c(-0.3,-.2),c(0.1,.1))
#alt
alpha_par_bef <- c(-.3,.2)
beta_par_bef <- list(c(0.6,-.1),c(-0.4,.1))

alpha_par_aft <- c(.01,-.05)
beta_par_aft <- list(c(0.05,.01),c(-0.03,-0.02))

sigma_bef <- 1
sigma_aft <- 1

sigma_bef <- .1
sigma_aft <- .1

sigma_bef <- .01
sigma_aft <- .01

###### ###### ###### ###### ###### ###### ###### ###### 
# Scenario with change only in residual variance:
###### ###### ###### ###### ###### ###### ###### ###### 

# set parameters
alpha_par_bef <- c(-.1,.2)
beta_par_bef <- list(c(-0.3,-.2),c(0.1,.1))

alpha_par_aft <- c(-.1,.2)
beta_par_aft <- list(c(-0.3,-.2),c(0.1,.1))

sigma_bef <- 1
sigma_aft <- .5

# simulate network
set.seed(seeds[1])
net_er <- erdos.renyi.game(20,p.or.m = 76,type = "gnm",directed = TRUE) 
V(net_er)$name <- as.character(seq(1,20,1))# name nodes (characters)
edgelist_er <- get.edgelist(net_er)
nedges_er <- ecount(net_er)

# Simulate time series on the edges for 1st interval
simdata_bef <- gnar_edge_sim(n=100,net=net_er,alphaParams=alpha_par_bef,
                         betaParams=beta_par_bef,
                         sigma = sigma_bef, meann=0,nedges=nedges_er,
                         data_edges = edgelist_er)
# Simulate time series on the edges for 1st interval
simdata_aft <- gnar_edge_sim(n=100,net=net_er,alphaParams=alpha_par_aft,
                             betaParams=beta_par_aft,
                             sigma = sigma_aft, meann=0,nedges=nedges_er,
                             data_edges = edgelist_er)

# run GNAR-edge
globalalpha <- TRUE
alphaOrder <- 2
betaOrder <- rep(2,2)

# under H1 alternative, change
gnar_bef <- gnar_edge_fit(simdata_bef,edgelist_er,alphaOrder = alphaOrder, betaOrder = betaOrder,net = net_er,lead_lag_mat = lead_lag_2sic,
                                        globalalpha = TRUE,lead_lag_weights = FALSE)
est_var_res_bef <- (summary(gnar_bef[[1]])$sigma)^2

gnar_aft <- gnar_edge_fit(simdata_aft,edgelist_er,alphaOrder = alphaOrder, betaOrder = betaOrder,net = net_er,lead_lag_mat = lead_lag_2sic,
                                        globalalpha = TRUE,lead_lag_weights = FALSE)

est_var_res_aft <- (summary(gnar_aft[[1]])$sigma)^2

# Under H0 , no change
gnar_whole <- gnar_edge_fit(cbind(simdata_bef,simdata_aft),edgelist_er,alphaOrder = alphaOrder, betaOrder = betaOrder,net = net_er,lead_lag_mat = lead_lag_2sic,
                                    globalalpha = TRUE,lead_lag_weights = FALSE)
est_sigma_Ho <- (summary(gnar_whole[[1]])$sigma)^2

# statistic
lambda <- 100*log(est_sigma_Ho/est_var_res_bef)+100*log(est_sigma_Ho/est_var_res_aft)
pchisq(lambda,7,lower.tail = FALSE) 

#### METHOD USING R PACKAGE VARDetect
datawhole <- t(cbind(simdata_bef,simdata_aft))
tbss(datawhole,method="group sparse",q=2)

#### METHOD USING MST 
# create list of adjacency matrices (no zeros so ultimately fixed network structure and changing weights only)
adj_list <- list()
for(i in 1:200){
  E(net_er)$weight <- datawhole[i,]
  #E(net_er)$weight <- simdata[,i]
  adj_list[[i]] <- get.adjacency(net_er,attr = 'weight')
}
# Obtain distance matrices for weighted graphs and calculate statistics
dist_mat_frob <- nd.edd(adj_list,out.dist = TRUE)
#dist_mat_spec <- nd.wsd(adj_list,out.dist = TRUE) # Normalized Laplacian matrix contains topological information of a corresponding network via its spectrum. nd.wsd adopts weighted spectral distribution of eigenvalues and brings about a metric via binning strategy.
#dist_mat_diff <- nd.gdd(adj_list,out.dist = TRUE) # Graph Diffusion Distance (nd.gdd) quantifies the difference between two weighted graphs of same size. It takes an idea from heat diffusion process on graphs via graph Laplacian exponential kernel matrices.
mst_edeglist <- mstree(dist_mat_frob$D)
gseg1(200,mst_edeglist,statistics = "all")

gseg2(200,mst_edeglist,statistics = "all")


###### ###### ###### ######
# Scenario with no change
###### ###### ###### ######

# set parameters
alpha_par <- c(-.1,.2)
beta_par <- list(c(-0.3,-.2),c(0.1,.1))

sigma <- 1

# simulate network
set.seed(seeds[1])
net_er <- erdos.renyi.game(20,p.or.m = 76,type = "gnm",directed = TRUE) 
V(net_er)$name <- as.character(seq(1,20,1))# name nodes (characters)
edgelist_er <- get.edgelist(net_er)
nedges_er <- ecount(net_er)

# Simulate time series on the edges 
simdata<- gnar_edge_sim(n=200,net=net_er,alphaParams=alpha_par,
                             betaParams=beta_par,
                             sigma = sigma, meann=0,nedges=nedges_er,
                             data_edges = edgelist_er)
simdata_bef <- simdata[,1:100]
simdata_aft <- simdata[,101:200]

# run GNAR-edge
globalalpha <- TRUE
alphaOrder <- 2
betaOrder <- rep(2,2)

# under H1 alternative, change
gnar_bef <- gnar_edge_fit(simdata_bef,edgelist_er,alphaOrder = alphaOrder, betaOrder = betaOrder,net = net_er,lead_lag_mat = lead_lag_2sic,
                          globalalpha = TRUE,lead_lag_weights = FALSE)
est_var_res_bef <- (summary(gnar_bef[[1]])$sigma)^2

gnar_aft <- gnar_edge_fit(simdata_aft,edgelist_er,alphaOrder = alphaOrder, betaOrder = betaOrder,net = net_er,lead_lag_mat = lead_lag_2sic,
                          globalalpha = TRUE,lead_lag_weights = FALSE)

est_var_res_aft <- (summary(gnar_aft[[1]])$sigma)^2

# Under H0 , no change
gnar_whole <- gnar_edge_fit(simdata,edgelist_er,alphaOrder = alphaOrder, betaOrder = betaOrder,net = net_er,lead_lag_mat = lead_lag_2sic,
                            globalalpha = TRUE,lead_lag_weights = FALSE)
est_sigma_Ho <- (summary(gnar_whole[[1]])$sigma)^2

# statistic
lambda <- 100*log(est_sigma_Ho/est_var_res_bef)+100*log(est_sigma_Ho/est_var_res_aft)
pchisq(lambda,7,lower.tail = FALSE) # p-val>.05 strong no reject H0


###### SUMMARISE RESULTS:
# Scenario with change both on autoregressive param and sigma, p-value=2.361566e-08, reject Ho, presence of change (lambda=48.91131)
# Scenario with change only on autoregressive param, p-value=0.450982, change not detected (lambda=6.790801)
# Scenario with change only on residual variance, p-value=1.143959e-07, reject Ho, presence of change (lambda=45.39931)
# Scenario with no change from generating time series from the same model param, p-value=0.9999975 (lambda=0.102028)

##### IMPORTANT: FOR THE NETWORK TIME SERIES WITH CHANGE POINT, CHECK WHETHER TEST IDENTIFIES ALSO OTHER STATISTICALLY SIGNIFICANT CHANGE POINTS ##### 

simwhole <- cbind(simdata_bef,simdata_aft)
# explore changes for whole time series τ=50,...,150


change_detect_gnaredge <- function(data,time_vec,data_edges,alphaOrder,betaOrder,graph_2dig,whole=TRUE,concat=FALSE){
  nomin_ll_mle <- nomin_ll_c <- c()
  est_var_res_bef <- est_var_res_aft <- c()
  est_sigma_Ho <- c()
  counter <- 1
  for (c in time_vec){# seq(27-2,81-2,6), seq(39-2,72-2,3)
    print(c("iteration ",c))
    if (whole==FALSE){
      data_bef <- data[,1:c]
      data_aft <- data[,(c+1):ncol(data)]
      # difference and normalise
      dtn_bef <- stationary_ts(data_bef)
      data_norm_bef <- t(dtn_bef[[1]])
      data_sd_bef <- dtn_bef[[2]] 
      
      dtn_aft <- stationary_ts(data_aft)
      data_norm_aft <- t(dtn_aft[[1]])
      data_sd_aft <- dtn_aft[[2]] 
      if(concat==TRUE){
        gnar_alpha8_stage1 <- gnar_edge_fit(cbind(data_norm_bef,data_norm_aft),data_edges,alphaOrder = alphaOrder, betaOrder = betaOrder,net = graph_2dig,lead_lag_mat = "none",
                                            globalalpha = TRUE,lead_lag_weights = FALSE)
        #ll_all <- logLik(gnar_alpha8_stage1$mod)
        #ll_all <- sum(dnorm(x = gnar_alpha8_stage1$y, mean = predict(gnar_alpha8_stage1$mod), sd = summary(gnar_alpha8_stage1$mod)$sigma,log = TRUE))
        est_sigma_Ho <- c(est_sigma_Ho,(summary(gnar_alpha8_stage1[[1]])$sigma)^2)
      }
    }else{
      if (counter==1){ 
        # difference and normalise
        dtn <- stationary_ts(data)
        data_norm <- t(dtn[[1]])
        data_sd <- dtn[[2]] # keep sd of data_train
        
        gnar_alpha8_stage1 <- gnar_edge_fit(data_norm,data_edges,alphaOrder = alphaOrder, betaOrder = betaOrder,net = graph_2dig,lead_lag_mat = "none",
                                            globalalpha = TRUE,lead_lag_weights = FALSE)
        #ll_all <- logLik(gnar_alpha8_stage1$mod)
        #ll_all <- sum(dnorm(x = gnar_alpha8_stage1$y, mean = predict(gnar_alpha8_stage1$mod), sd = summary(gnar_alpha8_stage1$mod)$sigma,log = TRUE))
        est_sigma_Ho <- (summary(gnar_alpha8_stage1[[1]])$sigma)^2
        
      }
      data_norm_bef <- data_norm[,1:(c-1)] # indices for time stamps change due to differencing (lost one value)
      data_norm_aft <- data_norm[,c:length(data_norm)]
      counter <- counter+1
    }
    # run GNAR-edge
    gnar_bef_alpha8_stage1 <- gnar_edge_fit(data_norm_bef,data_edges,alphaOrder = alphaOrder, betaOrder = betaOrder,net = graph_2dig,lead_lag_mat = lead_lag_2sic,
                                            globalalpha = TRUE,lead_lag_weights = FALSE)
    
    gnar_aft_alpha8_stage1 <- gnar_edge_fit(data_norm_aft,data_edges,alphaOrder = alphaOrder, betaOrder = betaOrder,net = graph_2dig,lead_lag_mat = lead_lag_2sic,
                                            globalalpha = TRUE,lead_lag_weights = FALSE)
    
    # log-likel for change point with logLik() function (for a model fitted by maximum likelihood)
    ll_bef <- logLik(gnar_bef_alpha8_stage1$mod)
    ll_aft <- logLik(gnar_aft_alpha8_stage1$mod)
    ll_mle <- ll_bef+ll_aft
    
    # log-likel for change point
    ll_bef_c <- sum(dnorm(x = gnar_bef_alpha8_stage1$y, mean = predict(gnar_bef_alpha8_stage1$mod), sd = summary(gnar_bef_alpha8_stage1$mod)$sigma,log = TRUE))
    ll_aft_c <- sum(dnorm(x = gnar_aft_alpha8_stage1$y, mean = predict(gnar_aft_alpha8_stage1$mod), sd = summary(gnar_aft_alpha8_stage1$mod)$sigma,log = TRUE))
    ll_c <- ll_bef_c+ll_aft_c
    
    nomin_ll_mle <- c(nomin_ll_mle,ll_mle)
    nomin_ll_c <- c(nomin_ll_c,ll_c)
    
    # estimated variance of residuals
    est_var_res_bef <- c(est_var_res_bef,(summary(gnar_bef_alpha8_stage1[[1]])$sigma)^2)
    est_var_res_aft <- c(est_var_res_aft,(summary(gnar_aft_alpha8_stage1[[1]])$sigma)^2)
  }
  if(whole==FALSE & concat==TRUE){
    return(list(est_var_res_bef,est_var_res_aft,est_sigma_Ho,nomin_ll_mle,nomin_ll_c))
  }else if(whole==TRUE){
    return(list(est_var_res_bef,est_var_res_aft,nomin_ll_mle,nomin_ll_c))
  }
}



##### IMPORTANT: TRY FOR BIGGER LAG SIZES AND STAGE NEIGHBOURS TO SEE IF SIMILAR SIZES OF LAMBDAS ARE OBTAINED 