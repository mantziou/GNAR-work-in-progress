########################################
######### MC sampling of edges ######### 
########################################

# LOAD DATA
load("/YOUR_DIR/data_allver_alignednet.RData")

nodes_only <- FALSE
gr <- TRUE
mc_nodes <- list()
data_vec_list <- list()
gdp_sum <- c()
pred_ver <- lapply(1:length(1:50), function(x) vector("list",9))
gdp_sum_mat <- lapply(1:length(1:50),function(x) matrix(nrow = 9,ncol = 4))
gdp_sum_rel_err <- lapply(1:length(1:50),function(x) matrix(nrow = 9,ncol = 4))
set.seed(100)
for (i in 1:50){
  ### RANDOM SAMPLE TO REMOVE
  removepaym_pos <- sample(seq(1:833),125)
  
  
  #### REMOVAL STEP
  data_edges3 <-  get(paste("prep_data_rev_5uni",sep = ""))$data_edges[-removepaym_pos,]
  # not in common nodes (after removing edges, removing also some nodes)
  removegdp_pos <- as.numeric(setdiff(get(paste("prep_data_rev_5uni",sep = ""))$data_nodes,unique(c(unique(data_edges3[,2]),unique(data_edges3[,1])))))
  
  if (length(removegdp_pos)!=0){
    CPA_node3 <- get(paste("prep_data_rev_5uni",sep = ""))$CPA_node[-removegdp_pos,]
    data_nodes3 <- get(paste("prep_data_rev_5uni",sep = ""))$data_nodes[-removegdp_pos]
  }else{
    CPA_node3 <- get(paste("prep_data_rev_5uni",sep = ""))$CPA_node
    data_nodes3 <- get(paste("prep_data_rev_5uni",sep = ""))$data_nodes
  }
  graph3 <-  graph_from_data_frame(as.data.frame(data_edges3),directed = TRUE,vertices = data_nodes3)
  nedges3 <-  ecount(graph3)
  nnodes3 <- vcount(graph3)
  n_edges_nodes3 <-  ecount(graph3)+vcount(graph3)
  mc_nodes[[i]] <- data_nodes3
  
  for (rev in 5:6){
    df_pay3 <- get(paste("prep_data_rev_",rev,"uni",sep = ""))$df_pay[-removepaym_pos,]
    if (length(removegdp_pos)==0){
      df_gva3 <- get(paste("prep_data_rev_",rev,"uni",sep = ""))$df_gva
      ts_pay_gdp3 <-  get(paste("prep_data_rev_",rev,"uni",sep = ""))$ts_pay_gdp[-removepaym_pos,]
    }else{
      df_gva3 <- get(paste("prep_data_rev_",rev,"uni",sep = ""))$df_gva[-removegdp_pos,]  
      ts_pay_gdp3 <-  get(paste("prep_data_rev_",rev,"uni",sep = ""))$ts_pay_gdp[-c(removepaym_pos,removegdp_pos+833),]
    }
    
    assign(paste("mc_",rev,sep = ""),list(df_pay=df_pay3,df_gva=df_gva3,ts_pay_gdp=ts_pay_gdp3,CPA_node=CPA_node3,graph=graph3,nnodes=nnodes3,nedges=nedges3,
                                                          n_edges_nodes=n_edges_nodes3,data_edges=data_edges3,data_nodes=data_nodes3))
    
  }
  
  ##### GET GROWTH RATES AND STL DECOMPOSITION
  data <- mc_5$ts_pay_gdp
  data <- growthrates(data)

  # stl
  res_stl <- prepr_stl(data,str_sub(colnames(data)[1],6,7),str_sub(tail(colnames(data),n=1),6,7),str_sub(colnames(data)[1],1,4),str_sub(tail(colnames(data),n=1),1,4),nts=1)
  data_train_stl <- res_stl$data_train_stl
  
  ###### TRUE VAL OF SUBSEQUENT RELEASE TO EVALUATE ERROR
  data2 <- mc_6$ts_pay_gdp
  # check that two consecutive data revisions, same graph
  ifelse(all(rownames(data2)==rownames(data)),print("same graph"),warning("graph not the same between 5 and rev 6"))
  data_vec <- data2[,ncol(data)+1]
  if(nodes_only==TRUE){
    data_vec <- data_vec[(mc_5$nedges+1):mc_5$n_edges_nodes]
  }
  data_vec_list[[i]] <- data_vec
  gdp_sum[i] <- sum(tail(data_vec,mc_6$nnodes))

  ##### RUN GNAR MODEL FOR VARIOUS LAG AND STAGE
  for (lagi in 1:9){
    ind_stag <- 1
    for (stag in 0:3){
      print(c("data iter",i,"Lag iter",lagi," stag iter ",stag))
      alphaOrder <- lagi
      betaOrder <- rep(stag,lagi)
      gammaOrder <- lagi
      deltaOrder <- rep(stag,lagi)
      
      # fit GNAR-x
      fit_train <- gnar_x_fit(data_train_stl,mc_5$nnodes,mc_5$nedges,mc_5$n_edges_nodes,mc_5$data_edges,mc_5$data_nodes,
                              alphaOrder,betaOrder,gammaOrder,deltaOrder,mc_5$graph,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now = FALSE,lag_0_sep = FALSE,
                              nodes_only = nodes_only)
      #print(summary(fit_train$mod)$coefficients[,4] < 0.05)
      # Predict
      pred <- gnar_x_predict(fit_train,data_train_stl,alphaOrder,betaOrder,
                             gammaOrder,deltaOrder,mc_5$n_edges_nodes,mc_5$nnodes,mc_5$nedges,1,fit_train$wei_mat,fit_train$data_loc_mat,
                             pay_now=FALSE,pay_data_now=NULL)
      
      pred_adj <- pred[,lagi+1]+res_stl$pred_trend+res_stl$seas_comp_nts
      
      if (gr==TRUE){
        pred_adj <- mc_5$ts_pay_gdp[,ncol(mc_5$ts_pay_gdp)]*(1+pred_adj)
      }
      
      if (nodes_only==TRUE){
        pred_adj <- pred_adj[(mc_5$nedges+1):mc_5$n_edges_nodes]
      }
      
      pred_ver[[i]][[lagi]][[ind_stag]] <- pred_adj
      gdp_sum_mat[[i]][lagi,ind_stag] <- sum(tail(as.vector(pred_adj),n=mc_5$nnodes)) 
      
      gdp_sum_rel_err[[i]][lagi,ind_stag] <- abs((gdp_sum_mat[[i]][lagi,ind_stag]-gdp_sum[i])/gdp_sum[i])
      
      ind_stag <- ind_stag+1
    }
  }
}

bestmod2 <- t(sapply(gdp_sum_rel_err2, function(x) which(x==min(x),arr.ind = TRUE)))
minbestmod2 <- sapply(gdp_sum_rel_err2, function(x) min(x))
properror <- data.frame(x=0.00001/minbestmod2)
ggplot(properror,aes(x=x))+geom_histogram(color="blue",fill="white")+xlab("error with Pearson's cutoff/error with random cutoff")
