################################################################################
######## CUTOFF CHOICE USING MODEL AVERAGING OF RELATIVE ERROR ##########
################################################################################

ind_stag <- 1
nodes_only <- FALSE
gr <- TRUE

pred_ver <- lapply(1:length(1:9), function(x) vector("list",9))
pred_trend_list <- lapply(1:9, function(x) list())
seas_comp_list <- lapply(1:9, function(x) list())
data_vec_list <- lapply(1:9, function(x) list())
ma_rel_err <- lapply(1:9,function(x)list())
ma_rel_err_all <- list()
gdp_sum_rel_err_min <- lapply(1:9,function(x)list())
gdp_sum_rel_err_min_all <- list()

thres_p_vec <- c(.3,.5,.6,.7)
gdp_sum_mat <- lapply(1:length(1:9),function(x) matrix(nrow = 9,ncol = 4))
gdp_sum_rel_err <- lapply(1:length(1:9),function(x) matrix(nrow = 9,ncol = 4))
for(th in 1:4){
  # load data
  load(paste("/YOUR_DIR/sparse_data_allver_nodedgerem_rev_pears",thres_p_vec[th],".RData"))
  for (rev in 1:9){
    prep_data <- get(paste("sparse_",rev,"_nodedgerem_thres",sep = ""))
    data <- prep_data$ts_pay_gdp_2

    # growth rates
    data <- growthrates(data)
    
    # stl
    res_stl <- prepr_stl(data,str_sub(colnames(data)[1],6,7),str_sub(tail(colnames(data),n=1),6,7),str_sub(colnames(data)[1],1,4),str_sub(tail(colnames(data),n=1),1,4),nts=1)
    data_train_stl <- res_stl$data_train_stl
    pred_trend_list[[rev]] <- res_stl$pred_trend
    seas_comp_list[[rev]] <- res_stl$seas_comp_nts
    simtrainvar <- t(data_train_stl)
    
    
    for (lagi in 1:9){
      ind_stag <- 1
      for (stag in 0:3){
        print(c("data rev",rev,"Lag iter",lagi," stag iter ",stag))
        alphaOrder <- lagi
        betaOrder <- rep(stag,lagi)
        gammaOrder <- lagi
        deltaOrder <- rep(stag,lagi)
        
        # fit GNAR-x
        fit_train <- gnar_x_fit(data_train_stl,prep_data$nnodes_2,prep_data$nedges_2,prep_data$n_edges_nodes_2,prep_data$data_edges_2,prep_data$data_nodes_2,
                                alphaOrder,betaOrder,gammaOrder,deltaOrder,prep_data$graph_2,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now = FALSE,lag_0_sep = FALSE
        )
        # Predict
        pred <- gnar_x_predict(fit_train,data_train_stl,alphaOrder,betaOrder,
                               gammaOrder,deltaOrder,prep_data$n_edges_nodes_2,prep_data$nnodes_2,prep_data$nedges_2,1,fit_train$wei_mat,fit_train$data_loc_mat,
                               pay_now=FALSE,pay_data_now=NULL)
        
        pred_adj <- pred[,lagi+1]+res_stl$pred_trend+res_stl$seas_comp_nts
        
        if (gr==TRUE){
          pred_adj <- prep_data$ts_pay_gdp_2[,ncol(prep_data$ts_pay_gdp_2)]*(1+pred_adj)
        }
        
        if (nodes_only==TRUE){
          pred_adj <- pred_adj[(prep_data$nedges_2+1):prep_data$n_edges_nodes_2]
        }
        
        pred_ver[[rev]][[lagi]][[ind_stag]] <- pred_adj
        data2 <- get(paste("sparse_",rev+1,"_nodedgerem_thres",sep = ""))$ts_pay_gdp_2
        # check that two consecutive data revisions, same graph
        ifelse(all(rownames(data2)==rownames(data)),print("same graph"),warning(paste("graph not the same between rev ",rev," and rev ",rev+1,sep="")))
        
        if (rev!=4){
          data_vec <- data2[,ncol(data)+1]
        }else{
          data2 <- get(paste("sparse_",rev+2,"_nodedgerem_thres",sep = ""))$ts_pay_gdp
          data_vec <- data2[,ncol(data)+1]
        }
        
        if(nodes_only==TRUE){
          data_vec <- data_vec[(prep_data$nedges_2+1):prep_data$n_edges_nodes_2]
        }
        data_vec_list[[rev]] <- data_vec
        
        gdp_sum <- sum(tail(data_vec_list[[rev]],n=prep_data$nnodes_2))
        gdp_sum_mat[[rev]][lagi,ind_stag] <- sum(tail(as.vector(pred_ver[[rev]][[lagi]][[ind_stag]]),n=prep_data$nnodes_2)) 
        gdp_sum_rel_err[[rev]][lagi,ind_stag] <- abs((gdp_sum_mat[[rev]][lagi,ind_stag]-gdp_sum)/gdp_sum)
        
        ind_stag <- ind_stag+1
      } # end stag
    } # end lagi
    
    ma <- mean(gdp_sum_mat[[rev]])
    ma_rel_err[[rev]] <- abs((ma-gdp_sum)/gdp_sum)
    gdp_sum_rel_err_min[[rev]] <- min(gdp_sum_rel_err[[rev]])
  } # end rev
  ma_rel_err_all[[th]] <- unlist(ma_rel_err)
  gdp_sum_rel_err_min_all[[th]] <- unlist(gdp_sum_rel_err_min)
}

gdp_sum_rel_err_min_all2 <- gdp_sum_rel_err_min_all

gdp_sum_rel_err_min_mat <- sapply(gdp_sum_rel_err_min_all2,rbind)
colnames(gdp_sum_rel_err_min_mat) <- thres_p_vec
rownames(gdp_sum_rel_err_min_mat) <- paste("rev",c(1:9))

ma_rel_err_mat <- sapply(ma_rel_err_all,rbind)
colnames(ma_rel_err_mat) <- thres_p_vec
rownames(ma_rel_err_mat) <- paste("rev",c(1:9))

plot(gdp_sum_rel_err_min_mat[i,],type = "l")
plot(ma_rel_err_mat[i,],type = "l")

## ADD THE RESULTS FOR CUTOFF 0.4
ma_rel_err_0.4 <- gdp_sum_rel_err_min_0.4 <- c()
for (rev in 1:9){
  gdp_sum <- sum(tail(data_vec_list[[rev]],32))
  ma <- mean(gdp_sum_mat[[rev]])
  ma_rel_err_0.4[rev] <- abs((ma-gdp_sum)/gdp_sum)
  gdp_sum_rel_err_min_0.4[rev] <- min(gdp_sum_rel_err[[rev]])
}

ma_rel_err_all_new <- append(ma_rel_err_all,list(ma_rel_err_0.4),after = 1)
gdp_sum_rel_err_min_all_new <- append(gdp_sum_rel_err_min_all,list(gdp_sum_rel_err_min_0.4),after = 1)

gdp_sum_rel_err_min_mat <- sapply(gdp_sum_rel_err_min_all_new,rbind)
colnames(gdp_sum_rel_err_min_mat) <- c(0.3,.4,.5,.6,.7)
rownames(gdp_sum_rel_err_min_mat) <- paste("rev",c(1:9))

ma_rel_err_mat <- sapply(ma_rel_err_all_new,rbind)
colnames(ma_rel_err_mat) <- c(0.3,.4,.5,.6,.7)
rownames(ma_rel_err_mat) <- paste("rev",c(1:9))

for (rev in 1:9){
  rmse_df <- data.frame(gdp_sum_rel_err_min_mat[rev,])
  colnames(rmse_df) <- c("error")

  assign(paste("p_",rev,sep = ""),ggplot(rmse_df, aes(x = c(0.3,.4,.5,.6,.7), y = error)) +  geom_line() + geom_point(size=1)+
           theme(axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),
                 axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.text = element_text(size=12),legend.title = element_text(size=12),
                 plot.title = element_text(hjust = 0.5))+
           scale_x_continuous("cutoff")+ylab("Min rel error")+ggtitle(paste("Release ",substr(colnames(dt)[9:18],start = 8,stop = 14)[rev],sep = "")))
}
library(patchwork)
p_1+p_2+p_3
p_4+p_5+p_6
p_7+p_8+p_9
library(gridExtra)
grid.arrange(p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,nrow=3)

