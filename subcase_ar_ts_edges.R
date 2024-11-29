library(igraph)
library(SparseM)
library(Matrix)
library(faux)
library(lqmm)
library(matrixcalc)
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######  ######  ######  ###### 
###### SUBCASE of Extended of GNAR model (ONLY NODAL TIME SERIES REGRESSION) ###### 
###### ###### ###### ###### ###### ###### ###### ###### ######  ###### ###### ######  ######  ###### 

response_vec<-function(ts_data,alphaOrder){
  # ts_data: dataframe with rows edges and columns time
  # alphaOrder: max lag
  # returns: vectorised response (concatenate ts for multiple time series)
  predt <- ncol(ts_data)-alphaOrder
  yvec <- NULL
  for(en in 1:nrow(ts_data)){
    yvec <- c(yvec, ts_data[en,((alphaOrder+1):(predt+alphaOrder))])
  }
  return(unlist(yvec))
}

design_mat<-function(ts_data,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,alphaOrder,betaOrder,gammaOrder,deltaOrder,
                     net,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now=TRUE){
  
  # ts_data: matrix with rows edges and columns time
  # nnodes: number of nodes
  # nedges: number of edges
  # n_edges_nodes: total number of rows of ts_data==nnodes+nedges
  # data_edges: dataframe with edge list, column one node from which edge starts (character type) column two node edge ends (character type)
  # data_nodes: vector of nodes
  # alphaOrder: max lag for autoregressive
  # betaOrder: vector of length as number of lags alphaOrder considered. Each entry give max stage neighbour to be considered for each lag
  # gammaOrder: max lag for cross-autoregressive
  # deltaOrder: vector of length as number of lags gammaOrder considered. Each entry give max stage neighbour to be considered for each cross-lag
  # net: network resulting from edge list data_edges
  # lead_lag_mat: matrix with rows and col corresponding to edge time series, entries showing leadingness of each time series over the other
  # globalalpha: if True, non edge-specific alpha parameters, but only lag specific alpha param
  # lead_lag_weights: whether lead-lag matrix used as weights or not
  
  if (pay_now==TRUE){
    ilagstart <- 0
  }else{
    ilagstart <- 1
  }
  
  n_edges_nodes <- nrow(ts_data)
  
  predt <- ncol(ts_data)-alphaOrder # lenght of time series observations considered when removing max lag
  # CREATE NAMES FOR COLS OF DESIGN MAT ACCORDING TO model PARAMETERS
  parNames <- NULL
  for (ilag in ilagstart:alphaOrder) {
    if (ilag==0){
      parNames <- c(parNames, paste("gamma", ilag, sep = ""))
      if (deltaOrder[ilag+1] > 0){ # ilag = 0 , but indexing in R starts at 1
        for (stag in 1:deltaOrder[ilag+1]){
          parNames <- c(parNames, paste("delta", ilag, ".stage",
                                        stag, sep = ""))
        }
      }
    }else{
      if(globalalpha){
        parNames <- c(parNames, paste("alpha", ilag, sep = ""))
      }else{
        for(edg in 1:nedges){
          parNames <- c(parNames, paste("alpha", ilag, "edge", edg, sep=""))
        }
      }
      if (betaOrder[ilag] > 0) {
        for (stag in 1:betaOrder[ilag]) {
          parNames <- c(parNames, paste("beta", ilag, ".stage",
                                        stag, sep = ""))
        }
      }
      parNames <- c(parNames, paste("gamma", ilag, sep = ""))
      if (pay_now==TRUE){
        ilagstag <- ilag+1
      }else{
        ilagstag <- ilag
      }
      if (deltaOrder[ilagstag] > 0){
        for (stag in 1:deltaOrder[ilagstag]){
          parNames <- c(parNames, paste("delta", ilag, ".stage",
                                        stag, sep = ""))
        }
      }
    }
  }
  
  # CREATE EMPTY DESIGN MATRIX WITH 0s
  dmat <- matrix(0, nrow = predt * n_edges_nodes, ncol = length(parNames), dimnames = list(NULL, parNames))
  
  # FILL DESIGN MATRIX COLUMNS CORRESPONDING TO ALPHA PARAMETERS WITH LAGGED TIME SERIES FOR EACH EDGE AND EACH NODE
  for (r in 1:n_edges_nodes) {
    for (ilag in 1:alphaOrder) {
      if(globalalpha){
        alphaLoc <- paste("alpha", ilag, sep = "")
        #}else{
        #  alphaLoc <- paste("alpha", ilag, "edge", r, sep = "")
      }
      dmat[((predt * (r - 1) + 1):(predt * r)), alphaLoc] <-  unlist(ts_data[r,((alphaOrder +1 - ilag):(predt + (alphaOrder - ilag)))])
    }
  }
  
  
  # FILL DESIGN MATRIX COLUMNS FOR BETA PARAMETERS
  if (sum(betaOrder) > 0) {
    wei_mat <- lapply(1:max(betaOrder),function(x)Matrix(0,n_edges_nodes,n_edges_nodes,doDiag = FALSE))
    
    
    for (al in 1:n_edges_nodes){
      #print(c("edge",al))
      
      # GET NEIGHBOUR EDGES and NODES
      if (al <= nedges){
        nodess <- c(data_edges[al,1],data_edges[al,2]) # get nodes involved in edg
        nei <- edge_neighbors(net,max(betaOrder),nodess) # for these nodess get the r-stage neighbours
      }else{
        nei <- node_neighbors(net,rownames(ts_data)[al],max(betaOrder))
      }
      
      
      
      if (!lead_lag_weights){# if not using lead-lag matrix, equally weight neighbour edges (mean)
        wei <- sapply(nei, length)
        wei <- lapply(wei, function(x) rep(1/x,x))
      }else{ # else use lead-lag matrix
        if (al <= nedges){ # lead lag option for edges only
          wei <- list()
          if ((!is.null(nei)) & (length(nei) > 0)) {
            for (ilagg in 1:alphaOrder){
              if (betaOrder[ilagg]>=1){
                for (stagg in 1:betaOrder[ilagg]){
                  edg_loc <- edge_matching(data_edges, nei[[stagg]])
                  wei[[stagg]] <- lead_lag_mat[edg_loc,al]
                }
              }
            }
          }
        }
      }
      
      
      
      # IF THERE ARE NEIGHBOURS GET IN IF
      if ((!is.null(nei)) & (length(nei) > 0)) { # check that are non empty lists
        
        
        
        
        for (ilag in 1:alphaOrder){
          if (betaOrder[ilag]>=1){
            
            for (stag in 1:betaOrder[ilag]){
              #print(c("ilag ",ilag, " stage ",stag) )
              betaLoc <- paste("beta", ilag, ".stage", stag, sep = "")
              
              if ((!is.null(nei[[stag]]))){
                
                # CASE WHERE stag-STAGE NEIGHBOURS MORE THAN ONE (TO NORMALISE RESPECTIVELY)
                if (length(nei[[stag]]) > 1) {
                  
                  # get indexing/locations of neighbours of current "al" to cut df accordingly
                  if (al <= nedges){
                    data_loc <- edge_matching(data_edges, nei[[stag]])  
                  }else{
                    #data_loc <- which(data_nodes %in% nei[[stag]])+nedges
                    data_loc <- which(V(net) %in% nei[[stag]])+nedges
                  }
                  
                  
                  # CUT TIME SERIES MAT FOR COLUMNS CORRESPONDING TO NEIGHBOUR NODES AND ROWS CORRESP TO LAGGED TIME SERIES
                  vts.cut <- ts_data[data_loc,((alphaOrder + 1 - ilag):(predt + (alphaOrder - ilag)))]
                  
                  # FOR CASE WHERE VTS.CUT CONTAINS NA VALUES -> IMPUTATION
                  # for each col in cut of time series matrix transform the multiple ts to single ts by getting their weighted sum
                  for (t in 1:ncol(vts.cut)) {
                    if (any(is.na(vts.cut[ ,t]))) {
                      if (all(is.na(vts.cut[ ,t]))) {
                        #if there are no neighbours left at any time point, set to zero
                        vts.cut[,t] <- 0
                      }
                      else {
                        new.wei <- wei[[stag]][which(!is.na(vts.cut[,t]))]
                        if (!lead_lag_weights){ # do the re-normalisation only if not lead-lag weights
                          new.wei <- new.wei/sum(new.wei) # re-normalise after removing weights for NA
                        }
                        sub.val <- vts.cut[which(!is.na(vts.cut[,t])),t] %*% new.wei # create from values of multiple time series without NAs at time t, a single value at time t by the weighted sum
                        vts.cut[which(is.na(vts.cut[,t])),t] <- sub.val # fill in NA with weighted sum of time series at time t without NAs
                      }
                    }
                  }
                  
                  # ASSIGN TO DESIGN MAT THE SUM OF WEIGHTED TIME SERIES FOR LAG
                  dmat[((predt * (al - 1) + 1):(predt * al)), betaLoc] <- as.vector(t(vts.cut) %*% wei[[stag]])
                  wei_mat[[stag]][data_loc,al] <- wei[[stag]]
                }
                
                # CASE WHERE ONLY ONE stag-STAGE NEIGHBOUR EDGE  (TO NORMALISE RESPECTIVELY)
                else{
                  
                  # get indexing/locations of neighbours of current "al" to cut df accordingly
                  if (al <= nedges){
                    data_loc <- edge_matching(data_edges, nei[[stag]])  
                  }else{
                    #data_loc <- which(data_nodes %in% nei[[stag]])+nedges
                    data_loc <- which(V(net) %in% nei[[stag]])+nedges
                  }
                  
                  # CASE WHERE stag-STAGE NEIGHBOUR EDGES IS EXACTLY ONE EDGE (NO NEED TO NORMALISE) AND NON NA VALUES
                  if ((length(nei[[stag]]) == 1) & (!is.na(nei[[stag]]))) {
                    vts.cut <- unlist(ts_data[data_loc,((alphaOrder + 1 - ilag):(predt + (alphaOrder - ilag)))])
                    #and if this is missing at any time point, set to zero (cannot fill in with mean value of neighbours as previously as no other neighbours)
                    vts.cut[is.na(vts.cut)] <- 0
                    dmat[((predt * (al - 1) + 1):(predt * al)), betaLoc] <- as.vector(t(vts.cut) * wei[[stag]])
                  }
                  
                  # CASE WHERE stag-STAGE NEIGHBOUR NAN VALUES
                  else {
                    # set all to 0
                    dmat[((predt * (al - 1) + 1):(predt * al)), betaLoc] <- 0
                  }
                  wei_mat[[stag]][data_loc,al] <- wei[[stag]]
                }
              } # end if stag specific neighbours existing
              
              # case where no neighbours of stag size
              else{
                dmat[((predt * (al - 1) + 1):(predt * al)), betaLoc] <- 0
              }
              
            } # end for stage
          } # end if check at least one neighbour considered
        } # end for ilag
        
        
        
      } # end if neighbours existing
      
      
      else{
        for (ilag in 1:alphaOrder){
          for (stag in 1:betaOrder[ilag]){
            betaLoc <- paste("beta", ilag, ".stage", stag, sep = "")
            dmat[((predt * (al - 1) + 1):(predt * al)),betaLoc] <- 0
          }
        }
      }
      
    } # END FOR LOOP OVER EDGES&NODES
  }else{ # END IF
    wei_mat <- NULL # NOTE : if stages of neighbours for beta are 0 then this null matrix will cause issues with delta parameters construction
  }
  
  
  # FILL DESIGN MATRIX COLUMNS CORRESPONDING TO GAMMA PARAMETERS WITH LAGGED TIME SERIES FOR EACH EDGE AND EACH NODE
  data_loc_mat <- Matrix(0,n_edges_nodes,n_edges_nodes,doDiag = FALSE)
  for (r in 1:n_edges_nodes) {
    #print(c("r edge/node ",r))
    for (ilag in ilagstart:gammaOrder) {
      gammaLoc <- paste("gamma", ilag, sep = "")
      if (r <= nedges){
        if (pay_now==TRUE){
          if (ilag==0){
            dmat[((predt * (r - 1) + 1):(predt * r)), gammaLoc] <- 0
          }else{
            nodeloc <- which(data_nodes %in% data_edges[r,])+nedges # location of nodes relevant to edge, alternative:which(rownames(ts_data) %in% data_edges[r,])
            vts.cut.edge <-  ts_data[nodeloc,((gammaOrder +1 - ilag):(predt + (gammaOrder - ilag)))]
            if (length(nodeloc)==1){
              data_loc_mat[nodeloc,r] <- 1
              dmat[((predt * (r - 1) + 1):(predt * r)), gammaLoc] <- as.vector(t(vts.cut.edge)) # c(.5,.5) equally weight the two nodes corresponding to edge r
            }else{
              data_loc_mat[nodeloc,r] <- c(0.5,0.5) 
              dmat[((predt * (r - 1) + 1):(predt * r)), gammaLoc] <- as.vector(t(vts.cut.edge) %*% c(.5,.5)) # c(.5,.5) equally weight the two nodes corresponding to edge r
            }
          }
        }else{
          nodeloc <- which(data_nodes %in% data_edges[r,])+nedges # location of nodes relevant to edge, alternative:which(rownames(ts_data) %in% data_edges[r,])
          vts.cut.edge <-  ts_data[nodeloc,((gammaOrder +1 - ilag):(predt + (gammaOrder - ilag)))]
          if (length(nodeloc)==1){
            data_loc_mat[nodeloc,r] <- 1
            dmat[((predt * (r - 1) + 1):(predt * r)), gammaLoc] <- as.vector(t(vts.cut.edge)) # c(.5,.5) equally weight the two nodes corresponding to edge r
          }else{
            data_loc_mat[nodeloc,r] <- c(0.5,0.5) 
            dmat[((predt * (r - 1) + 1):(predt * r)), gammaLoc] <- as.vector(t(vts.cut.edge) %*% c(.5,.5)) # c(.5,.5) equally weight the two nodes corresponding to edge r
          }
        }
      }else{
        edgeloc <- which(apply(data_edges,1,function(x) rownames(ts_data)[r] %in% x))
        if (length(edgeloc)>1){
          edgewei <- rep(1/length(edgeloc),length(edgeloc)) # equally weight edges related to node r
          data_loc_mat[edgeloc,r] <- edgewei
          vts.cut.node <-  ts_data[edgeloc,((gammaOrder +1 - ilag):(predt + (gammaOrder - ilag)))]
          dmat[((predt * (r - 1) + 1):(predt * r)), gammaLoc] <- as.vector(t(vts.cut.node) %*% edgewei)
        }else if (length(edgeloc)==1){
          data_loc_mat[edgeloc,r] <- 1
          vts.cut.node <-  ts_data[edgeloc,((gammaOrder +1 - ilag):(predt + (gammaOrder - ilag)))]
          dmat[((predt * (r - 1) + 1):(predt * r)), gammaLoc] <- as.vector(t(vts.cut.node))
        }else{
          dmat[((predt * (r - 1) + 1):(predt * r)), gammaLoc] <- 0
        }
      }
    }
  }
  
  
  # FILL DESIGN MATRIX COLUMNS FOR DELTA PARAMETERS
  if (sum(deltaOrder) > 0) {
    
    for (al in 1:n_edges_nodes){
      
      
      # USE wei_mat OBTAINED FOR BETA PARAM THAT CONTAINS INFO ABOUT NEIGHBOURS AND WEIGHTS
      
      for (ilag in ilagstart:alphaOrder){
        if (pay_now==TRUE){
          ilagstag <- ilag+1
        }else{
          ilagstag <- ilag
        }
        if (deltaOrder[ilagstag]>=1){
          
          for (stag in 1:deltaOrder[ilagstag]){
            deltaLoc <- paste("delta", ilag, ".stage", stag, sep = "")
            
            # get indexing/locations of neighbours of current "al" to cut df accordingly
            if (al <= nedges){
              if (pay_now==TRUE){
                if (ilag==0){
                  nei_loc <- NULL # set it to NULL to get in else and assign 0 to dmat
                }else{
                  data_loc <- which(data_nodes %in% data_edges[al,])+nedges # location of nodes relevant to edge, alternative:which(rownames(ts_data) %in% data_edges[r,])
                  if (length(data_loc)>1){ # case where more than one nodes involved for the al edge
                    nei_loc <- which(apply(wei_mat[[stag]][(nedges+1):n_edges_nodes,data_loc],1,sum)!=0) + nedges
                    nei_loc <- setdiff(nei_loc,data_loc) # exclude the nodes that are involved in data_loc from current al
                  }else if (length(data_loc)==1) { # if exactly one node involved
                    nei_loc <- which(wei_mat[[stag]][(nedges+1):n_edges_nodes,data_loc]!=0) + nedges
                  }
                }
              }else{
                data_loc <- which(data_nodes %in% data_edges[al,])+nedges # location of nodes relevant to edge, alternative:which(rownames(ts_data) %in% data_edges[r,])
                if (length(data_loc)>1){ # case where more than one nodes involved for the al node
                  nei_loc <- which(apply(wei_mat[[stag]][(nedges+1):n_edges_nodes,data_loc],1,sum)!=0) + nedges
                  nei_loc <- setdiff(nei_loc,data_loc) # exclude the nodes that are involved in data_loc from current al
                }else if (length(data_loc)==1) { # if exactly one node involved
                  nei_loc <- which(wei_mat[[stag]][(nedges+1):n_edges_nodes,data_loc]!=0) + nedges
                }
              }
            }else{
              data_loc <- which(apply(data_edges,1,function(x) rownames(ts_data)[al] %in% x))
              if (length(data_loc)>1){ # case where more than one edge involved for the al node
                nei_loc <- which(apply(wei_mat[[stag]][1:nedges,data_loc],1,sum)!=0) 
                nei_loc <- setdiff(nei_loc,data_loc) # exclude the edges that are involved in data_loc from current al
              }else if (length(data_loc)==1) { # if exactly one edge involved
                nei_loc <- which(wei_mat[[stag]][1:nedges,data_loc]!=0) 
              }else{ # if node is not involved in any edge, also no neighbouring edges
                nei_loc <- NULL
              }
            }
            
            if (length(nei_loc)>1){ # if there are neighbours (at least one)
              
              wei_nei_loc <- rep(1/length(nei_loc),length(nei_loc))
              wei_mat[[stag]][nei_loc,al] <- wei_nei_loc
              
              # CUT TIME SERIES MAT FOR COLUMNS CORRESPONDING TO NEIGHBOUR NODES AND ROWS CORRESP TO LAGGED TIME SERIES
              vts.cut <- ts_data[nei_loc,((alphaOrder + 1 - ilag):(predt + (alphaOrder - ilag)))]
              
              
              # FOR CASE WHERE VTS.CUT CONTAINS NA VALUES -> IMPUTATION
              # for each col in cut of time series matrix transform the multiple ts to single ts by getting their weighted sum
              for (t in 1:ncol(vts.cut)) {
                if (any(is.na(vts.cut[ ,t]))) {
                  if (all(is.na(vts.cut[ ,t]))) {
                    #if there are no neighbours left at any time point, set to zero
                    vts.cut[,t] <- 0
                  }
                  else {
                    new.wei <- wei_nei_loc[which(!is.na(vts.cut[,t]))]
                    if (!lead_lag_weights){ # do the re-normalisation only if not lead-lag weights
                      new.wei <- new.wei/sum(new.wei) # re-normalise after removing weights for NA
                    }
                    sub.val <- vts.cut[which(!is.na(vts.cut[,t])),t] %*% new.wei # create from values of multiple time series without NAs at time t, a single value at time t by the weighted sum
                    vts.cut[which(is.na(vts.cut[,t])),t] <- sub.val # fill in NA with weighted sum of time series at time t without NAs
                  }
                }
              }
              
              # ASSIGN TO DESIGN MAT THE SUM OF WEIGHTED TIME SERIES FOR LAG
              dmat[((predt * (al - 1) + 1):(predt * al)), deltaLoc] <- as.vector(t(vts.cut) %*% wei_nei_loc)
              #wei_mat_2[[stag]][nei_loc,al] <- wei_nei_loc 
            }
            else if (length(nei_loc)==1){
              
              wei_nei_loc <- rep(1/length(nei_loc),length(nei_loc))
              wei_mat[[stag]][nei_loc,al] <- wei_nei_loc
              
              vts.cut <- unlist(ts_data[nei_loc,((alphaOrder + 1 - ilag):(predt + (alphaOrder - ilag)))])
              #and if this is missing at any time point, set to zero (cannot fill in with mean value of neighbours as previously as no other neighbours)
              vts.cut[is.na(vts.cut)] <- 0
              dmat[((predt * (al - 1) + 1):(predt * al)), deltaLoc] <- as.vector(t(vts.cut) * wei_nei_loc)
            }
            else{ # if there are no neighbours
              dmat[((predt * (al - 1) + 1):(predt * al)), deltaLoc] <- 0
            }
          } # end for stage
        } # end if 
      } # end for ilag
    } # END FOR LOOP OVER EDGES&NODES
  }
  out <- list(dmat=dmat,wei_mat=wei_mat,data_loc_mat=data_loc_mat)
  return(out)
}

gnar_x_fit <- function(ts_data,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,alphaOrder,betaOrder,gammaOrder,deltaOrder,
                       net,lead_lag_mat,globalalpha=TRUE,lead_lag_weights=FALSE,pay_now=TRUE,lag_0_sep=FALSE,nodes_only=FALSE){
  # ts_data: matrix with rows edges and columns time
  # nnodes: number of nodes
  # nedges: number of edges
  # n_edges_nodes: total number of rows of ts_data==nnodes+nedges
  # data_edges: dataframe with edge list, column one node from which edge starts (character type) column two node edge ends (character type)
  # data_nodes: vector of nodes
  # alphaOrder: max lag for autoregressive
  # betaOrder: vector of length as number of lags alphaOrder considered. Each entry give max stage neighbour to be considered for each lag
  # gammaOrder: max lag for cross-autoregressive
  # deltaOrder: vector of length as number of lags gammaOrder considered. Each entry give max stage neighbour to be considered for each cross-lag
  # net: network resulting from edge list data_edges
  # lead_lag_mat: matrix with rows and col corresponding to edge time series, entries showing leadingness of each time series over the other
  # globalalpha: if True, non edge-specific alpha parameters, but only lag specific alpha param
  # lead_lag_weights: whether lead-lag matrix used as weights or not
  
  # Returns: list with 1st element fitted model, 2nd element response vec, 3rd element design matrix, 4th element weight matrix used
  
  des <- design_mat(ts_data,nnodes,nedges,n_edges_nodes,data_edges,data_nodes,alphaOrder,betaOrder,gammaOrder,deltaOrder,net,
                    lead_lag_mat,globalalpha=globalalpha,lead_lag_weights=lead_lag_weights,pay_now=pay_now)
  if (nodes_only){
    predt <- ncol(ts_data)-alphaOrder # lenght of time series observations considered when removing max lag
    
    dmat <- tail(des$dmat,n=predt*nnodes)
    yvec <- tail(response_vec(ts_data = ts_data,alphaOrder = alphaOrder),n=predt*nnodes)
  }else{
    dmat <- des$dmat
    yvec <- response_vec(ts_data = ts_data,alphaOrder = alphaOrder)
  }
  # check if NA in response to remove, and respectively remove this rows from design mat
  if(sum(is.na(yvec))>0){
    yvec2 <- yvec[!is.na(yvec)]
    dmat2 <- dmat[!is.na(yvec),]
    if (lag_0_sep==FALSE){
      modNoIntercept <- lm(yvec2~dmat2+0)
      out <- list(mod=modNoIntercept, y=yvec2, dd=dmat2, wei_mat=des$wei_mat, data_loc_mat=des$data_loc_mat)
    }else{
      # fix design mat for two separate models
      dmat_lag <- dmat2[,-which(grepl("0", colnames(dmat2), fixed = TRUE))]
      dmat_lag0 <- dmat2[,which(grepl("0", colnames(dmat2), fixed = TRUE))]
      dmat_lag0 <- tail(dmat_lag0,n=(ncol(ts_data)-alphaOrder)*nnodes) # n: number of entries of y response that are non-zero
      # fix response vec for two separate models
      yvec_lag0 <- tail(yvec2,n=(ncol(ts_data)-alphaOrder)*nnodes)
      # fit two models
      modNoIntercept_lag0 <- lm(yvec_lag0~dmat_lag0+0)
      modNoIntercept <- lm(yvec2~dmat_lag+0)
      out <- list(modlag0=modNoIntercept_lag0,mod=modNoIntercept, y=yvec2, dd=dmat2, wei_mat=des$wei_mat, data_loc_mat=des$data_loc_mat)
    }
  }else{
    if (lag_0_sep==FALSE){
      modNoIntercept <- lm(yvec~dmat+0)
      out <- list(mod=modNoIntercept, y=yvec, dd=dmat, wei_mat=des$wei_mat, data_loc_mat=des$data_loc_mat)
    }else{
      # fix design mat for two separate models
      dmat_lag <- dmat[,-which(grepl("0", colnames(dmat), fixed = TRUE))]
      dmat_lag0 <- dmat[,which(grepl("0", colnames(dmat), fixed = TRUE))]
      dmat_lag0 <- tail(dmat_lag0,n=(ncol(ts_data)-alphaOrder)*nnodes) # n: number of entries of y response that are non-zero
      # fix response vec for two separate models
      yvec_lag0 <- tail(yvec,n=(ncol(ts_data)-alphaOrder)*nnodes)
      # fit two models
      modNoIntercept_lag0 <- lm(yvec_lag0~dmat_lag0+0)
      modNoIntercept <- lm(yvec~dmat_lag+0)
      out <- list(modlag0=modNoIntercept_lag0,mod=modNoIntercept, y=yvec, dd=dmat, wei_mat=des$wei_mat, data_loc_mat=des$data_loc_mat)
    }
  }
  return(out)
}


edge_matching <- function(edge_set, neig_edge_set){
  # edge_set: two column set of edges for data set
  # neig_edge_set: igraph object, neighbouring edges of some specified edge
  # Returns: list with indices/location of edges in original data/data edges of ts
  ids <- as_ids(neig_edge_set)
  split_ids <- strsplit(ids,split='|', fixed=TRUE)
  ids_loc <- sapply(split_ids, function(x) which((edge_set[,1]==x[1]) & (edge_set[,2]==x[2])))
  return(ids_loc)
}


edge_neighbors <- function(net, maxstage, stage_nodes){
  # net: igraph object 
  # maxstage: max stage of neighbour edges to be considered
  # stage_nodes: list of nodes corresponding to edge for which we seek neighbour edges
  # Returns: list of lists, with r element containing neighbour edges for stage r
  if (!is_connected(net)){ # if graph is not connected, identify connected component involving stage_nodes to find r-stage neighb
    net_comp <- decompose(net)
    for (i in 1:length(net_comp)){ # find subgraph (component) involving stage_nodes of interest
      if ((stage_nodes[1] %in% vertex_attr(net_comp[[i]])$name) & (stage_nodes[2] %in% vertex_attr(net_comp[[i]])$name)){
        net <- net_comp[[i]]
        break
      }
    }
  }
  neighb_edges <- vector(mode = "list", length = maxstage)
  net <- net - edge(paste(stage_nodes[1],"|",stage_nodes[2],sep = ""))
  if (length(E(net))!=0){
    for (stag in 1:maxstage){
      neighb_edges_loc=list()
      for (node in stage_nodes){
        if (length(neighb_edges_loc)==0){
          neighb_edges_loc<-incident(net, node, mode = "all")#all, out, in
        }else{
          neighb_edges_loc<-igraph::union(neighb_edges_loc,incident(net, node, mode = "all"))#all, out, in
        }
      }
      neighb_edges[[stag]] <- neighb_edges_loc
      # create new list of nodes to get their incident edges for stages>1
      stage_nodes <- setdiff(as.vector(ends(net,neighb_edges_loc)),stage_nodes)
      # remove previous neighbour edge from previous stage from graph, not be included in next set of edges
      net <- net - neighb_edges_loc
      if (length(E(net))==0) break
    }
    return(neighb_edges)
  }else{
    neighb_edges <- NULL
    return(neighb_edges)
  }
}

node_neighbors <- function(net, node, max.stage){
  tot.nei <- vector(mode="list", length=max.stage)
  if(!is.null(neighbors(net,node,mode="all"))){
    tot.nei[[1]] <- unique(neighbors(net,node,mode="all"))
  }
  
  if(max.stage>1){
    for(stag in 2:max.stage){
      if(!is.null(tot.nei[[stag-1]][1])){
        if(!is.na(tot.nei[[stag-1]][1])){
          tmp.nei <- NULL
          for(nod in 1:length(tot.nei[[stag-1]])){
            #get neighbours and weights
            tmp.nei <- c(tmp.nei, neighbors(net,tot.nei[[stag-1]][nod],mode="all"))
          }
          
          #remove node from list
          if(node%in%tmp.nei){
            posrem <- which(tmp.nei==node)
            tmp.nei <- tmp.nei[-posrem]
          }
          
          #remove nodes in previous lists
          if(sum(tmp.nei%in%unlist(tot.nei))>0){
            posrem <- which(tmp.nei%in%unlist(tot.nei))
            tmp.nei <- tmp.nei[-posrem]
          }
          
          #find minimum where have different paths to new node
          if(length(unique(tmp.nei))==0){
            break
          }else {
            tot.nei[[stag]] <- unique(tmp.nei)
          }
          
        }
      }
    }
  }
  return(tot.nei)
}



gnar_x_predict <- function(fit_mod,ts_data,alphaOrder,betaOrder,gammaOrder,deltaOrder,n_edges_nodes,nnodes,nedges,npred, wei_mat,data_loc_mat,
                           pay_now=TRUE, pay_data_now=NULL, set.noise=NULL, allcoefs=FALSE,globalalpha=TRUE){
  # fit_mod: fitted model
  # ts_data: data used to train the model
  # alphaOrder: max lag
  # betaOrder: vector of length as number of lags alphaOrder considered. Each entry give max stage neighbour to be considered for each lag
  # nedges: number of edges
  # npred: number of time stamps to be predicted
  # wei_mat: list with i element matrix with each column weights for each edge, at stage i
  if(!is.null(set.noise)){
    sig <- set.noise
  }else{
    sig <- sigma(fit_mod$mod)
  }
  if(!allcoefs){
    # keep only coefficients that are stat significant
    nas <- is.na(fit_mod$mod$coefficients)
    pvs <- summary(fit_mod$mod)$coefficients[,4] < 0.05
    vals <- rep(0, length(pvs))
    vals[pvs] <- summary(fit_mod$mod)$coefficients[pvs,1]
    coefvec <- rep(0, length(nas))
    coefvec[(!nas)] <- vals
  }else{
    # keep all coeff
    coefvec <- fit_mod$mod$coefficients
    coefvec[is.na(coefvec)] <- 0
  }
  
  if(globalalpha){
    # global alpha has one alpha per time lag
    # alphaout: list with i element the list of hat(alphas) for lag i
    # betaout: list with i element the hat(betas) for lag i
    if (pay_now==TRUE){
      alphaout <-  vector(mode="list", length=alphaOrder)
      betaout <- as.list(rep(0,length=alphaOrder))
      gammaout <-  vector(mode="list", length=(gammaOrder+1))
      deltaout <- as.list(rep(0,length=(gammaOrder+1)))
      
      count <- 1
      for (ilag in 0:max(alphaOrder,gammaOrder)){
        if (ilag==0){
          gammaout[[ilag+1]] <- coefvec[count]
          if(deltaOrder[ilag+1]>0){
            deltaout[[ilag+1]] <- coefvec[(count+1):(count+deltaOrder[ilag+1])]
          }
          count <- count + deltaOrder[ilag+1] + 1
        }else{
          alphaout[[ilag]] <- coefvec[count]
          if(betaOrder[ilag]>0){
            betaout[[ilag]] <- coefvec[(count+1):(count+betaOrder[ilag])]
          }
          gammaout[[ilag+1]] <- coefvec[count+betaOrder[ilag]+1]
          if(deltaOrder[ilag+1]>0){
            deltaout[[ilag+1]] <- coefvec[(count+betaOrder[ilag]+1+1):(count+betaOrder[ilag]+1+deltaOrder[ilag+1])]
          }
          count <- count+betaOrder[ilag]+1+deltaOrder[ilag+1] + 1
        }
      }
    }else{
      alphaout <-  vector(mode="list", length=alphaOrder)
      betaout <- as.list(rep(0,length=alphaOrder))
      gammaout <-  vector(mode="list", length=(gammaOrder))
      deltaout <- as.list(rep(0,length=(gammaOrder)))
      
      count <- 1
      for (ilag in 1:max(alphaOrder,gammaOrder)){
        alphaout[[ilag]] <- coefvec[count]
        if(betaOrder[ilag]>0){
          betaout[[ilag]] <- coefvec[(count+1):(count+betaOrder[ilag])]
        }
        gammaout[[ilag]] <- coefvec[count+betaOrder[ilag]+1]
        if(deltaOrder[ilag]>0){
          deltaout[[ilag]] <- coefvec[(count+betaOrder[ilag]+1+1):(count+betaOrder[ilag]+1+deltaOrder[ilag])]
        }
        count <- count+betaOrder[ilag]+1+deltaOrder[ilag] + 1
      }
    }
  }
  # else{
  #   #multiple alphas per time lag
  #   alphaout <-  vector(mode="list", length=alphaOrder)
  #   betaout <- as.list(rep(0,length=alphaOrder))
  #   count <- 1
  #   for(ilag in 1:alphaOrder){
  #     alphaout[[ilag]] <- coefvec[count:(count+nedges-1)]
  #     if(betaOrder[ilag]>0){
  #       betaout[[ilag]] <- coefvec[(count+nedges):(count+nedges+betaOrder[ilag]-1)]
  #     }
  #     count <- count + nedges + betaOrder[ilag]
  #   }
  # }
  
  #print(c("alphaout",alphaout))
  #print(c("betaout",betaout))
  max.nei <- max(betaOrder)
  
  # part of data used for training, to be needed for prediction of t+1
  xx.init <- ts_data[,(ncol(ts_data)-alphaOrder+1):ncol(ts_data)]
  
  # matrix with first columns data for training from xx.init and last columns the predictions
  xx.gen <- matrix(0, ncol=npred+alphaOrder, nrow=n_edges_nodes)
  xx.gen[,1:alphaOrder] <- as.matrix(xx.init)
  if (pay_now==TRUE){
    xx.gen[(1:nedges),(alphaOrder+1):(alphaOrder+npred)] <- pay_data_now
    
    # fill in the last columns (predict) of xx.gen with predicted values using xx.init
    for(tpred in (alphaOrder+1):(npred+alphaOrder)){
      
      # lags of same edge
      for(ilag in 1:alphaOrder){
        if(ilag==1){
          time.forecast <- alphaout[[ilag]]*xx.gen[((nedges+1):n_edges_nodes),(tpred-ilag)]
          #print(time.forecast)
          #print(c("time forecast: ",time.forecast))
        }else{
          tmpi <- alphaout[[ilag]]*xx.gen[((nedges+1):n_edges_nodes),(tpred-ilag)]
          time.forecast <- time.forecast + tmpi
          #print(c("time forecast: ",time.forecast))
        }
      }
      #print(c("time forecast: ",time.forecast))
      
      # for beta parameters neighbors
      if (sum(betaOrder)>0){
        nei.forecast.nodes  <- 0
        for(ilag in 1:alphaOrder){
          bb <- length(betaout[[ilag]])
          if(bb>0){
            for(dd in 1:bb){
              #nei.forecast.edges <- nei.forecast.edges + betaout[[ilag]][dd]*(xx.gen[1:nedges,(tpred-ilag)]%*%wei_mat[[dd]][1:nedges,1:nedges])
              nei.forecast.nodes <- nei.forecast.nodes + betaout[[ilag]][dd]*(xx.gen[(nedges+1):n_edges_nodes,(tpred-ilag)]%*%wei_mat[[dd]][(nedges+1):n_edges_nodes,(nedges+1):n_edges_nodes])
              #print(c("nei forecast: ",nei.forecast))
            }
          }
        }
      }else{
        nei.forecast.nodes  <- rep(0,nnodes)
      }
      
      # for gamma parameters 
      for(ilag in 0:gammaOrder){
        if(ilag==0){
          time.forecast.cross <- gammaout[[ilag+1]]*(xx.gen[(1:nedges),(tpred-ilag)]%*%data_loc_mat[(1:nedges),((nedges+1):n_edges_nodes)])
          #print(time.forecast.cross)
          #print(c("time forecast: ",time.forecast))
        }else{
          tmpi <- gammaout[[ilag+1]]*(xx.gen[(1:nedges),(tpred-ilag)]%*%data_loc_mat[(1:nedges),((nedges+1):n_edges_nodes)])
          time.forecast.cross <- time.forecast.cross + tmpi
          #print(time.forecast.cross)
          #print(c("time forecast: ",time.forecast))
        }
      }
      
      # for delta parameters neighbors
      if (sum(deltaOrder)>0){
        nei.forecast.nodes.cross <- 0
        for(ilag in 0:gammaOrder){
          bb <- length(deltaout[[ilag+1]])
          if(bb>0){
            for(dd in 1:bb){
              nei.forecast.nodes.cross <- nei.forecast.nodes.cross + deltaout[[ilag+1]][dd]*(xx.gen[1:nedges,(tpred-ilag)]%*%wei_mat[[dd]][1:nedges,(nedges+1):n_edges_nodes])
              # print(c("nei forecast: ",nei.forecast))
            }
          }
        }
      }else{
        nei.forecast.nodes.cross <- rep(0,nnodes)
      }
      
      xx.gen[(nedges+1):n_edges_nodes,tpred] <- time.forecast+as.vector(nei.forecast.nodes)+as.vector(time.forecast.cross)+as.vector(nei.forecast.nodes.cross)
    }
  }else{
    # fill in the last columns (predict) of xx.gen with predicted values using xx.init
    for(tpred in (alphaOrder+1):(npred+alphaOrder)){
      
      # lags of same edge
      for(ilag in 1:alphaOrder){
        if(ilag==1){
          time.forecast <- alphaout[[ilag]]*xx.gen[,(tpred-ilag)]
          #print(time.forecast)
          #print(c("time forecast: ",time.forecast))
        }else{
          tmpi <- alphaout[[ilag]]*xx.gen[,(tpred-ilag)]
          time.forecast <- time.forecast + tmpi
          #print(c("time forecast: ",time.forecast))
        }
      }
      #print(c("time forecast: ",time.forecast))
      
      # for beta parameters neighbors
      nei.forecast.nodes  <- nei.forecast.edges <- 0
      if (sum(betaOrder)>0){
        for(ilag in 1:alphaOrder){
          bb <- length(betaout[[ilag]])
          if(bb>0){
            for(dd in 1:bb){
              nei.forecast.edges <- nei.forecast.edges + betaout[[ilag]][dd]*(xx.gen[1:nedges,(tpred-ilag)]%*%wei_mat[[dd]][1:nedges,1:nedges])
              nei.forecast.nodes <- nei.forecast.nodes + betaout[[ilag]][dd]*(xx.gen[(nedges+1):n_edges_nodes,(tpred-ilag)]%*%wei_mat[[dd]][(nedges+1):n_edges_nodes,(nedges+1):n_edges_nodes])
              #print(c("nei forecast: ",nei.forecast))
            }
          }
        }
        nei.forecast <- cbind(nei.forecast.edges,nei.forecast.nodes)
      }
      
      # for gamma parameters 
      for(ilag in 1:gammaOrder){
        if(ilag==1){
          time.forecast.cross <- gammaout[[ilag]]*(xx.gen[,(tpred-ilag)]%*%data_loc_mat)
          #print(time.forecast.cross)
          #print(c("time forecast: ",time.forecast))
        }else{
          tmpi <- gammaout[[ilag]]*(xx.gen[,(tpred-ilag)]%*%data_loc_mat)
          time.forecast.cross <- time.forecast.cross + tmpi
          #print(time.forecast.cross)
          #print(c("time forecast: ",time.forecast))
        }
      }
      
      # for delta parameters neighbors
      nei.forecast.nodes.cross <- nei.forecast.edges.cross <- 0
      if (sum(deltaOrder)>0){
        for(ilag in 1:gammaOrder){
          bb <- length(deltaout[[ilag]])
          if(bb>0){
            for(dd in 1:bb){
              nei.forecast.edges.cross <- nei.forecast.edges.cross + deltaout[[ilag]][dd]*(xx.gen[(nedges+1):n_edges_nodes,(tpred-ilag)]%*%wei_mat[[dd]][(nedges+1):n_edges_nodes,1:nedges])
              nei.forecast.nodes.cross <- nei.forecast.nodes.cross + deltaout[[ilag]][dd]*(xx.gen[1:nedges,(tpred-ilag)]%*%wei_mat[[dd]][1:nedges,(nedges+1):n_edges_nodes])
              # print(c("nei forecast: ",nei.forecast))
            }
          }
        }
        nei.forecast.cross <- cbind(nei.forecast.edges.cross,nei.forecast.nodes.cross)
      }
      
      if (sum(deltaOrder)==0 & sum(betaOrder)==0){
        xx.gen[,tpred] <- time.forecast+as.matrix(time.forecast.cross)
      }else{
        xx.gen[,tpred] <- time.forecast+as.matrix(nei.forecast)+as.matrix(time.forecast.cross)+as.matrix(nei.forecast.cross)
      }
    }
  }
  
  return(xx.gen)
}

newdat <- function(fit_mod,ts_data,alphaOrder,betaOrder,gammaOrder,deltaOrder,n_edges_nodes,nnodes,nedges,npred, wei_mat,data_loc_mat){
  
  xx.init <- ts_data[,(ncol(ts_data)-alphaOrder+1):ncol(ts_data)]
  # matrix with first columns data for training from xx.init and last columns the predictions
  xx.gen <- matrix(0, ncol=npred+alphaOrder, nrow=n_edges_nodes)
  xx.gen[,1:alphaOrder] <- as.matrix(xx.init)
  
  # create data frame for new predictors to forecast next time stamp
  xx.new <- matrix(0,ncol = length(fit_mod$mod$coefficients),nrow = n_edges_nodes)
  colnames(xx.new) <- colnames(fit_mod$dd)#names(fit_mod$mod$coefficients)
  
  # fill in the last columns (predict) of xx.gen with predicted values using xx.init
  for(tpred in (alphaOrder+1):(npred+alphaOrder)){
    
    # lags of same edge
    for(ilag in 1:alphaOrder){
      xx.new[,paste("alpha",ilag,sep = "")] <- xx.gen[,(tpred-ilag)]
    }
    
    # for beta parameters neighbors
    if (sum(betaOrder)>0){
      for(ilag in 1:alphaOrder){
        bb <- length(which(grepl(paste("dmatbeta",ilag,sep = ""),names(fit_mod$mod$coefficients),fixed = TRUE)))#length(betaout[[ilag]])
        if(bb>0){
          for(dd in 1:bb){
            nei.forecast.edges <- xx.gen[1:nedges,(tpred-ilag)]%*%wei_mat[[dd]][1:nedges,1:nedges]
            nei.forecast.nodes <- xx.gen[(nedges+1):n_edges_nodes,(tpred-ilag)]%*%wei_mat[[dd]][(nedges+1):n_edges_nodes,(nedges+1):n_edges_nodes]
            xx.new[,paste("beta",ilag,".stage",dd,sep = "")] <- as.vector(cbind(nei.forecast.edges,nei.forecast.nodes))
          }
        }
      }
    }
    
    # for gamma parameters
    for(ilag in 1:gammaOrder){
      xx.new[,paste("gamma",ilag,sep = "")] <- as.vector(xx.gen[,(tpred-ilag)]%*%data_loc_mat)
    }
    
    # for delta parameters neighbors
    if (sum(deltaOrder)>0){
      for(ilag in 1:gammaOrder){
        bb <- length(which(grepl(paste("dmatdelta",ilag,sep = ""),names(fit_mod$mod$coefficients),fixed = TRUE)))#length(deltaout[[ilag]])
        if(bb>0){
          for(dd in 1:bb){
            nei.forecast.edges.cross <- xx.gen[(nedges+1):n_edges_nodes,(tpred-ilag)]%*%wei_mat[[dd]][(nedges+1):n_edges_nodes,1:nedges]
            nei.forecast.nodes.cross <- xx.gen[1:nedges,(tpred-ilag)]%*%wei_mat[[dd]][1:nedges,(nedges+1):n_edges_nodes]
            xx.new[,paste("delta",ilag,".stage",dd,sep = "")] <- as.vector(cbind(nei.forecast.edges.cross,nei.forecast.nodes.cross))
          }
        }
      }
    }
  }
  return(xx.new)
}


gnar_x_sim <- function(n=200, net, alphaParams, betaParams, gammaParams, deltaParams, sigma=1, meann=0, nedges,nnodes,n_edges_nodes,
                       data_edges,data_nodes,noise="n",vt=3,rcor=0.5,pay_now=TRUE){
  # n: number of time stamps
  # net: network on whose edges are simulated
  # alphaParams: alpha parameters
  # betaParams: beta parameters
  # nedges: number of edges
  # data_edges: dataframe with edge list, column one node from which edge starts (character type) column two node edge ends (character type)
  # noise: categorical ("n","t")
  # vt: degrees of freedom for t distribution
  
  max.nei <- max(unlist(lapply(betaParams, length)),unlist(lapply(deltaParams, length)))
  
  #create weight matrices for neighbours
  nei.mats <- nei_wei_mat(net,data_edges,data_nodes,max.nei,nedges,nnodes,n_edges_nodes)
  # create matrix with associated edges/nodes for nodes/edges respectively
  data_loc_mat <- node_edge_assoc(n_edges_nodes,nedges,nnodes,data_nodes,data_edges)
  
  #seed the process from normal dist with mean 0 and sigma as given
  #do this for as many time points as needed for alpha
  if (pay_now==TRUE){
    lags <- max(length(alphaParams),length(gammaParams)-1)
  }else{
    lags <- max(length(alphaParams),length(gammaParams))
  }
  
  # initialise time series to use lags according to max lag size, to produce next time steps
  xx.init <- matrix(rnorm(n_edges_nodes*lags, mean=0, sd=sigma), nrow=n_edges_nodes, ncol=lags)
  
  xx.gen <- matrix(0, ncol=n+50, nrow=n_edges_nodes) # rows edges, columns time series 
  xx.gen[,1:lags] <- xx.init
  
  if (pay_now==FALSE){
    # create matrix of correlated errors across time (dependent Norm across time with correlation rcor)
    if (noise=="corr_t"){
      corr.err <- rnorm_multi(n_edges_nodes,mu = rep(0,(n+50-lags)), sd = rep(1,(n+50-lags)), r = rcor)
    }
    ind.tt <- 1 # index for column of error for correlated errors
    
    for(tt in (lags+1):(n+50)){
      
      # for alpha parameters
      for(ilag in 1:lags){
        if(ilag==1){
          time.forecast <- alphaParams[[ilag]]*xx.gen[,(tt-ilag)]
        }else{
          tmpi <- alphaParams[[ilag]]*xx.gen[,(tt-ilag)]
          time.forecast <- time.forecast + tmpi
        }
      }
      
      # for beta parameters neighbors
      nei.forecast.nodes <- nei.forecast.edges <- 0
      for(ilag in 1:lags){
        bb <- length(betaParams[[ilag]])
        if(bb>0){
          for(dd in 1:bb){
            nei.forecast.edges <- nei.forecast.edges + betaParams[[ilag]][dd]*(xx.gen[1:nedges,(tt-ilag)]%*%nei.mats[[dd]][1:nedges,1:nedges])
            nei.forecast.nodes <- nei.forecast.nodes + betaParams[[ilag]][dd]*(xx.gen[(nedges+1):n_edges_nodes,(tt-ilag)]%*%nei.mats[[dd]][(nedges+1):n_edges_nodes,(nedges+1):n_edges_nodes])
            #print(c("nei forecast: ",nei.forecast))
          }
        }
      }
      nei.forecast <- cbind(nei.forecast.edges,nei.forecast.nodes)
      
      # for gamma parameters 
      for(ilag in 1:lags){
        if(ilag==1){
          time.forecast.cross <- gammaParams[[ilag]]*(xx.gen[,(tt-ilag)]%*%data_loc_mat)
          #print(c("time forecast: ",time.forecast))
        }else{
          tmpi <- gammaParams[[ilag]]*(xx.gen[,(tt-ilag)]%*%data_loc_mat)
          time.forecast.cross <- time.forecast.cross + tmpi
          #print(c("time forecast: ",time.forecast))
        }
      }
      
      # for delta parameters neighbors
      nei.forecast.nodes.cross <- nei.forecast.edges.cross <- 0
      for(ilag in 1:lags){
        bb <- length(deltaParams[[ilag]])
        if(bb>0){
          for(dd in 1:bb){
            nei.forecast.edges.cross <- nei.forecast.edges.cross + deltaParams[[ilag]][dd]*(xx.gen[(nedges+1):n_edges_nodes,(tt-ilag)]%*%nei.mats[[dd]][(nedges+1):n_edges_nodes,1:nedges])
            nei.forecast.nodes.cross <- nei.forecast.nodes.cross + deltaParams[[ilag]][dd]*(xx.gen[1:nedges,(tt-ilag)]%*%nei.mats[[dd]][1:nedges,(nedges+1):n_edges_nodes])
            #print(c("nei forecast: ",nei.forecast))
          }
        }
      }
      nei.forecast.cross <- cbind(nei.forecast.edges.cross,nei.forecast.nodes.cross)
      
      if (noise=="n"){
        xx.gen[,tt] <- time.forecast+as.matrix(nei.forecast)+as.matrix(time.forecast.cross)+as.matrix(nei.forecast.cross)+rnorm(n_edges_nodes, mean=meann, sd=sigma)
      }else if (noise=="t"){
        xx.gen[,tt] <- time.forecast+as.matrix(nei.forecast)+as.matrix(time.forecast.cross)+as.matrix(nei.forecast.cross)+rt(n_edges_nodes, df=vt)
      }else if (noise=="corr_t"){
        xx.gen[,tt] <- time.forecast+as.matrix(nei.forecast)+as.matrix(time.forecast.cross)+as.matrix(nei.forecast.cross)+corr.err[,ind.tt]
      }
      ind.tt <- ind.tt+1
    }
  }else{
    # create matrix of correlated errors across time (dependent Norm across time with correlation rcor)
    if (noise=="corr_t"){
      corr.err.edg <- rnorm_multi(nedges,mu = rep(0,(n+50-lags)), sd = rep(1,(n+50-lags)), r = rcor)
      corr.err.nod <- rnorm_multi(nnodes,mu = rep(0,(n+50-lags)), sd = rep(1,(n+50-lags)), r = rcor)
    }
    ind.tt <- 1 # index for column of error for correlated errors
    
    for(tt in (lags+1):(n+50)){
      ########################
      # FIRST OBTAIN tt FOR EDGES
      ########################
      # for alpha parameters
      for(ilag in 1:lags){
        if(ilag==1){
          time.forecast <- alphaParams[[ilag]]*xx.gen[(1:nedges),(tt-ilag)]
        }else{
          tmpi <- alphaParams[[ilag]]*xx.gen[(1:nedges),(tt-ilag)]
          time.forecast <- time.forecast + tmpi
        }
      }
      
      # for beta parameters neighbors
      nei.forecast.edges <- 0
      for(ilag in 1:lags){
        bb <- length(betaParams[[ilag]])
        if(bb>0){
          for(dd in 1:bb){
            nei.forecast.edges <- nei.forecast.edges + betaParams[[ilag]][dd]*(xx.gen[1:nedges,(tt-ilag)]%*%nei.mats[[dd]][1:nedges,1:nedges])
            #nei.forecast.nodes <- nei.forecast.nodes + betaParams[[ilag]][dd]*(xx.gen[(nedges+1):n_edges_nodes,(tt-ilag)]%*%nei.mats[[dd]][(nedges+1):n_edges_nodes,(nedges+1):n_edges_nodes])
            #print(c("nei forecast: ",nei.forecast))
          }
        }
      }
      #nei.forecast <- cbind(nei.forecast.edges,nei.forecast.nodes)
      
      # for gamma parameters 
      for(ilag in 1:lags){
        if(ilag==1){
          time.forecast.cross <- gammaParams[[ilag+1]]*(xx.gen[((nedges+1):n_edges_nodes),(tt-ilag)]%*%data_loc_mat[((nedges+1):n_edges_nodes),(1:nedges)])
          #print(c("time forecast: ",time.forecast))
        }else{
          tmpi <- gammaParams[[ilag+1]]*(xx.gen[((nedges+1):n_edges_nodes),(tt-ilag)]%*%data_loc_mat[((nedges+1):n_edges_nodes),(1:nedges)])
          time.forecast.cross <- time.forecast.cross + tmpi
          #print(c("time forecast: ",time.forecast))
        }
      }
      
      # for delta parameters neighbors
      nei.forecast.edges.cross <- 0
      for(ilag in 1:lags){
        bb <- length(deltaParams[[ilag+1]])
        if(bb>0){
          for(dd in 1:bb){
            nei.forecast.edges.cross <- nei.forecast.edges.cross + deltaParams[[ilag+1]][dd]*(xx.gen[(nedges+1):n_edges_nodes,(tt-ilag)]%*%nei.mats[[dd]][(nedges+1):n_edges_nodes,1:nedges])
            #nei.forecast.nodes.cross <- nei.forecast.nodes.cross + deltaParams[[ilag+1]][dd]*(xx.gen[1:nedges,(tt-ilag)]%*%nei.mats[[dd]][1:nedges,(nedges+1):n_edges_nodes])
            #print(c("nei forecast: ",nei.forecast))
          }
        }
      }
      #nei.forecast.cross <- cbind(nei.forecast.edges.cross,nei.forecast.nodes.cross)
      
      xx.gen[1:nedges,tt] <- time.forecast+as.matrix(nei.forecast.edges)+as.matrix(time.forecast.cross)+as.matrix(nei.forecast.edges.cross)
      
      if (noise=="n"){
        xx.gen[1:nedges,tt] <- xx.gen[1:nedges,tt]+rnorm(nedges, mean=meann, sd=sigma)
      }else if (noise=="t"){
        xx.gen[1:nedges,tt] <- xx.gen[1:nedges,tt]+rt(nedges, df=vt)
      }else if (noise=="corr_t"){
        xx.gen[1:nedges,tt] <- xx.gen[1:nedges,tt]+corr.err.edg[,ind.tt]
      }
      
      ########################
      # SECOND OBTAIN tt FOR NODES
      ########################
      
      # for alpha parameters
      for(ilag in 1:lags){
        if(ilag==1){
          time.forecast <- alphaParams[[ilag]]*xx.gen[((nedges+1):n_edges_nodes),(tt-ilag)]
        }else{
          tmpi <- alphaParams[[ilag]]*xx.gen[((nedges+1):n_edges_nodes),(tt-ilag)]
          time.forecast <- time.forecast + tmpi
        }
      }
      
      # for beta parameters neighbors
      nei.forecast.nodes <- 0
      for(ilag in 1:lags){
        bb <- length(betaParams[[ilag]])
        if(bb>0){
          for(dd in 1:bb){
            nei.forecast.nodes <- nei.forecast.nodes + betaParams[[ilag]][dd]*(xx.gen[(nedges+1):n_edges_nodes,(tt-ilag)]%*%nei.mats[[dd]][(nedges+1):n_edges_nodes,(nedges+1):n_edges_nodes])
          }
        }
      }
      
      # for gamma parameters 
      for(ilag in 0:lags){
        if(ilag==0){
          time.forecast.cross <- gammaParams[[ilag+1]]*(xx.gen[(1:nedges),(tt-ilag)]%*%data_loc_mat[(1:nedges),((nedges+1):n_edges_nodes)])
          #print(c("time forecast: ",time.forecast))
        }else{
          tmpi <- gammaParams[[ilag+1]]*(xx.gen[(1:nedges),(tt-ilag)]%*%data_loc_mat[(1:nedges),((nedges+1):n_edges_nodes)])
          time.forecast.cross <- time.forecast.cross + tmpi
          #print(c("time forecast: ",time.forecast))
        }
      }
      
      # for delta parameters neighbors
      nei.forecast.nodes.cross <- 0
      for(ilag in 0:lags){
        bb <- length(deltaParams[[ilag+1]])
        if(bb>0){
          for(dd in 1:bb){
            nei.forecast.nodes.cross <- nei.forecast.nodes.cross + deltaParams[[ilag+1]][dd]*(xx.gen[1:nedges,(tt-ilag)]%*%nei.mats[[dd]][1:nedges,(nedges+1):n_edges_nodes])
          }
        }
      }
      
      xx.gen[(nedges+1):n_edges_nodes,tt] <- time.forecast+as.matrix(nei.forecast.nodes)+as.matrix(time.forecast.cross)+as.matrix(nei.forecast.nodes.cross)
      
      if (noise=="n"){
        xx.gen[(nedges+1):n_edges_nodes,tt] <- xx.gen[(nedges+1):n_edges_nodes,tt]+rnorm(nnodes, mean=meann, sd=sigma)
      }else if (noise=="t"){
        xx.gen[(nedges+1):n_edges_nodes,tt] <- xx.gen[(nedges+1):n_edges_nodes,tt]+rt(nnodes, df=vt)
      }else if (noise=="corr_t"){
        xx.gen[(nedges+1):n_edges_nodes,tt] <- xx.gen[(nedges+1):n_edges_nodes,tt]+corr.err.nod[,ind.tt]
      }
      ind.tt <- ind.tt+1
    }
    
  }
  return(xx.gen[,51:(n+50)])
}

nei_wei_mat <- function(net,data.edges,data.nodes,max.stage,nedges,nnodes,n_edges_nodes){
  # net: igraph network
  # data.edges: matrix with edge list for net
  # max.stage: (scalar) max r-stage neighbour edges
  # nedges: (scalar) number of edges
  # Rerurns: list of r-stage matrices, each matrix of size nedges x nedges, 
  # with each column j corresponding to neighbours edges of edge j with equal weights 
  wei_mat <- lapply(1:max.stage,function(x)matrix(0,nrow=n_edges_nodes,ncol=n_edges_nodes))
  for (edg in 1:n_edges_nodes){
    #print(c("edges",edg))
    
    # FIRST: GET NEIGHBOUR EDGES FOR EDGES, AND NEIGHBOUR NODES FOR NODES 
    # I.E. FILL IN THE 1:NEDGESx1:NEDGES AND 1:NNODES X 1:NNODES OF WEI_MAT
    if (edg <= nedges){
      nodess <- c(data.edges[edg,1],data.edges[edg,2]) # get nodes involved in edg
      nei <- edge_neighbors(net,max.stage,nodess) # for these nodess get the r-stage neighbours
    }else{
      nei <- node_neighbors(net,data.nodes[edg-nedges],max.stage)
    }
    
    # equally weight neighbour edges (mean)
    wei <- sapply(nei, length)
    wei <- lapply(wei, function(x) rep(1/x,x))
    
    
    
    for (stag in 1:max.stage){
      if (!is.null(nei[[stag]])){
        if (edg <= nedges){
          data_loc <- edge_matching(data.edges, nei[[stag]])  
        }else{
          data_loc <- which(data.nodes %in% nei[[stag]])+nedges
        }
        wei_mat[[stag]][data_loc,edg] <- wei[[stag]]
      }else{
        #wei_mat[[stag]] <- NULL
        warning("beta order too large for network, neighbour set ",stag," is empty")
      }
    }
  }
  
  # SECOND: GET NEIGHBOUR EDGES FOR NODES, AND NEIGHBOUR NODES FOR EDGES 
  # I.E. FILL IN THE 1:NEDGES x (NEDGES+1):(NEDGES+NNODES) AND (NEDGES+1):(NNODES+NEDGES) X 1:NEDGES OF WEI_MAT USING ALREADY OBTAINED WEI_MAT
  
  for (edg in 1:n_edges_nodes){
    for (stag in 1:max.stage){
      if (edg <= nedges){
        data_loc <- which(data.nodes %in% data.edges[edg,])+nedges # location of nodes relevant to edge, alternative:which(rownames(ts_data) %in% data_edges[r,])
        # data_loc for edges will always involve two nodes, thus apply below always applies
        nei_loc <- which(apply(wei_mat[[stag]][(nedges+1):n_edges_nodes,data_loc],1,sum)!=0) + nedges
        nei_loc <- setdiff(nei_loc,data_loc) # exclude the nodes that are involved in data_loc from current al
      }else{
        data_loc <- which(apply(data.edges,1,function(x) data.nodes[edg - nedges] %in% x))
        if (length(data_loc)>1){ # case where more than one edge involved for the al node
          nei_loc <- which(apply(wei_mat[[stag]][1:nedges,data_loc],1,sum)!=0) 
          nei_loc <- setdiff(nei_loc,data_loc) # exclude the edges that are involved in data_loc from current al
        }else if (length(data_loc)==1) { # if exactly one edge involved
          nei_loc <- which(wei_mat[[stag]][1:nedges,data_loc]!=0) 
        }else{ # if node is not involved in any edge, also no neighbouring edges
          nei_loc <- NULL
        }
      }
      if (!is.null(nei_loc)){ # if there are neighbours (at least one)
        wei_nei_loc <- rep(1/length(nei_loc),length(nei_loc))
        wei_mat[[stag]][nei_loc,edg] <- wei_nei_loc
      }else{
        warning("delta order too large for network, neighbour set ",stag," is empty")
      }
    }
  }
  return(wei_mat)
}

node_edge_assoc <- function(n_edges_nodes,nedges,nnodes,data_nodes,data_edges){
  data_loc_mat <- Matrix(0,n_edges_nodes,n_edges_nodes,doDiag = FALSE)
  for (r in 1:n_edges_nodes) {
    if (r <= nedges){
      nodeloc <- which(data_nodes %in% data_edges[r,])+nedges # location of nodes relevant to edge, alternative:which(rownames(ts_data) %in% data_edges[r,])
      if (length(nodeloc)==1){
        data_loc_mat[nodeloc,r] <- 1
      }else{
        data_loc_mat[nodeloc,r] <- c(0.5,0.5) 
      }
    }else{
      edgeloc <- which(apply(data_edges,1,function(x) data_nodes[r-nedges] %in% x))
      if (length(edgeloc)>0){
        edgewei <- rep(1/length(edgeloc),length(edgeloc)) # equally weight edges related to node r
        data_loc_mat[edgeloc,r] <- edgewei
      }
    }
  }
  return(data_loc_mat)
}

gnarxbic <- function(fitmod,timetrain,globalalpha=TRUE,makepd=FALSE,nodes_only=FALSE,no_nodes=NULL){
  # fitmod: fitted GNAR-x model
  # timetrain: timestamps of training sample data set
  # returns BIC information criterion
  if (nodes_only){
    resmat <- matrix(fitmod$mod$residuals,  ncol=no_nodes, byrow=FALSE)
  }else{
    resmat <- matrix(fitmod$mod$residuals,  ncol=nrow(fitmod$data_loc_mat), byrow=FALSE)
  }
  resmat[is.na(resmat)] <- 0 
  covresmat <- (t(resmat) %*% resmat)/timetrain # covariance matrix of residuals
  #stopifnot(det(covresmat) != 0)
  if (det(covresmat)<=0){
    if (makepd==FALSE){
      for (i in 2:60){
        if (det(covresmat*(2^i))>0){
          #if (det(covresmat*(i))>0){  
          #const <- i
          const <- 2^i
          tmp1 <- log(det(covresmat*const))
          tmp2 <- (ncol(fitmod$dd) * log(timetrain)) / timetrain
          return(tmp1 + tmp2 - dim(covresmat)[1]*const)
        }
      }
      tmp1 <- log(det(covresmat))
      tmp2 <- (ncol(fitmod$dd) * log(timetrain)) / timetrain
      return(tmp1 + tmp2)
    }else if (makepd==TRUE & !is.positive.definite(covresmat)){
      covresmat <- make.positive.definite(covresmat)
      tmp1 <- sum(log(eigen(covresmat,only.values = TRUE)$values))
      tmp2 <- (ncol(fitmod$dd) * log(timetrain)) / timetrain
      return(tmp1 + tmp2)
    }else if (makepd==TRUE & is.positive.definite(covresmat)){
      tmp1 <- sum(log(eigen(covresmat,only.values = TRUE)$values))
      tmp2 <- (ncol(fitmod$dd) * log(timetrain)) / timetrain
      return(tmp1 + tmp2)
    }
  }else{
    tmp1 <- log(det(covresmat))
    tmp2 <- (ncol(fitmod$dd) * log(timetrain)) / timetrain
    return(tmp1 + tmp2)
  }
  
  
  # if(globalalpha){
  #   tmp2 <- (ncol(fitmod$dd) * log(timetrain)) / timetrain
  # }else{
  #   tmp2 <- ( ncol(tmp.resid) * alphas.in + sum(betas.in) ) * log(tot.time) / tot.time
  # }
  # return(tmp1 + tmp2 - dim(covresmat)[1]*const)
}


gnarxaic <- function(fitmod,timetrain,globalalpha=TRUE){
  
  resmat <- matrix(fitmod$mod$residuals,  ncol=nrow(fitmod$data_loc_mat), byrow=FALSE)
  resmat[is.na(resmat)] <- 0 
  covresmat <- (t(resmat) %*% resmat)/timetrain # covariance matrix of residuals
  stopifnot(det(covresmat) != 0)
  tmp1 <- log(det(covresmat))
  
  if(globalalpha){
    tmp2 <- (ncol(fitmod$dd) * 2) / timetrain
  }#else{
  #  tmp2 <- ( ncol(tmp.resid) * alphas.in + sum(betas.in) ) * k / tot.time
  #}
  
  return(tmp1 + tmp2)
}
