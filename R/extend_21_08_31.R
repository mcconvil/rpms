#######################################################################################################
#
# Code to produce random forest from rpms function 
# 
#  functions: forest, predict.forest
#  
#  Last updated: 4/13/2017
#
########################################################################################################


##############################################################################################################
#                       internal function f_rpms Fast rpms for use with rpms_forest 
################################################################################################################

f_rpms<-function(rp_equ, data, weights, strata, clusters, e_equ, e_fn, l_fn, bin_size,
                 y, mX, X, vX, cat_vec, des_ind, n_cats, cat_table, modfit0){
  
  ##################################  Recursive Partitioning  ####################################################
  #          calls C++ funtions get_node and split_rpms 
  ################################################################################################################
  
  # # ----- randomize tree growth ------------------------------
  # N <- nrow(data)
  # #--randomize data set ----- 
  # b <- sample(N, size=N, replace = TRUE) #bootstap sample
  # 
  # # -----randomize variables used -------------------------
  # # only consider these variables for splitting
  # #randv=sample(seq(vX), size=sample((length(vX)-1):length(vX), size=1))
  # randv=seq(vX) #use all the variables
  # randv=randv - 1 # to make index start at 0 for C++
  # 
  # 
  # #-------end randomize ------------------------------------
  # 
  # frame <-
  #   rbind_splits(get_node(node=1, cat=as.integer(NA), vname="Root", y=as.matrix(y[b]), weights=weights[b], mxval=as.matrix(NA), s=as.matrix(0), 
  #                         modfit=modfit0), 
  #                random_split(node=1, y=y[b], mX=mX, X=X[b,], vnames=vX, cat_vec=cat_vec,
  #                             weights=weights[b], strata=strata[b], clusters=clusters[b], des_ind=des_ind,
  #                             bin_size=bin_size, randv))
  

    # -----randomize variables used -------------------------
  # only consider these variables for splitting
  #randv=sample(seq(vX), size=sample((length(vX)-1):length(vX), size=1))
#  randv=seq(vX) #use all the variables
#  randv=randv - 1 # to make index start at 0 for C++
  
  
  #-------end randomize ------------------------------------

  frame <-
    rbind_splits(get_node(node=1, cat=as.integer(NA), vname="Root", y=as.matrix(y), weights=weights, mxval=as.matrix(NA), s=as.matrix(0), 
                          modfit=modfit0), 
                 random_split2(node=1, y=y, mX=mX, X=X, vnames=vX, cat_vec=cat_vec,
                            weights=weights, strata=strata, clusters=clusters, des_ind=des_ind,
                            bin_size=bin_size))
                 # rand_split(node=1, y=y, mX=mX, X=X, vnames=vX, cat_vec=cat_vec,
                 #            weights=weights, strata=strata, clusters=clusters, des_ind=des_ind,
                 #            bin_size=bin_size, gridpts=gridpts))

  

  
  ################################################################################################################
  
  #------ Make frame nice ------- 
  frame <-  make_nice(frame) # format resulting tree frame for R 
  frame <- mark_ends(frame) # puts an 'E' after each end node on the frame

  
  #======================================= return to original the factor levels in the data ==============================
  if(n_cats>0){
    for(i in 1:n_cats){
      #---- return to original levels -------
      levels(data[,vX[which(cat_vec)[i]]]) <- unlist(cat_table[[i]]) #needed to use names apparently
    }
  }
  
  # put original level names in the frame
  cat_splits<-which(frame$cat==1)
  
  for(i in cat_splits){
    vis <- which(vX[cat_vec]==frame$var[[i]]) #variable location in cat_table
    frame$xval[i] <- list(cat_table[[vis]][unlist(frame$xval[i])]) #replace numbers with factor names
  }
  #======================================================================================================================
  
  partition<-get_partition(frame)
  row.names(partition)<-NULL
  #partition$value<-as.numeric(partition$value)
  #partition$value<-as.numeric(partition$cvar)



  
  #----- endnodes is a vector containing the node number containing that observation, for each observation in data -----------
  #---- n-vector -----
 # endnodes<-apply(covariates(partition[,"splits"], data), 1, function(x) partition$node[which(x==1)])
  endnodes<-partition$node[as.integer(covariates(partition$splits, data) %*% seq(length(partition$splits)))]

  #----- fm_ends is vector containing the row index of the frame corresponding to each partion end node
  #--- k-vector -----
  fm_ends <- match(partition$node, frame$node)
  
  N_hat <- sapply(partition$node, function(x) sum(weights[endnodes==x]))
  partition <- cbind(partition, N_hat)
 
  # y_var <- sapply(partition$node, function(x) return(sum((y[which(endnodes==x)]-frame$mean[which(frame$node==x)])^2)/(frame$n[which(frame$node==x)]-1)))
  # partition <- cbind(partition, y_var)
  # y_sig2 <- frame$loss[fm_ends]
  # y_sig2_0 <- frame$loss[1]
  
  y_sig2 <- frame$loss[fm_ends]/(partition$N_hat-1)
  y_sig2_0 <- frame$loss[1]/(sum(weights)-1)
  
  partition$y_sig2 <- y_sig2
  partition$lambda <- 1/(y_sig2 + 1)
  
  r_2 <- 1 - sum(y_sig2)/y_sig2_0

  t1 <- list(callvals=NULL, rp_equ=rp_equ, e_equ=e_equ, r_2=r_2,  frame=frame, partition=partition)
  class(t1)<-c("rpms")
  return(t1)
  
}
#############################################  End f_rpms ####################################################


###############################################################################################
#                                                                                             #
#                                            forest                                           #
#                                                                                             #
###############################################################################################
#' rpms_forest
#' 
#' @description produces a random forest using rpms to create the individual
#'              trees.
#'
#' @param rp_equ formula containing all variables for partitioning
#' @param data data.frame that includes variables used in rp_equ, e_equ, 
#'        and design information
#' @param weights formula or vector of sample weights for each observation 
#' @param strata formula or vector of strata labels
#' @param clusters formula or vector of cluster labels
#' @param e_fn string name of function to use for modeling 
#'        (only "survLm" is operational)
#' @param l_fn loss function (ignored)
#' @param bin_size numeric minimum number of observations in each node 
#' @param f_size integer specifying the number of trees in the forest        
#' @param cores integer number of cores to use in parallel if > 1 
#'        (doesn't work with Windows operating systems)
#' 
#' @return object of class "rpms"
#' 
#' @export
rpms_forest <- function(rp_equ, data, weights=~1, strata=~1, clusters=~1, 
                        e_fn="survLm", l_fn=NULL, bin_size=5, f_size=500, cores=1){

  #=============  Check paramaters and format data set ====================================================================
  #
  #check_parmas CHANGES paramater values 
  # returns list with y, X, mX, vX, des_ind, cat_vec, n_cats, cat_table
  #
  # ok <- check_params(rp_equ=rp_equ, data=data, e_equ=e_equ, weights=weights, strata=strata, clusters=clusters, 
  #                    bin_size=bin_size, random=random, gridpts=gridpts, perm_reps=perm_reps, pval=pval)
  
  if(is.null(bin_size)) bin_size <- 5
  else 
    if(bin_size<=2) stop("bin_size must be set > 1")  
  
  
  #------- check variables in equations ------
  if(length(all.vars(rp_equ))==0) 
    stop("no variable for recursive partition given in formula") 
  
  
  
  #============= format data set ================================================
  varlist <- unique(c(all.vars(rp_equ), all.vars(weights), 
                      all.vars(strata), all.vars(clusters)))  
  
  
  
  #---------- check all variables are in data set --------------
  if(!all(varlist %in% names(data))){
    e1ind <- which(!(varlist %in% names(data)))[1]
    stop(paste("Variable", varlist[e1ind], "not in dataset"))
  }
  
  #-----check all variables are numeric or categorical
  if(length(which(sapply(data[,varlist], 
                         function(x) !(is.factor(x)|is.numeric(x)))))>0){
    stop("RPMS works only on numeric or factor data types.")
  }
  
  #-----------------------------------------------------------------------
  
  #------ recurisive partitioning variables ------
  vX=all.vars(rp_equ)[-1]
  
  ##########################   handle categorical variables for C++  ########################
  
  # ------------------identify categorical variable ------
  # need to handle length 1 separately
  if(length(vX)==1) cat_vec <- is.factor(data[,vX]) 
  else
    cat_vec <- sapply(data[, vX], FUN=is.factor)
  
  n_cats<-length(which(cat_vec))
  
  #---------- There are categorical variables ------------
  if(n_cats>0){
    
    # ----- function to turn NA into ? category
    fn_q <- function(x){
      #x<-as.integer(x)
      #x[which(is.na(x))]<- (min(x, na.rm=TRUE)-1)
      nas <- which(is.na(x))
      if(length(nas>0)){
        x <- factor(x, levels=c(levels(x), "?"))
        x[nas]<- factor("?")
      }
      return(as.factor(x))
    } # end internal function fn
    # 
    
    # ---------- apply function to turn each NA into ? category
    if(n_cats==1) {
      data[,vX[which(cat_vec)]] <- fn_q(data[,vX[which(cat_vec)]])
    }
    else{
      data[,vX[which(cat_vec)]] <- lapply(data[,vX[which(cat_vec)]], fn_q)
    }
    
    # ----- store original levels for each categorical variable ------------
    cat_table <- list()
    
    # ----- function to turn categories into integers for C++ 
    for(i in 1:n_cats){
      #print(paste("variable ", vX[which(cat_vec)[i]], " has ", length(levels(data[,vX[which(cat_vec)[i]]])), "levels"))
      cat_table[[i]] <- levels(data[,vX[which(cat_vec)[i]]]) # --store original levels in cat_table
      
      #---- replace original levels with integer -------
      levels(data[,vX[which(cat_vec)[i]]]) <- (1:length(levels(data[,vX[which(cat_vec)[i]]])))
      #print(paste("variable ", vX[which(cat_vec)[i]], "now has ", length(levels(data[,vX[which(cat_vec)[i]]])), "levels"))
    } #end for loop
    
  } # done with if categorical variables
  ######################################################################################################
  
  #--------- reduce data to only y with observations and required variables and no more missing values ---
  n_o<-nrow(data)
  if(n_o <= bin_size) stop("Number of observations must be greater than bin_size")
  data <- na.omit(data[,varlist]) #remove any other missing
  nas<-abs(n_o -nrow(data))
  if(nas > 0)
    warning(paste0("Data had ", nas, " incomplete rows, which have been removed."))
  
  if(nrow(data) <= bin_size) stop("Number of observations must be greater than bin_size")
  
  #===================== Design Information =================================
  #capture design information if equations for use in graphing
  design <- c(weights=NA, strata=NA, clusters=NA)
  # des<-list(weights=~1, strata=~1, clusters=~1)
  
  #------ Sample Weights ----------------------------
  if(is.numeric(weights)) { 
    if(length(weights)!=nrow(data)) stop("Number of design weights != rows of data")
    design[1] <- TRUE 
  } else
    if(length(all.vars(weights))==0) {
      weights <- rep(1, nrow(data))
    } else 
      if(all.vars(weights)[1] %in% names(data)){
        #  des$weights=weights
        weights <- as.numeric(data[,all.vars(weights)[1]])
        if(var(weights)>0) {
          design[1] <- all.vars(weights)[1]
        } 
      } else {stop(paste("Problem with weights:",
                         all.vars(weights)[1], "not in data set"))}
  
  
  #------ Strata Labels ----------------------------
  if(is.numeric(strata) | is.factor(strata)){
    strata<-as.integer(strata)
    design[2] <- TRUE
    if(length(strata)!=nrow(data)) 
      stop("Number of strata labels != rows of data")
  } else
    if(length(all.vars(strata))==0) {strata <- rep(1L, nrow(data))} else
      if(all.vars(strata)[1] %in% names(data)) {
        # des$strata=strata
        design[2] <- all.vars(strata)[1]
        strata <- as.integer(data[,all.vars(strata)[1]])} else 
          stop(paste("Problem with strata:",
                     all.vars(strata)[1], "not in data set"))
  
  #------ Cluster Labels ---------------------------- 
  if(is.numeric(clusters) | is.factor(clusters)){
    clusters<-as.integer(clusters)
    design[3] <- TRUE
    if(length(clusters)!=nrow(data)) 
      stop("Number of cluster labels != rows of data")
  } else
    if(length(all.vars(clusters))==0) {clusters <- seq(1L:nrow(data))} else
      if(all.vars(clusters)[1] %in% names(data)) {
        #      des$clusters <- clusters
        design[3] <- all.vars(clusters)[1]
        clusters <- as.integer(data[,all.vars(clusters)[1]])} else
          stop(paste("Problem with clusters:",
                     all.vars(clusters)[1], "not in data set")) 
  
  des_ind<- !is.na(design)
  
  #================= finished design variables ===============================================  
  
  #-----------------------------------------------------------------------
  # get model matrix and variable data in matrix form
  #-----------------------------------------------------------------------
  e_equ<-formula(paste(all.vars(rp_equ)[1], 1, sep="~"))
  mX<-model.matrix(e_equ, data)
  if(det(t(mX)%*%mX)==0) stop("Model matrix is singular")
  
  X<- data.matrix(as.data.frame(data[,vX]))
  y <- data[,all.vars(e_equ)[1]]
  
  #_________________________________ Done Processing Data ___________________________________________________
  #################################################################################################################
 


  #--- linear modle on root node ----
  modfit0 <- survLm_model(y=as.matrix(y), X=mX, weights=weights, strata=strata, clusters=clusters)

  #===================================== internal function get_trees ==========================================
  get_trees<-function(tnum){

    treeobj<-f_rpms(rp_equ=rp_equ, data=data, weights=weights, strata=strata, clusters=clusters,
                  e_equ=e_equ, e_fn=e_fn, l_fn=l_fn, bin_size=bin_size,
                  y=y, mX=mX, X=X, vX=vX, cat_vec=cat_vec, des_ind=des_ind, n_cats=n_cats, 
                  cat_table=cat_table, modfit0=modfit0)
    treeobj$num<-tnum
    
    return(treeobj)
    
    }
  #===================================== end get_trees ==========================================
  
  if(cores > 1){
    avcores<-parallel::detectCores()
    
    if(cores > avcores){
      cores <- (avcores-1) 
      warning(paste0("Cores set to ", cores))}
    
    tree<-parallel::mclapply(1:f_size, get_trees, mc.cores = cores)
  } 
  else{
    #--- empty vector of trees ----
    # tree <- vector(mode="list", length=f_size) # forest
    # for(i in 1:f_size) tree[[i]] <- get_trees(i)
    tree<-lapply(1:f_size, get_trees)
    }


  callvals <-list(design=design, data.size=length(y), bin_size=bin_size)
  

  f1<-list(callvals=callvals, rp_equ=rp_equ, tree=tree)
 # f1<-list(callvals=callvals, rp_equ=rp_equ, tree=tree, eba=eba)
  class(f1)<-c("rpms_forest")
  
  return(f1)
  
} #end forest

##############################################  END FOREST  ##################################################################



###################################### predict.rpms_forest ###########################################################
#' predict.rpms_forest
#' 
#' @param object  Object inheriting from  \code{rpms_forest}
#' @param newdata data frame with variables to use for predicting new values. 
#' @param ...	further arguments passed to or from other methods.
#' 
#' @description  Gets predicted values given new data based on \code{rpms_forest} model.
#' 
#' @return vector of predicticed values for each row of newdata
#'
#' @export
predict.rpms_forest <- function(object, newdata, ...){
  
  #------------------------------------- process newdata --------------------------------------
  vX=all.vars(object$rp_equ)[-1]
  
  newdata<-newdata[,vX, drop=FALSE]
  
  # ------------------identify categorical variable ------
  # need to handle length 1 separately
  if(length(vX)==1) {
    cat_vec <- is.factor(newdata[,vX])
  } 
  else
    cat_vec <- sapply(newdata[,vX], FUN=is.factor)
  
  n_cats<-sum(cat_vec)
  
  #---------- There are categorical variables ------------
  if(n_cats>0){
    
    # ----- function to turn NA into ? category
    fn_q <- function(x){
      #x<-as.integer(x)
      #x[which(is.na(x))]<- (min(x, na.rm=TRUE)-1)
      nas <- is.na(x)
      if(sum(nas)>0){
        x <- factor(x, levels=c(levels(x), "?"))
        x[nas]<- factor("?")
      }
      return(as.factor(x))
    } # end internal function fn
    
    # ---------- now turn each NA into ? category
    if(sum(cat_vec)==1) {
      newdata[,vX[cat_vec]] <- fn_q(newdata[,vX[cat_vec]])
    }
    else{
      newdata[,vX[cat_vec]] <- lapply(newdata[,vX[cat_vec]], fn_q)
    }
  } #end if n_cats>0
  
  
  #------------------------done processing newdata ---------------------------------------------

  n<-nrow(newdata)
  M<-length(object$tree)
  
  params <-list(...)
  
  if("mode" %in% names(params)) mode <- params$mode else mode="wav"

  # Best estimate predictor worked terrible 

  
  # Weighted average estimator with bias correction
  if(mode=="wav") {  
    
  #---------------------------------------------------------------------------------------------
  # for each tree, find which box the observation falls in for each row of new data
  # boxes is an (n x fsize) matrix where each row contains box index for each tree for observation i
  #
  # estM is an (n x fsize) matrix where each row i contains estimated mean by each tree for observation i

  #index of box and estimated bias correction for each observation
  box_i <- eba <- vector("numeric", n)

    estM <- wghtM <- matrix(vector("numeric", n*M), nrow=n, ncol=M) #makes code faster
    for(j in seq(M)){
      
     # box_i <- apply(covariates(object$tree[[j]]$partition$splits, newdata), 1, function(x) match(1, x))
       box_i <- as.integer(covariates(object$tree[[j]]$partition$splits, newdata) %*% seq(length(object$tree[[j]]$partition$splits)))
       box_i[box_i==0] <- NA
      
      #get the matrix of estimated values for tree i -- nx fsize matrix
      estM[,j] <- as.numeric(object$tree[[j]]$partition$mean)[box_i]
      
      #get the matrix of weights for each value for tree i -- nx fsize matrix
      wghtM[,j] <- object$tree[[j]]$partition$lambda[box_i]
      
    } #end for loop
    
    
    # set sum of weights to 1 
    wghtM[is.na(wghtM)] <- 0
    wghtM <- t(apply(wghtM, 1, function(x) x/sum(x)))
    
    # get bias correction for each value = M*cov
    for(i in seq(n))
      eba[i] <- sum(!is.na(estM[i,]))*cov(wghtM[i,], estM[i,], use="complete.obs")
    
    return(rowSums(estM*wghtM, na.rm=TRUE)-eba)
    
  } # end if mode is wav 
  
  
} #end predict

########################### end predict.rpms_forest ############################




###################################### print.rpms_forest ###########################################################
#' print.rpms_forest
#' 
#' @param x  Object inheriting from  \code{rpms_forest}
#' @param ...	further arguments passed to or from other methods.
#' 
#' @description  Prints information for a given \code{rpms_forest} model.
#' 
#' @return vector of predicticed values for each row of newdata
#'
#' @export
print.rpms_forest <- function(x, ...){
  
  cat("\n")
  cat("    rpms Forest Model     \n")
  cat("##################################################")
  cat("\n")
  cat("---  -----\n \n")
  paste0("Model with ", length(x$tree), " trees")  
  cat("\n")
} 

########################### end print.rpms_forest ############################



################################ boost ##########################################
#
#                   estimate a boosted rpms models             
##############################################################################
#
#' rpms_boost  
#' 
#' @description function for producing boosted rpms models (trees or random forests)
#'
#' @param rp_equ formula containing all variables for partitioning
#' @param data data.frame that includes variables used in rp_equ, e_equ, 
#'        and design information
#' @param weights formula or vector of sample weights for each observation 
#' @param strata formula or vector of strata labels
#' @param clusters formula or vector of cluster labels
#' @param e_equ formula for modeling data in each node
#' @param bin_size numeric minimum number of observations in each node
#' @param gridpts integer number of middle points to do in search
#' @param perm_reps integer specifying the number of thousands of permuation 
#'        replications to use to estimate p-value
#' @param pval numeric p-value used to reject null hypothesis in permutation 
#'        test
#' @param f_size integer specifying the number of trees in the forest (only used 
#'        if model_type is "forest") 
#' @param model_type string: one of "tree" or "forest"
#' 
#' @param times integer specifying number of boosting levels to try.
#' 
#' @return object of class "rpms_boost"
#' 
#' 
#' @examples
#' {
#' # model mean of retirement contributions with a binary tree while accounting 
#' # for clusterd data and sample weights.
#' 
#' rpms_boost(IRAX~EDUCA+AGE+BLS_URBN, data = CE,  weights=~FINLWT21, clusters=~CID, pval=.01)
#'
#' 
#' }
#' @export
rpms_boost <- function(rp_equ, data, weights=~1, strata=~1, clusters=~1, 
                         e_equ=~1, bin_size = NULL, gridpts=3, perm_reps=100L, pval=.05, 
                         f_size=200L, model_type="tree", times=2L){
  
  times<-round(times, 0)
  if(times <= 0) return(NULL)
  
  if("rpms_resids" %in% names(data)) 
    stop("Variable named rpms_resids can't be in dataset please rename")

  y <- all.vars(rp_equ)[1] 
  
  level <- vector(mode="list", length=times) # levels of boosted model
  if(model_type=="tree")
    level[[1]]<-rpms(rp_equ=rp_equ, data=data, weights=weights, strata=strata, clusters=clusters, 
                     e_equ=e_equ, bin_size = bin_size, gridpts=gridpts, perm_reps=perm_reps, pval=pval)
  else
    level[[1]]<-rpms_forest(rp_equ=rp_equ, data=data, weights=weights, strata=strata, clusters=clusters, 
                            bin_size = bin_size, f_size=f_size) 
  
if(times >1){
  
  bvars <- paste(all.vars(rp_equ)[-1], collapse="+")
  b_equ <- as.formula(paste0("rpms_resids~", "~", bvars))
 
  data$rpms_resids <- data[,y]-predict.rpms(level[[1]], data)
  
  if(model_type=="tree"){
    
    for(i in c(2:times)){
      
      est=0
      for(j in 1:(i-1)){est<-est+predict.rpms(level[[j]], data)}
      data$rpms_resids <- data[,y]- est
      
      level[[i]]<-rpms(rp_equ=b_equ, data=data, weights=weights, strata=strata, clusters=clusters, 
                       e_equ=e_equ, bin_size = bin_size, gridpts=gridpts, perm_reps=perm_reps, pval=pval)
    
    }
  } #end tree method
  
  else{
    for(i in c(2:times)){
    
    est=rep(0, nrow(data))
    for(j in 1:i-1){est<-est+predict.rpms_forest(level[[i]], data)}

    data$rpms_resids <- data[,y]- est
    
    level[[i]]<-rpms_forest(rp_equ=b_equ, data=data, weights=weights, strata=strata, clusters=clusters,
                            bin_size=bin_size, f_size=f_size)
    
  } #end for loop

  } #end forest method
  
} #end if times >1 

t1=list(level=level)  
class(t1)<-"rpms_boost"
  
return(t1)

} 
################################## end boost ################################################

##################################### predict.boost_rpms ################################################################
#' predict.rpms_boost
#' 
#' @param object  Object inheriting from  \code{rpms_boost}
#' @param newdata data frame with variables to use for predicting new values. 
#' @param ...	further arguments passed to or from other methods.
#' 
#' @description  Predicted values based on \code{rpms_boost} object
#' 
#' @return vector of predicticed values for each row of newdata
#' 
#'
#' @export
predict.rpms_boost<-function(object, newdata, ...){
  
  levels<-length(object$level)
  
  value<-rep(0, nrow(newdata))

  for(i in 1:levels){
    if(class(object$level[[i]])=="rpms")
    value<-value+predict.rpms(object$level[[i]], newdata)
    
    if(class(object$level[[i]])=="rpms_forest")
      value<-value+predict.rpms_forest(object$level[[i]], newdata)
  } 
  
  return(value)
}



######################## rpms_zinf ############################
#         rpms model when data is zero inflated             
##############################################################
#' rpms_zinf 
#' 
#' @description main function producing a regression tree using variables
#'  from rp_equ to partition the data and fit the model e_equ on each node.
#'  Currently only uses data with complete cases.
#'
#' @param rp_equ formula containing all variables for partitioning
#' @param data data.frame that includes variables used in rp_equ, e_equ, 
#'        and design information
#' @param weights formula or vector of sample weights for each observation 
#' @param strata formula or vector of strata labels
#' @param clusters formula or vector of cluster labels
#' @param e_equ formula for modeling data in each node
#' @param e_fn string name of function to use for modeling 
#'        (only "survLm" is operational)
#' @param l_fn loss function (does nothing yet)
#' @param bin_size numeric minimum number of observations in each node
#' @param gridpts integer number of middle points to do in search
#' @param perm_reps integer specifying the number of thousands of permuation 
#'        replications to use to estimate p-value
#' @param pval numeric p-value used to reject null hypothesis in permutation 
#'        test 
#' 
#' @return object of class "rpms"
#' 
#' @export
rpms_zinf<-function(rp_equ, data, weights=~1, strata=~1, clusters=~1, 
               e_equ=~1, e_fn="survLm", l_fn=NULL, 
               bin_size=NULL, gridpts=3, perm_reps=1000L, pval=.05){
  
  yvar<-all.vars(rp_equ)[1]

  
  #-----get propensity model: p_mod of P(y=1| x) ------------------
  zeros<-which(data[,yvar]==0)
  
  if(length(zeros)<1) stop("There are no zero values in the data")

  newdat<-data
  newdat$I_one <- 0
  newdat$I_one[-zeros] <- 1

  prE <- formula(paste("I_one",paste(all.vars(rp_equ)[-1], collapse="+"), sep="~"))

  if(length(all.vars(e_equ))==0)
    peE <- formula(paste("I_one", 1, sep="~")) else
      peE <- formula(paste("I_one",paste(all.vars(e_equ)[-1], collapse="+"), sep="~"))

  
  p_mod <- rpms(rp_equ=prE, data=newdat, weights=weights, strata=strata, clusters=clusters, 
               e_equ=peE, bin_size=bin_size, perm_reps=perm_reps, pval=pval)
  
  #---------------------- got p_mod --------------------------------------
  
   #-----get model: y_mod of E(y | x, y>0) -------------------------
  
  y_mod <- rpms(rp_equ, data=data[-zeros,], weights=weights, strata=strata, clusters=clusters, 
               e_equ=e_equ, e_fn, l_fn, bin_size, gridpts=gridpts, perm_reps=perm_reps, pval=pval)
  
  t1 <- list(pmod=p_mod, y_mod=y_mod)
  class(t1) <- "rpms_zinf"
  
  return(t1)
  
}

################################## print.rpms_zinf ###################################################################
#' print.rpms_zinf
#' 
#' @param x \code{rpms_zinf} object
#' @param ...	further arguments passed to or from other methods.
#' 
#' @description  print method for class \code{rpms_zinf}
#' 
#'
#' @export
print.rpms_zinf<-function(x, ...){
  
  cat("\n")
  cat("    Zero Inflated RPMS Model     \n")
  cat("##################################################")
  cat("\n\n")
  print(x[[1]], showEnv=FALSE)  
  cat("\n")
  cat("----- RPMS Model for E[ y | x, y >0] --------\n \n")
  
  print(x[[2]], showEnv = FALSE)
  cat("\n")
  
}

###################################### predict.rpms_zinf ################################################################
#' predict.rpms_zinf
#' 
#' @param object  Object inheriting from  \code{rpms_zinf}
#' @param newdata data frame with variables to use for predicting new values. 
#' @param ...	further arguments passed to or from other methods.
#' 
#' @description  Predicted values based on \code{rpms_zinf} model
#' 
#' @return vector of predicticed values for each row of newdata
#' 
#'
#' @export
predict.rpms_zinf<-function(object, newdata, ...){
  if(class(object)!="rpms_zinf") stop("object needs to be of type rpms_zinf")
  if(class(object[[1]])=="rpms")
  return(predict.rpms(object[[1]], newdata)*predict(object[[2]], newdata))
}


#############################################################################
#
#                     Other fun functions
#
###############################################################################

################### r2 function ########################
#' r2
#'
#' @param t1  Object inheriting from  \code{rpms} \code{rpms_forest}
#'            \code{rpms_boost} or \code{rpms_zinf} 
#' @param data data frame with variables used to estimate model
#' @param adjusted TRUE/FALSE whether to compute adjusted R^2
#' 
#' @description  Returns the estimated R^2 statistic for determining
#'               the fit of the given model to the data
#'
#' @return R^2 statistic computed using the model and provided data
#'
#' @export
r2stat<-function(t1, data, adjusted=TRUE){
 
  if(class(t1) %in% c("rpms", "rpms_proj")){
    y <- data[,all.vars(t1$rp_equ)[1]]
    p <- length(all.vars(t1$rp_equ))-1
  }
  else
    if(class(t1) == "rpms_forest"){
      y <- data[,all.vars(t1$tree[[1]]$rp_equ)[1]]
      p <- length(all.vars(t1$tree[[1]]$rp_equ))-1
    }
    else
      if(class(t1) == "rpms_zinf"){
        y <- data[,all.vars(t1[[2]]$rp_equ)[1]]
        p <- length(all.vars(t1[[2]]$rp_equ))-1
      }
      else
      if(class(t1) == "rpms_boost"){
        if(class(t1$level[[1]]) == "rpms"){
          y <- data[,all.vars(t1$level[[1]]$rp_equ)[1]]
          p <- length(all.vars(t1$level[[1]]$rp_equ))-1
        }
        else
          if(class(t1$level[[1]]) == "rpms_forest"){
            y <- data[,all.vars(t1$level[[1]]$tree[[1]]$rp_equ)[1]]
            p <- length(all.vars(t1$level[[1]]$tree[[1]]$rp_equ))-1
          }
        }
      
      else stop(paste0("r2 does not support objects of class ", class(t1)))
      
    residmse <- mean((predict(t1, data) - y)^2, na.rm = TRUE)
    n <- nrow(data)
    
    #definition from Wikipedia
    r2est<-1-(residmse/var(y, na.rm = TRUE))
    
    if(adjusted==TRUE){
      return(1-(1-(r2est*((n-1)/(n-p-1)))))
    }
    else{
      return(r2est)
    }
 
   #kelly's version: return(1-(residmse/var(y, na.rm = TRUE))*(n/(n-1)))
  
}
################### end r2 #############################


##############  rpms_proj ############################################
#' rpms_proj
#'
#' @param object  Object inheriting from  \code{rpms} 
#' @param newdata data frame with variables used to estimate model
#' @param weights vector containing survey weights used in model fit
#' @param weights formula or vector of sample weights for each observation 
#' @param strata formula or vector of strata labels
#' @param clusters formula or vector of cluster labels
#' 
#' @description  Returns a survLm_fit object with coeficients projecting 
#'               new data onto splits from the given rpms model. 
#'
#' @return survLm_fit object
#'
#' @export
rpms_proj <-function(object, newdata, weights=~1, strata=~1, clusters=~1){
  
  #------ Sample Weights ----------------------------
  if(is.numeric(weights)) { 
    if(length(weights)!=nrow(newdata)) stop("Number of design weights != rows of data")
  } else
    if(length(all.vars(weights))==0) {
      weights <- rep(1L, nrow(newdata))
    } else 
      if(all.vars(weights)[1] %in% names(newdata))
      {
        weights <- as.numeric(newdata[,all.vars(weights)[1]])
        if(var(weights)>0) {
        } 
      } else {stop(paste("Problem with weights:",
                         all.vars(weights)[1], "not in data set"))}
  
  #------ Strata Labels ----------------------------
  if(is.numeric(strata) | is.factor(strata)){
    strata<-as.integer(strata)
    if(length(strata)!=nrow(newdata)) 
      stop("Number of strata labels != rows of data")
  } else
    if(length(all.vars(strata))==0) {strata <- rep(1L, nrow(newdata))} else
      if(all.vars(strata)[1] %in% names(newdata)) {
        strata <- as.integer(newdata[,all.vars(strata)[1]])} else 
          stop(paste("Problem with strata:",
                     all.vars(strata)[1], "not in data set"))
  
  #------ Cluster Labels ---------------------------- 
  if(is.numeric(clusters) | is.factor(clusters)){
    clusters<-as.integer(clusters)
    if(length(clusters)!=nrow(newdata)) 
      stop("Number of cluster labels != rows of data")
  } else
    if(length(all.vars(clusters))==0) {clusters <- seq(1L:nrow(newdata))} else
      if(all.vars(clusters)[1] %in% names(newdata)) {
        clusters <- as.integer(newdata[,all.vars(clusters)[1]])} else
          stop(paste("Problem with clusters:",
                     all.vars(clusters)[1], "not in data set")) 
  #===================================================================================  

  
  proj<-linearize(object, newdata, weights = weights, strata = strata, clusters = clusters)
  # LX <- covariates(ln_mod$splits, newdata)
  # 
  # y <- newdata[,all.vars(object$rp_equ)[1]]
  # 
  # t1<-survLm_model(y=y, X=LX, weights=weights, strata=strata, clusters=clusters)
  # 
  #   #survlm_model has "coefficients" "residuals"  
  # 
  # # sptab<-as.data.frame(cbind(object$ln_split, as.numeric(t1$coefficients)))
  # # colnames(sptab)<-c("Splits", "Coefficients")
 
   proj$rp_equ<-object$rp_equ

  class(proj)<-"rpms_proj"
  
  return(proj)
    
}

###################################### predict.rpms_project ################################################################
#' predict.rpms_proj
#' 
#' @param object  Object inheriting from  \code{rpms_zinf}
#' @param newdata data frame with variables to use for predicting new values. 
#' @param ...	further arguments passed to or from other methods.
#' 
#' @description  Predicted values based on \code{rpms_zinf} model
#' 
#' @return vector of predicticed values for each row of newdata
#' 
#'
#' @export
predict.rpms_proj<-function(object, newdata, ...){
  
  LX <- covariates(object$splits, newdata)
  return(LX%*%object$ln_coef)  
}
