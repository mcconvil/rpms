#############################################################################################
# File Containing All support functions for RPM that are intended to be made public
#
# Version 0.6 
#
# Contains: survLm, end_nodes, in_node, node_plot, predict, print, (plot)
#
#
###########################################################################################

#-------------------- Generics ----------------------------



#-----------------------------------------------------------------

############### Length of tree ################################
#' length.rpms
#'
#' @param x rpms tree model object
#'
#' @return number of end-nodes of the tree
#'
#' @aliases length
#'
#' @export
length.rpms <-function(x){
  
  length(x$partition$splits)
}


#################  Front End to survLm_model  ###################
#' survLm
#' 
#' @param e_equ formula representing the equation to fit
#' @param data data.frame 
#' @param weights formula or vector of sample weights for each observation 
#' @param strata formula or vector of strata labels
#' @param clusters formula or vector of cluster labels
#' 
#' @return list containing coefficients, covariance matrix and the residuals
#' 
#' @description  wrapper function for the C++ function survLm_model 
#' 
#' @aliases survLm
#'
#' @keywords internal
survLm<-function(e_equ, data, weights=rep(1, nrow(data)), strata=rep(1L, nrow(data)), clusters=(1L:nrow(data))){
  
  if(length(all.vars(e_equ))> nrow(data)) stop("Number of parameters p > n")
  
  y<-data[,all.vars(e_equ)[1]]
  mX<-model.matrix(e_equ, data)
  
  survLm_model(y=as.matrix(y), X=mX, weights=weights, strata=strata, clusters=clusters)
  
}
####################################################################################################################

# #################################  covariates ##################################
# #---gets regressors for simple function
#
# covariates <- function(splits, data){
#   x<-as.matrix(rep(0,length=nrow(data))) #first column = x satifies no splits
#
#
#   if(length(splits)<=0) return(x) else
#     for(i in 1:length(splits))
#       x<-cbind(x, #if x satisfies split 1 in that column else 0 in that colum
#                eval(parse(text=paste("ifelse(", splits[i], ", 1, 0)")), data),
#                deparse.level=0)
#
#
# #  x <- sapply(splits, function(sp) eval(parse(text=paste("ifelse(", sp, ", 1, 0)")), data))
#
#     #colnames(x)<-paste("E", seq(0:(ncol(x)-1)), sep="")
#     return(x[,-1, drop=FALSE])
# }
# ##################################### End covariates ##########################

################################  covariates ##################################
##---gets regressors for simple function

covariates <- function(splits, data){

  x<-matrix(0, nrow = nrow(data), ncol=length(splits))
  
  #case of no splits
  if(length(splits)==1){
    x[,1]<-rep(1,nrow(data))
    return(x)
  } 

  #if x satisfies split i in that column else 0 in that colum
  for(i in seq(splits)){
    x[eval(parse(text=paste(splits[i])), data),i] <- 1
  }
   
    return(x)
}

##################################### End covariates ##########################



#---gets regressors for simple function
#################################  box_ind #########################################################################################
#' box_ind
#' 
#' @param x \code{rpms} object
#' @param newdata dataframe containing the variables used for the recursive 
#'       partitioning. 
#' 
#' @description  For each row of data, returns a vector indicators whether
#'               observation is in that box or not
#' 
#' 
#' @return Matrix where each row is a vector of indicators whether
#'         observation is in box or not. 
#' 
#' 
#' 
#' @export 
box_ind <- function(x, newdata){
  
  vX=all.vars(x$rp_equ)[-1]
  newdata<-newdata[,vX, drop=FALSE]
  
  
  # ------------------identify categorical variable ------
  # need to handle length 1 separately
  if(length(vX)==1) {
    cat_vec <- is.factor(newdata[,vX, drop=FALSE])
  } 
  else
    cat_vec <- sapply(newdata[,vX, drop=FALSE], FUN=is.factor)
  
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
    
    # ---------- now turn each NA into ? category
    if(length(which(cat_vec))==1) {
      newdata[,vX[which(cat_vec)]] <- fn_q(newdata[,vX[which(cat_vec)]])
    }
    else{
      newdata[,vX[which(cat_vec)]] <- lapply(newdata[,vX[which(cat_vec)]], fn_q)
    }
  } #end if n_cats>0
  
  boxes<-covariates(x$partition[,"splits"], newdata)
  # names(boxes)<-
  return(boxes)
  
}



###################################### prune,rpms ################################################################
#' prune_rpms
#' 
#' @param x \code{rpms} object
#' @param node number of node to prune to.
#' 
#' @description  prune rpms tree to given node  
#' 
#' @return subtree ending clipping off any splits after given node.
#' 
#' @export
prune_rpms<-function(x, node){

  #===========get new frame =============== 

    t_nodes<-x$frame$node
   
    #internal function get_indx finds rows in frame for all children of given node
    get_indx<-function(node){
    c(NULL,
      if(2*node %in% t_nodes) c(which(t_nodes==2*node), get_indx(2*node)),
      if((2*node+1) %in% t_nodes) c(which(t_nodes==(2*node+1)), get_indx(2*node+1)))
  } #end get indx
  
  indx<-get_indx(node)
  
  if(length(indx)==0) return(x)
  
  x$frame<-mark_ends(x$frame[-indx,])

  # done getting new frame  

  x$callvals$pruned <- c(x$callvals$pruned, paste0("at node ", node))
  x$partition<-get_partition(x$frame)
  
  return(x)
   
}
# ################################### End prune_rpms #################################################################


######################################## grow.rpms ################################################
#' grow_rpms
#' 
#' @param x \code{rpms} object
#' @param node node from which to grow tree further 
#' @param data data.frame that includes variables used in rp_equ, e_equ, 
#'        and design information
#' @param weights formula or vector of sample weights for each observation 
#' @param strata formula or vector of strata labels
#' @param clusters formula or vector of cluster labels
#' @param pval numeric p-value used to reject null hypothesis in permutation test 
#' @param bin_size numeric minimum number of observations in each node
#' 
#' @description  grow an rpms tree from a given node 
#' 
#' @return rpms tree expanded from node.
#' 
#' @export
grow_rpms <- function(x, node, data, weights=~1, strata=~1, clusters=~1, pval=NA, bin_size=NA){
  
  
  #============= format data set ================================================
  
  des_ind<- !is.na(x$callvals$design)
  if(sum(des_ind)==0) 
    varlist <- unique(c(all.vars(x$rp_equ), all.vars(x$e_equ)))
  else 
    varlist <- unique(c(all.vars(x$rp_equ), all.vars(x$e_equ), x$callvals$design[which(des_ind)]))
  
  #---------- check all variables are in data set --------------
  if(!all(varlist %in% names(data))){
    e1ind <- which(!(varlist %in% names(data)))[1]
    stop(paste("Variable", varlist[e1ind], "not in dataset"))
  }
  
  #find rows of data in the starting node
  if(!(node %in% x$partition$nodes)) x<-prune_rpms(x, node)
  s=which(end_nodes(x, data)==node)
  data<-data[s,]
  
  #================= finished design variables ===============================================  
  
  #-----------------------------------------------------------------------
  # get model matrix and variable data in matrix form
  #-----------------------------------------------------------------------
  mX<-model.matrix(x$e_equ, data)
  if(det(t(mX)%*%mX)==0) stop("Model matrix is singular")
  
  #-----check all variables are numeric or categorical
  if(length(varlist)<2) stop("rpms model needs at least two variables")
  if(length(which(sapply(data[,varlist], 
                         function(x) !(is.factor(x)|is.numeric(x)))))>0){
    stop("RPMS works only on numeric or factor data types.")
  }
  
  #-----------------------------------------------------------------------
  
  #------ recurisive partitioning variables ------
  vX=all.vars(x$rp_equ)[-1]
  
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
    
    # ---------- turn each NA into ? category
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
  data <- na.omit(data[,varlist]) #remove any other missing
  nas<-abs(n_o -nrow(data))
  if(nas > 0)
    warning(paste0("Data had ", nas, " incomplete rows, which have been removed."))
  
  
  #===================== Design Information =================================
  #capture design information if equations for use in graphing
  design <- c(weights=NA, stratum=NA, clusters=NA)
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
  
  
  #------ recurisive partitioning variables ------
  X<- data.matrix(as.data.frame(data[,vX]))
  y <- data[,all.vars(x$e_equ)[1]]
  
  if(is.na(pval)) pval<-x$callvals$pval
  if(is.na(bin_size)) bin_size <- x$callvals$bin_size
  
  
  frame2 <-
    rbind_splits(null_split(), 
                 split_rpms(node=node, y=y, mX=mX, X=X, vnames=vX, cat_vec=cat_vec, 
                            weights=weights, strata=strata, clusters=clusters, des_ind=des_ind, 
                            bin_size=bin_size, gridpts=x$callvals$gridpts, perm_reps=x$callvals$perm_reps, pval=pval))
  
  #if there were splits make new frame and partition otherwise just return original tree
  if(length(frame2$node)>1){
 
    #======================================= return to original the factor levels in the data ==============================
    if(n_cats>0){
      for(i in 1:n_cats){
        #---- return to original levels -------
        levels(data[,vX[which(cat_vec)[i]]]) <- unlist(cat_table[[i]])
      }
    }
    
    # put original level names in the frame
    cat_splits<-which(frame2$cat==1)
    
    for(i in cat_splits){
      vis <- which(vX[cat_vec]==frame2$var[[i]]) #variable location in cat_table
      frame2$xval[i] <- list(cat_table[[vis]][unlist(frame2$xval[i])]) #replace numbers with factor names
    }
    #======================================================================================================================
    
    frame2 <-  make_nice(frame2) # format resulting tree frame for R 
    # puts an 'E' after each end node on the frame
    
    # get first half of frame without ends marked
    frame1 <- x$frame[,-which(names(x$frame)=="end")]
    #new frame is combinded old and new
    frame<-rbind(frame1, frame2)
    frame <- mark_ends(frame)
    
    x$frame <- frame
    
    partition<-get_partition(frame)
    row.names(partition)<-NULL
    x$partition <- partition
    
  }# else return the original tree

  return(x)
  
}
######################################## end grow_rpms #################################################





#---gets regressors for simple function
#################################  boxes #########################################################################################
#' boxes
#' 
#' @param x \code{rpms} object
#' 
#' @description  returns end boxes that partition the data
#' 
#' 
#' @return data.frame including end_node, sample size, splits, and values for each end node
#' 
#' 
#' 
#' @export 
boxes <- function(x){
  
  return(x$partition)
  
}

##################################### End boxes ##########################


#---gets regressors for simple function
#################################  linearize #########################################################################################
#' linearize
#' 
#' @param x \code{rpms} object
#' @param data data.frame 
#' @param weights formula or vector of sample weights for each observation 
#' @param strata formula or vector of strata labels
#' @param clusters formula or vector of cluster labels
#' @param type is on of "part" or "lin"
#' 
#' @description  returns a linerized version of the splits. The coefficients represent 
#'               the effect that each split has on the mean
#' 
#' @return data.frame including splits and estimates for the coefficient and their standard errors
#' 
#' @export 
########################################## Linearize ###########################################
#takes frame (f1) from tree and outputs linear logical splits
linearize<-function(x, data, weights=~1, strata=~1, clusters=~1, type="part"){

  if(class(x)!="rpms") return("rpms::linearize is for rpms objects only")
  
  #============= format data set ================================================
  varlist <- unique(c(all.vars(x$rp_equ), all.vars(weights), 
                      all.vars(strata), all.vars(clusters)))  
  
  
  #---------- check all variables are in data set --------------
  if(!all(varlist %in% names(data))){
    e1ind <- which(!(varlist %in% names(data)))[1]
    stop(paste("Variable", varlist[e1ind], "not in dataset"))
  }
  
  #-----check all variables are numeric or categorical
  if(length(which(sapply(data[,varlist, drop=FALSE], 
                         function(x) !(is.factor(x)|is.numeric(x)))))>0){
    stop("RPMS works only on numeric or factor data types.")
  }
  
  #------ recurisive partitioning variables ------
  vX=all.vars(x$rp_equ)[-1]
  
  
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
   #------------------------------------------------------ 
    
    # ---------- apply function to turn each NA into ? category ---------
    if(n_cats==1) {
      data[,vX[which(cat_vec)]] <- fn_q(data[,vX[which(cat_vec)]])
    }
    else{
      data[,vX[which(cat_vec)]] <- lapply(data[,vX[which(cat_vec)]], fn_q)
    }
    
    
  } # done with if categorical variables
  ######################################################################################################
  
  #===================== Design Information =================================
  
  #------ Sample Weights ----------------------------
  if(is.numeric(weights)) { 
    if(length(weights)!=nrow(data)) stop("Number of design weights != rows of data")
  } else
    if(length(all.vars(weights))==0) {
      weights <- rep(1, nrow(data))
    } else 
      if(all.vars(weights)[1] %in% names(data))
      {
        #  des$weights=weights
        weights <- as.numeric(data[,all.vars(weights)[1]])
      } else {stop(paste("Problem with weights:",
                         all.vars(weights)[1], "not in data set"))}
  
  
  #------ Strata Labels ----------------------------
  if(is.numeric(strata) | is.factor(strata)){
    strata<-as.integer(strata)
    if(length(strata)!=nrow(data)) 
      stop("Number of strata labels != rows of data")
  } else
    if(length(all.vars(strata))==0) {strata <- rep(1L, nrow(data))} else
      if(all.vars(strata)[1] %in% names(data)) {
        strata <- as.integer(data[,all.vars(strata)[1]])} else 
          stop(paste("Problem with strata:",
                     all.vars(strata)[1], "not in data set"))
  
  #------ Cluster Labels ---------------------------- 
  if(is.numeric(clusters) | is.factor(clusters)){
    clusters<-as.integer(clusters)
    if(length(clusters)!=nrow(data)) 
      stop("Number of cluster labels != rows of data")
  } else
    if(length(all.vars(clusters))==0) {clusters <- seq(1L:nrow(data))} else
      if(all.vars(clusters)[1] %in% names(data)) {
        clusters <- as.integer(data[,all.vars(clusters)[1]])} else
          stop(paste("Problem with clusters:",
                     all.vars(clusters)[1], "not in data set")) 
  
  #===================================================================================  
  
  #--------- reduce data to only those with observations and required variables 
  n_o<-nrow(data)
  data <- na.omit(data[,varlist]) #remove any other missing
  nas<-abs(n_o -nrow(data))
  if(nas > 0)
    warning(paste0("Data had ", nas, " incomplete rows, which have been removed."))

  y <- data[,all.vars(x$rp_equ)[1]]

if(type=="lin")
  lt<-rbind("1", lin_splits(x$frame))
else
  lt<-x$partition$splits

   # row.names(lt)<-NULL
  
  
  fit<-survLm_model(y, covariates(lt, data), weights, strata, clusters)
  
  
  return(list(splits=lt, ln_coef=as.numeric(fit$coefficients),ln_coef_cov=fit$covM))
  
}#end function
##################################### End Linearize #########################################



###################################### end_nodes ################################################################
#' end_nodes
#' 
#' @param object \code{rpms} object
#' @param newdata data.frame
#' 
#' @description  Either a vector of end-node labels for each opbservation in newdata or 
#'               a vector of the endnodes in the tree model if newdata is not provided.
#' 
#' @return vector of end_node labels
#' 
#' 
#' @examples 
#' {
#' # model mean of retirement account value for households with reported 
#' # retirment account values > 0 using a binary tree while accounting for 
#' # clusterd data and sample weights.
#' 
#' s1<- which(CE$IRAX > 0)
#' r1 <-rpms(IRAX~EDUCA+AGE+BLS_URBN, data = CE[s1,],  weights=~FINLWT21, clusters=~CID) 
#'  
#'  end_nodes(r1)
#' }
#' 
#' @export
end_nodes <- function(object, newdata=NULL){
  
 if(is.null(newdata))
   return(object$partition[,"node"]) 
  
  else {
    
    vX=all.vars(object$rp_equ)[-1]
    varlist <- unique(c(all.vars(object$rp_equ[-1]), all.vars(object$e_equ)[-1]))  
    
    newdata<-newdata[,varlist, drop=FALSE]
    
    # ------------------identify categorical variable ------
    # need to handle length 1 separately
    if(length(vX)==1) {
      cat_vec <- is.factor(newdata[,vX])
    } 
    else
      cat_vec <- sapply(newdata[,vX], FUN=is.factor)
    
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
      
      # ---------- now turn each NA into ? category
      if(length(which(cat_vec))==1) {
        newdata[,vX[which(cat_vec)]] <- fn_q(newdata[,vX[which(cat_vec)]])
      }
      else{
        newdata[,vX[which(cat_vec)]] <- lapply(newdata[,vX[which(cat_vec)]], fn_q)
      }
    } #end if n_cats>0
    
    #newdata <- na.omit(newdata) #remove any other missing
    if(length(which(is.na(newdata)))>0) stop("Remove any records with missing numeric predictor variables")
    
    
    return(apply(covariates(object$partition[,"splits"], newdata), 1, function(x) object$partition$node[which(x==1)] ))
 } #end if new data
   
  
}
# ################################### End end_nodes #################################################################


###################################### in_node ################################################################
#' in_node
#' 
#' @param x \code{rpms} object
#' @param node integer label of the desired end-node.
#' @param data dataframe containing the variables used for the recursive 
#'       partitioning. 
#' 
#' @description  Get index of elements in dataframe that are in the specified 
#'               end-node of an \code{rpms} object.  A "which" function for end-nodes.
#' 
#' @return vector of indexes for observations in the end-node. 
#' 
#' 
#' @examples
#' {
#' # model mean of retirement account value for households with reported 
#' # retirment account values > 0 using a binary tree while accounting for 
#' # clusterd data and sample weights.
#' 
#' s1<- which(CE$IRAX > 0)
#' r1 <-rpms(IRAX~EDUCA+AGE+BLS_URBN, data = CE[s1,],  weights=~FINLWT21, clusters=~CID) 
#' 
#' # Get summary statistics of CUTENURE for households in end-nodes 7 and 8 of the tree
#'
#' if(7 %in% end_nodes(r1)) 
#'   summary(CE$CUTENURE[in_node(node=7, r1, data=CE[s1,])])
#' if(8 %in% end_nodes(r1)) 
#'   summary(CE$CUTENURE[in_node(node=8, r1, data=CE[s1,])])
#' }
#' 
#' @export 
in_node<-function(x, node, data){
  
  if(!(node %in% x$partition$nodes)) x<-prune_rpms(x, node)
    
    return(which(end_nodes(x, data)==node))
  
}

###################################### End in_node ############################

################################## node_plot ##################################
#' node_plot
#' 
#' @param object \code{rpms} object
#' @param node integer label of the desired end-node. 
#' @param data data.frame that includes variables used in rp_equ, e_equ, 
#'        and design information
#' @param variable string name of variable in data to use as x-axis in plot
#' @param ...	further arguments passed to plot function.      
#' 
#' 
#' @description  plots end-node of object of class \code{rpms}
#' 
#' @import graphics
#' 
#' @examples{
#' 
#' # model mean of retirement account value for households with reported 
#' # retirment account values > 0 using a binary tree while accounting for 
#' # clusterd data and sample weights.
#' 
#' s1<- which(CE$IRAX > 0)
#' r1 <-rpms(IRAX~EDUCA+AGE+BLS_URBN, data = CE[s1,],  weights=~FINLWT21, clusters=~CID) 
#' 
#' # plot node 6 if it is an end-node of the tree
#' if(6 %in% end_nodes(r1))
#'   node_plot(object=r1, node=6, data=CE[s1,])
#' 
#' # plot node 6 if it is an end-node of the tree
#' if(8 %in% end_nodes(r1))
#'   node_plot(object=r1, node=8, data=CE[s1,])
#' 
#' }
#'
#' @export
#'
node_plot<-function(object, node, data, variable=NA, ...){
  
  pars<-list(...)
  
  if("col" %in% names(pars)) col=pars$col else col="blue"
  if("lwd" %in% names(pars)) lwd=pars$lwd else lwd=3
  if("lty" %in% names(pars)) lty=pars$lty else lty=1
  
  #if no variable provided use first variable in estimating equation
  #as x axis to plot on
  if(is.na(variable)){
    if(length(all.vars(object$e_equ))>1)
      variable <- all.vars(object$e_equ)[2]
    else variable <- all.vars(object$rp_equ)[2]
    
  } #if variable provided check to make sure it is in the dataset
  else  if(!(variable %in% names(data) && is.character(variable)))
    stop(paste("variable", paste(variable, collapse=""), 
               "is not a name in the data set"))

  
  # y variable is the y from the estimating equation 
  yvariable<-all.vars(object$e_equ)[1]
  
  #which element has coefficients for that node
  nind <- which(object$partition[,"node"]==node)
  
  # if the variable chosen is numeric produce the best fit line in the plot
  if(is.numeric(data[, variable])) coline <-TRUE else coline <-FALSE
  
  #get index of observations contained in end node
  eindx<-in_node(node=node, x=object, data=data)
  plot(data[eindx, c(variable, yvariable)])
  title(main=paste0("Node ", node), ylab=yvariable, xlab=variable)
  
  #find variable in estimating equation being used to graph against
  vindx <- match(variable, all.vars(object$e_equ))

  if(coline && !is.na(vindx) && vindx>1)
      abline(coef=unlist(object$partition$value[nind])[c(1,vindx)], col=col, lwd=lwd, lty=lty)
  
  
} #end node_plot
################################### End plot.rpms #################################################################



###################################### predict.rpms ################################################################
#' predict.rpms
#' 
#' @param object  Object inheriting from  \code{rpms}
#' @param newdata data frame with variables to use for predicting new values. 
#' @param ...	further arguments passed to or from other methods.
#' 
#' @description  Predicted values based on \code{rpms} object
#' 
#' @return vector of predicticed values for each row of newdata
#' 
#'
#' @examples{
#' 
#' # get rpms model of mean Soc Security income for families headed by a 
#' # retired person by several factors
#' r1 <-rpms(SOCRRX~EDUCA+AGE+BLS_URBN+REGION, 
#'           data=CE[which(CE$INCNONWK==1),], clusters=~CID) 
#' 
#' r1
#' 
#' # first 10 predicted means
#' predict(r1, CE[10:20, ])
#' 
#' }
#'
#'@export
predict.rpms<-function(object, newdata, ...){
  
  # pars<-list(...)  #currently does not take other parameters

  #----------------------- Prepare data for predict -----------------------------------  
  
  new_equ <- object$e_equ[-2]
  vX=all.vars(object$rp_equ)[-1]
  varlist <- unique(c(all.vars(object$rp_equ[-1]), all.vars(object$e_equ)[-1]))  
  
  newdata<-newdata[,varlist, drop=FALSE]
  
  # ------------------identify categorical variable ------
  # need to handle length 1 separately
  if(length(vX)==1) {
    cat_vec <- is.factor(newdata[,vX])
  } 
  else
    cat_vec <- sapply(newdata[,vX], FUN=is.factor)
  
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
    
    # ---------- now turn each NA into ? category
    if(length(which(cat_vec))==1) {
      newdata[,vX[which(cat_vec)]] <- fn_q(newdata[,vX[which(cat_vec)]])
    }
    else{
      newdata[,vX[which(cat_vec)]] <- lapply(newdata[,vX[which(cat_vec)]], fn_q)
    }
  } #end if n_cats>0
  
  if(length(which(is.na(newdata)))>0) stop("Remove any records with missing numeric predictor variables")
  #-------------------------------------------- Done preparing data --------------------------------------------
  
   
  #coefficient values for each end node
  
  #beta_mat<-do.call(rbind, object$partition$value)
 
    beta_mat<-matrix(unlist(object$partition$value),
                     nrow=length(object$partition$value), byrow = TRUE)

  if(length(all.vars(new_equ))==0){
    #no linear function just return beta
    return(as.vector(apply(covariates(object$partition[,"splits"], newdata), 1, function(x) beta_mat[match(1, x)])))
  }
  else{ # linear model
    
    #get model matrix X
    X <- model.matrix(new_equ, newdata) #x values
    
    #-- beta values for that each observation
    b<-as.matrix(apply(covariates(object$partition[,"splits"], newdata), 1, function(x) beta_mat[match(1, x), ,drop=FALSE]))

    return(colSums(t(X)*b))
    
  }# end if lin model
  
}

################################### End predict #################################################################


################################# print.rpms ###################################################################
#' print.rpms
#' 
#' @param x \code{rpms} object
#' @param ...	further arguments passed to or from other methods.
#' 
#' @description  print method for class \code{rpms}
#' 
#'
#' @export
print.rpms<-function(x, ...){
  
  if(class(x)!="rpms") stop("argument must be of class rpms") else t1 <- x
  
  if(!is.null(x$rp_equ)){
    cat("\n")
    cat("RPMS Recursive Partitioning Equation \n")
    print(t1$rp_equ, showEnv=FALSE)
    cat("\n")
    cat("Estimating Equation \n")
    print(t1$e_equ, showEnv=FALSE)
    cat("\n")
  }

  #callvals not in btree
  if("callvals" %in% names(x)){
    
    if(!is.null(x$callvals)){
      des_ind<- !is.na(x$callvals$design)
      des_string <- if(sum(des_ind)==0) "Simple Random Sample"
      else paste(c("unequal probability of selection", "stratified", "clustered")[which(des_ind)], "sample design", sep=", ")
      
      cat("\n")
      print(des_string)
      
      if(!is.null(t1$callvals$pruned)){
        print(paste0("Pruned rpms tree  ", t1$callvals$pruned), showEnv=FALSE)
        cat("\n")
        
      }# end if pruned tree
    } # end if callvals null
  } # end if object has "callvals"/ btree trees don't
  

    
    #make matrix of estimated coefficients (colums) for each end node (row)
    ends<-which(t1$frame$end=="E")
    
    if("value" %in% names(frame)){
      if(attr(terms(t1$e_equ), "intercept")==1)
        coef_names<-list(node=t1$frame$node[ends], coefficients=c("1", all.vars(t1$e_equ)[-1]))
      else
        coef_names<-list(node=t1$frame$node[ends], coefficients=all.vars(t1$e_equ)[-1])
      coef_mat<-matrix(data=unlist(t1$frame$value[ends]),
                       nrow=length(ends), 
                       ncol=length(coef_names$coefficients), 
                       byrow=TRUE, dimnames=coef_names)
    }else
      coef_mat<-NULL
   
    
  
    if("r2" %in% names(t1)){print(paste0("R-squared of model: ", t1$r_2))}
    
    sptab<-lin_splits(t1$frame)
   # if(dim(sptab)==0) sptab <- as.matrix(sptab)
    colnames(sptab)<-c("Splits")
    cat("\n")
    cat("===================== Tree Model =================== \n \n")
    
    print(sptab, quote=F, zero.print=".")
    cat("\n")
    
    if(!is.null(coef_mat)){
      print(coef_mat, quote=F, zero.print=".")
      cat("\n \n")
    }

}#end print.rpms

################################### End print.rpms #################################################################


################################### plot.rpms #################################################################
#
#  Workhorse function of plot function
#
get_tree_plot <- function(frame, cnode, x, y, xstep, ystep, yname) {
  if(frame[frame$node==cnode, ]$end=="E") {
    if("value" %in% names(frame)){
      text(x, y-.02, paste0(yname, " = ", round(unlist(frame[frame$node==cnode, ]$value), 2)), cex = 0.8)
    } else{text(x, y-.02, yname, cex = 0.8) }
  } else {
    if(frame[frame$node==2*cnode, ]$cat == 0) {
      text(x, y, paste(frame[frame$node==2*cnode, ]$var, "<=", round(unlist(frame[frame$node==2*cnode, ]$xval), 2)), cex = 0.8)
      lines(c(x, x - xstep), c(y-.02, y + ystep +.02))
      lines(c(x, x + xstep), c(y-.02, y + ystep +.02))
      get_tree_plot(frame, cnode=2*cnode, x - xstep, y + ystep, xstep*(2/3), ystep, yname)
      get_tree_plot(frame, cnode=(2*cnode +1), x + xstep, y + ystep, xstep*(2/3), ystep, yname)
    } else {
      categories <- paste(frame[frame$node==2*cnode, ]$xval, collapse = "\n")
      text(x, y, paste(frame[frame$node==2*cnode, ]$var, "=", categories), cex = 0.8)
      lines(c(x, x - xstep), c(y-.02, y + ystep +.02))
      lines(c(x, x + xstep), c(y-.02, y + ystep +.02))
      get_tree_plot(frame, cnode=2*cnode, x - xstep, y + ystep, xstep*(2/3), ystep, yname)
      get_tree_plot(frame, cnode=(2*cnode+1), x + xstep, y + ystep, xstep*(2/3), ystep, yname)
    }
  }
}  # End work horse for plot program

################################# plot.rpms ###################################################################
#' plot.rpms
#'
#' Function to plot rpms tree structure using frame information
#'
#' @param x \code{rpms} object
#' @param ...	further arguments passed to or from other methods.
#'
#' @description  plot method for class \code{rpms}
#'
#' @aliases plot
#' 
#' @export
plot.rpms <- function(x, ...) {
  
  if (!inherits(x, "rpms")) {
    stop("Input must be a rpms model.")
  }
  
  rpms_tree <- x
  rm(x)
  yname <- all.vars(rpms_tree$rp_equ)[1] 
  
  pars<-list(...)
  
  if("x_0" %in% names(pars)){x <- as.numeric(pars$x_0)}else{x <- 0.5} 
  if("y_0" %in% names(pars)){y <- as.numeric(pars$y_0)}else{y <- 0.95} 
  if("xstep" %in% names(pars)){xstep <- as.numeric(pars$xstep)}else{xstep <- 0.2} 
  if("ystep" %in% names(pars)){ystep <- as.numeric(pars$ystep)}else{ystep <- -0.15} 
  
  tree_frame <- rpms_tree$frame
  
  # Create an empty plot
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "")
  
  # Call recursive function to plot tree
  get_tree_plot(tree_frame, cnode=1, x, y, xstep, ystep, yname)
}

# 
# ################################### End plot.rpms #################################################################






########################################### function ################################################################################

##################################################################################################################################

########################################### function ################################################################################

##############################################################################################################################

########################################### function ################################################################################

##############################################################################################################################


