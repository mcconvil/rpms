#Write a Recursive Partition Algorithm with Modeling
#Version 0.6.0
#Daniell Toth 09/04/2015
#Handles one-stage stratified and clusters samples with unequal weights
#Does not have plot function
#Uses C++ functions through Rcpp
#Uses variable select; M=1000 permutations 

######################################################################################
#                                                                                    #
#                                       Programs for RPMS                            #
#                                                                                    #
######################################################################################

#returns the row number of every row with a missing value
na.rows<-function(x){
  if(ncol(x)<=1) return(which(is.na(x))) else 
    return(which(as.vector(apply(x, 1, function(x){any(is.na(x))}))))
}


##############################################################################################
#
#   functions to make output of rpms more readable
#
##############################################################################################

#------------------------ gives linear version of splits for output -----------------------------
# f1 is frame of an rpms object
#recursive function starts with node number 1 
lin_splits<-function(f1, n_num=1, sp=NULL){
  
  if(length(f1$node)==1) return(as.matrix("1"))
  
  if(!(n_num %in% f1$node)) return(NULL) #end
  
  i <-which(f1$node==n_num)
  and<-ifelse(is.null(sp),"", "&") #No & for first split
  
  if(n_num==1) #first split
    return(rbind(lin_splits(f1, 2), lin_splits(f1, 3))) 
  
  if((n_num %% 2)==0) #left side 
    if(f1$cat[i]==0){
      #  print(paste("sp is", sp, "in cat"))
      sp<-paste(sp, and, paste(f1$var[i], "<=", f1$xval[i], sep=" "))
      return(rbind(sp, lin_splits(f1, 2*n_num, sp), lin_splits(f1, 2*n_num+1, sp)))
    }  
  else {
    # print(paste("sp is", sp, "in cat"))
    sp<-paste(sp, and, paste(f1$var[i], " %in% c('", paste(f1$xval[i][[1]], collapse="','"), "')", sep=""))
    return(rbind(sp, lin_splits(f1, 2*n_num, sp), lin_splits(f1, 2*n_num+1, sp)))
  }
  else #right side
    if((2*n_num) %in% f1$node)
      if(f1$cat[i]==0){
        sp<-paste(sp, and,paste(f1$var[i], ">", f1$xval[i], sep=" "))
        return(rbind(lin_splits(f1, 2*n_num, sp), lin_splits(f1, 2*n_num+1, sp)))
      }        
  else {
    sp<-paste(sp, and, paste(f1$var[i], " %in%  c('", paste(f1$xval[i][[1]], collapse="','"), "')", sep=""))
    return(rbind(lin_splits(f1, 2*n_num, sp), lin_splits(f1, 2*n_num+1, sp)))
  }
  else return(NULL)
} #end function get lin_splits
#-------------------------------------------------------------------------------------------------------------


####################### Mark the end nodes as ends in Frame ######################
#returns the frame with "E" at the end of each frame
mark_ends<-function(frame, node=1){
  frame$end <- "" #sets new colum to blank string vector
  
  get_ends<-function(frame, n_num=node){
    if((2*n_num) %in% frame$node) 
      return(c(get_ends(frame, 2*n_num), get_ends(frame, 2*n_num+1))) 
    else return(n_num)
  }
  
  frame$end[which(frame$node %in% get_ends(frame))]="E"
  
  return(frame)
}

############################ get_partition   ##############################################
# f1 is frame of an rpms object
#recursive function starts with node number 1 
get_partition<-function(f1, n_num=1, sp=NULL){
  
  if(length(f1$node)==1) return(cbind(f1[1, c("node", "n")], splits="1", f1[1, c("value", "cvar", "mean")]))
  
  i <-which(f1$node==n_num)
  
  if(n_num==1) {
    if(f1$end[i]=="E") return(NULL) 
    else return(rbind(get_partition(f1, 2*n_num, sp), get_partition(f1, 2*n_num+1, sp)))
  }
  
  and<-ifelse(is.null(sp), "", " & ") #No & for first split
  
  #update the split 
  
  #left side 
  if((n_num %% 2)==0){
    if(f1$cat[i]==0) #numeric
      sp<-paste0(sp, and, paste(f1$var[i], "<=", f1$xval[i]))
    else #cat
      sp<-paste0(sp, and, paste0(f1$var[i], " %in% c('", paste(f1$xval[i][[1]], collapse="','"), "')"))
  } #end left side 
  #right side
  else{
    if(f1$cat[i]==0) #numeric
      sp<-paste0(sp, and, paste(f1$var[i], ">", f1$xval[i]))
    else 
      sp<-paste0(sp, and, paste0(f1$var[i], " %in%  c('", paste(f1$xval[i][[1]], collapse="','"), "')"))
  } #right side
  
  #done updating split
  
  if(f1$end[i]=="E") {
   return(cbind(f1[i, c("node", "n")], splits=sp, f1[i, c("value", "cvar", "mean")]))
  }
  else
    return(rbind(get_partition(f1, 2*n_num, sp), get_partition(f1, 2*n_num+1, sp)))
    
} #end funciton  
############################### end get_partition   ##########################################


#################################################################################################
#
#  make nice puts the List returned from Cpp into a data.frame and puts nice names on the labels
#
#################################################################################################
make_nice<-function(frame){
  
  nice.frame<-as.data.frame(frame[-c(4,7,8)])
  names(nice.frame)<-c("node", "cat", "var", "n", "loss", "mean")
  
  nice.frame$xval=frame$xval
  
  nice.frame$value=frame$value
  names(nice.frame$value)<-"value"
  
  nice.frame$cvar=frame$cvar
  names(nice.frame$cvar)<-"cvar"
  
  return(nice.frame[c("node", "cat", "var", "xval", "n", "loss", "value", "cvar", "mean")])
  
}#end function 
######################  end make nice #############################################################




#####################################################################################################
################################################################################
#                                                                              #
#                               MAIN Function rpm                              #
#                                                                              #
################################################################################
#' rpms  
#' 
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
#' @param l_fn loss function (ignored)
#' @param bin_size integer specifying minimum number of observations in each node
#' @param gridpts integer number of middle points to do in search; set to n for 
#'        categorical variables when e_equ is used.
#' @param perm_reps integer specifying the number of thousands of permutation 
#'        replications to use to estimate p-value
#' @param pval numeric p-value used to reject null hypothesis in permutation 
#'        test 
#' 
#' @return object of class "rpms"
#' 
#' 
#' @description main function producing a regression tree using variables
#'  from rp_equ to partition the data and fit the model e_equ on each node.
#'  Currently only uses data with complete cases of continuous variables.
#' 
#' @examples
#' {
#' # model mean of retirement account value for households with reported 
#' # retirment account values > 0 using a binary tree while accounting for 
#' # clusterd data and sample weights.
#' 
#' s1<- which(CE$IRAX > 0)
#' rpms(IRAX~EDUCA+AGE+BLS_URBN, data=CE[s1,], weights=~FINLWT21, clusters=~CID)
#'
#'                  
#' # model linear fit between retirement account value and amount of income
#' # conditioning on education and accounting for clusterd data for households 
#' # with reported retirment account values > 0
#' 
#' rpms(IRAX~EDUCA, e_equ=IRAX~FINCBTAX, data=CE[s1,], weights=~FINLWT21, clusters=~CID)    
#' 
#' }
#' @export
#' @aliases rpms
#' 
rpms<-function(rp_equ, data, weights=~1, strata=~1, clusters=~1, 
               e_equ=~1, e_fn="survLm", l_fn=NULL, 
               bin_size=NULL, gridpts=3, perm_reps=1000L, pval=.05){
  
#--------------------------------------------------------------------------------------#
#                         check parameters and format variables                        #
#--------------------------------------------------------------------------------------#
  
  # Unnecessary l_fn is ignored ##########
  # if(is.null(l_fn)) l_fn <- function(x){sum(x^2)}
  # else if(!(class(l_fn)=="function")) stop("l_fn is not of type function")
  

  if(is.null(bin_size)) bin_size <- ceiling(nrow(data)^(11/20))
  else 
    if(bin_size<=2) stop("bin_size must be set > 0")  
  
  
  perm_reps=round(perm_reps)
  if(pval <= 0) stop("pval must be set > 0")
  if(pval > 1) stop("pval should not be greater than 1")
  if(perm_reps < (.5/pval)) stop("perm_reps are too low for stated pval")
  
  
  #------- check variables in equations ------
  if(length(all.vars(rp_equ))==0) 
    stop("no variable for recursive partition given in formula") 
  
  if(length(all.vars(e_equ))==0)
    e_equ<-formula(paste(all.vars(rp_equ)[1], 1, sep="~")) else
      if(all.vars(e_equ)[1]!=all.vars(rp_equ)[1]) 
        stop("Dependent variable must be same for rp_equ as e_equ")
  
  
#============= format data set ================================================
  varlist <- unique(c(all.vars(rp_equ), all.vars(e_equ), all.vars(weights), 
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
        design[1] <- all.vars(weights)[1]
        weights <- as.numeric(data[,all.vars(weights)[1]])
        if(var(weights)==0) {
          design[1] <- NA
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
  mX<-model.matrix(e_equ, data)
  if(det(t(mX)%*%mX)==0) stop("Model matrix is singular")
  
  X<- data.matrix(as.data.frame(data[,vX]))
  y <- data[,all.vars(e_equ)[1]]
  
#_________________________________ Done Processing Data ___________________________________________________


  ##################################  Recursive Partitioning  ####################################################
  #          calls C++ funtions get_node and split_rpms 
  ################################################################################################################
  frame <-
    rbind_splits(get_node(node=1, cat=as.integer(NA), vname="Root", y=as.matrix(y), weights=weights, mxval=as.matrix(NA), s=as.matrix(0), 
                          modfit=survLm_model(y=as.matrix(y), X=mX, weights=weights, strata=strata, clusters=clusters)), 
                 split_rpms(node=1, y=y, mX=mX, X=X, vnames=vX, cat_vec=cat_vec, 
                            weights=weights, strata=strata, clusters=clusters, des_ind=des_ind, 
                            bin_size=bin_size, gridpts=gridpts, perm_reps=perm_reps, pval=pval))

  
  #####################################################################################################################
  
  #--- Make frame nice ---- 
  frame <-  make_nice(frame) # format resulting tree frame for R 
  frame <- mark_ends(frame) # puts an 'E' after each end node on the frame
  
#======================================= return to original the factor levels in the data ==============================
  if(n_cats>0){
    for(i in 1:n_cats){
      #---- return to original levels -------
      levels(data[,vX[which(cat_vec)[i]]]) <- unlist(cat_table[[i]])
    }
  }

  # put original level names in the frame
  cat_splits<-which(frame$cat==1)
  
  for(i in cat_splits){
    vis <- which(vX[cat_vec]==frame$var[[i]]) #variable location in cat_table
    frame$xval[i] <- list(cat_table[[vis]][unlist(frame$xval[i])]) #replace numbers with factor names
  }
#======================================================================================================================
  

  callvals <-list(design=design, data.size=length(y), perm_reps=perm_reps, pval=pval, gridpts=gridpts, bin_size=bin_size)
  
  #------------- Get partition -------------------               
  partition<-get_partition(frame)
  row.names(partition)<-NULL
  
  #--- add N_hat estimate to partition -----
  #endnodes<-apply(covariates(partition[,"splits"], data), 1, function(x) partition$node[which(x==1)])
  endnodes<-partition$node[as.integer(covariates(partition$splits, data) %*% seq(length(partition$splits)))]
  N_hat <- sapply(partition$node, function(x) sum(weights[endnodes==x]))
  partition <- cbind(partition, N_hat)
  
  #add estimated variance of prediction to partition if no linear model
  if(length(all.vars(e_equ))==1){
    fm_ends <- match(partition$node, frame$node)
    
    #                                y_sig2                         var_y_bar  
    partition$p_se= sqrt(frame$loss[fm_ends]/partition$N_hat + as.numeric(frame$cvar[fm_ends]))
  }

  # Get tatal estimated R^2 for tree model
  r_2 <- 1 - sum(frame$loss[match(partition$node, frame$node)])/frame$loss[1]
  
  t1 <- list(callvals=callvals, rp_equ=rp_equ, e_equ=e_equ, r_2=r_2, frame=frame, partition=partition)
  class(t1)<-c("rpms")
  return(t1)
  
}
#############################################  End rpms ####################################################
