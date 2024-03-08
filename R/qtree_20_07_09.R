

get.qtree<-function(f1, n.num, start_node, title=NA, digits, s_size=TRUE){
  # takes a frame f1, split number n.num, root name, and rounding digits
  i <-which(f1$node==n.num)
  if(n.num==start_node) #first split
    if(f1$end[i]=="E"){
      val=round(f1$mean[i], digits=digits)
      return(cat("\\", "Tree [.{", title, "} ", "{", 
                 ifelse(s_size, paste0("\\", "fbox{n = ", f1$n[i], "} ", "\\", "\\ "),""),
                 "\\", "fbox{", paste(round(val, digits)), "}} ]", sep="")) 
    }
  else #first split not end
    return(cat("\\", "Tree [.{", title, "} ", 
               get.qtree(f1, n.num=2*start_node, start_node=start_node, digits=digits, s_size=s_size), 
               get.qtree(f1, n.num=2*start_node+1, start_node=start_node, digits=digits, s_size=s_size), "]", sep="")) 
  
  
  if((n.num %% 2)==0) # left hand split (LHS)
    if(f1$end[i]=="E"){ # end
      val=round(unlist(f1$value[i], use.names=FALSE), digits=digits)
      if(f1$cat[i]==0) # numeric split
        return(paste0("[.{$", f1$var[i], " \\", "leq ", 
                     round(unlist(f1$xval[i]),digits), " $} ",
                     "{", "\\", "fbox{node ", n.num, "} ", "\\", "\\ ",
                     ifelse(s_size, paste0("(n = ", f1$n[i], ")", "\\", "\\ "),""),
                     if(length(val)==1){
                       paste0("$", "\\mu = ", paste(round(val[1], digits)), "$} ]")
                     }
                     else{
                       paste0("$", "\\beta_0 = ", paste(round(val[1], digits)), "$",  "\\", "\\", 
                              "$", "\\beta_1= ", paste(round(val[2], digits)), "$} ]")
                     }
        ))
      else  #categorical split
        return(paste0("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                     paste(unlist(f1$xval[i]),collapse=","), "\\", "}$} ",
                     "{", "\\", "fbox{node ", n.num, "}", " \\", "\\ ",
                     ifelse(s_size, paste0("(n = ", f1$n[i], ")", "\\", "\\ "),""),
                     if(length(val)==1){
                       paste0("$", "\\mu = ", paste(round(val[1], digits)), "$} ]")
                     }
                     else{
                       paste0("$", "\\beta_0 = ", paste(round(val[1], digits)), "$",  "\\", "\\", 
                              "$", "\\beta_1= ", paste(round(val[2], digits)), "$} ]")
                     }
        ))
    }  
  else #LHS not end node
    if(f1$cat[i]==0)
      return(paste("[.{$", f1$var[i], " \\", "leq ", 
                   round(unlist(f1$xval[i]),digits), " $} ",
                   get.qtree(f1, n.num=2*n.num, start_node=start_node, digits=digits, s_size=s_size), 
                   get.qtree(f1, n.num=2*n.num+1, start_node=start_node, digits=digits, s_size=s_size), "]", sep="")) 
  else  #categorical
    return(paste("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                 paste(unlist(f1$xval[i]),collapse=","), "\\", "}$} ",
                 get.qtree(f1, n.num=2*n.num, start_node=start_node, digits=digits, s_size=s_size), 
                 get.qtree(f1, n.num=2*n.num+1, start_node=start_node, digits=digits, s_size=s_size), "]", sep="")) 
  
  else #right hand split (RHS)
    if(f1$end[i]=="E"){ #end
      val=round(unlist(f1$value[i], use.names=FALSE), digits=digits)
      
      if(f1$cat[i]==0) #numeric split
        return(paste0("[.{$", f1$var[i], " > ", 
                     round(unlist(f1$xval[i]),digits), " $} ",
                     "{", "\\", "fbox{node ", n.num, "} ", "\\", "\\ ",
                     ifelse(s_size, paste0("(n = ", f1$n[i], ") ", "\\", "\\ "),""),
                     if(length(val)==1){
                       paste0("$", "\\mu = ", paste(round(val[1], digits)), "$} ]") 
                     }
                     else{
                       paste0("$", "\\beta_0 = ", paste(round(val[1], digits)), "$",  " \\", "\\ ", 
                              "$",  "\\beta_1= ", paste(round(val[2], digits)), "$} ]") 
                     }
        ))
      else  #categorical split
        return(paste0("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                     paste(unlist(f1$xval[i]),collapse=","), "\\", "}$} ",
                     "{", "\\", "fbox{node ", n.num, "}", " \\", "\\ ",
                     ifelse(s_size, paste0("(n = ", f1$n[i], ")", "\\", "\\ "),""),
                     if(length(val)==1){
                       paste0("$", "\\mu = ", paste(round(val[1], digits)), "$} ]")
                     }
                     else{
                       paste0("$", "\\beta_0 = ", paste(round(val[1], digits)), "$",  "\\", "\\", 
                              "$", "\\beta_1= ", paste(round(val[2], digits)), "$} ]")
                     }
        ))
      
    }
  else #RHS not end node
    if(f1$cat[i]==0)
      return(paste("[.{$", f1$var[i], " > ", round(unlist(f1$xval[i]),digits), " $} ",
                   get.qtree(f1, n.num=2*n.num, start_node=start_node, digits=digits, s_size=s_size), 
                   get.qtree(f1, n.num=2*n.num+1, start_node=start_node, digits=digits, s_size=s_size), "]", sep=""))  
  else #RHS is categorical
    return(paste("[.{$", f1$var[i], " \\", "in ", "\\", "{", 
                 paste(unlist(f1$xval[i]), collapse=","), "\\","}$} ",
                 get.qtree(f1, n.num=2*n.num, start_node=start_node, digits=digits, s_size=s_size), 
                 get.qtree(f1, n.num=2*n.num+1, start_node=start_node, digits=digits, s_size=s_size), "]", sep="")) 
  
  return("Error in qtree")
  
}


# Returns a frame of a tree starting at node and going down 
# for use with qtree
get_subframe <- function(f1, node){
    
    #===========get new frame =============== 
    
    root_indx <- which(f1$node==node)
    
    if(is.null(root_indx)) stop(paste0(paste0("Node ", node), " is not in tree"))
    
    #internal function get_indx finds rows in frame of children of given node
    get_indx<-function(node){
      c(NULL,
        if(2*node %in% f1$node) c(which(f1$node==2*node), get_indx(2*node)),
        if((2*node+1) %in% f1$node) c(which(f1$node==(2*node+1)), get_indx(2*node+1)))
    } #end get indx
    
    #attach first node to children node
    indx<-c(root_indx, get_indx(node))
    
    # get new subframe 
    frame<-mark_ends(f1[indx,], node)
     
    return(frame)

}


################################################################################
#' qtree
#' 
#' Code to write a latex qtree plot
#' takes a rpm frame and returns latex code to produce qtree
#' uses linearize as a guide
#' Produces text code to produce tree structure in tex document
#' Requires using LaTex packages and the following commands in preamble of 
#' LaTex doc: 
#' \\usepackage\{lscape\} and 
#' \\usepackage\{tikz-qtree\}
#'
#' @param t1 rpms object created by rpms function
#' @param title string for the top node of the tree
#' @param label string used for labeling the tree figure
#' @param caption string used for caption
#' @param digits integer number of displayed digits
#' @param s_size boolean indicating whether or not to include sample size
#' @param scale numeric factor for scaling size of tree
#' @param lscape boolean to display tree in landscape mode
#' @param subnode starting node of subtree to plot
#' @examples
#' {
#' # model mean of retirement account value for households with reported 
#' # retirment account values > 0 using a binary tree while accounting for 
#' # clusterd data and sample weights.
#' 
#' s1<- which(CE$IRAX > 0)
#' r1 <-rpms(IRAX~EDUCA+AGE+BLS_URBN, data = CE[s1,],  weights=~FINLWT21, clusters=~CID) 
#' 
#' # get Latex code
#' qtree(r1)
#' 
#' }
#' @export

qtree<-function(t1, title=NULL, label=NA, caption="", digits=2, s_size=TRUE, scale=1, lscape=FALSE, subnode=1){
  
 
  if(subnode>1){
    f1<- get_subframe(t1$frame, subnode)
    if(is.null(title)) title=paste0("node ", subnode)} 
  else f1<-t1$frame
  
  if(is.null(title)) title=all.vars(t1$e_equ)[1]
  
  if(lscape) cat("\\", "begin{landscape} \n", sep="")

  cat("\\", "begin{figure}[htb] \n", sep="")
  cat("\\", "centering \n", sep="")
    
   cat("\\", "begin{tikzpicture}[scale=", scale, ", ] \n", sep="")
   cat("\\", "tikzset{every tree node/.style={align=center,anchor=north}} \n", sep="")
     get.qtree(f1, n.num=subnode, start_node=subnode, title=title, digits=digits, s_size=s_size)
     
     cat("\n")   
   cat("\\", "end{tikzpicture} \n", sep="")
    
  cat("\\", "caption{","\\", "small ", caption, "} \n", sep="")
  if(!is.na(label))
    cat("\\", "label{", label, "} \n", sep="")

  cat("\\", "end{figure} \n", sep="")
  if(lscape) cat("\\", "end{landscape} \n", sep="")
}

###################################### end qtree  ##############################################################################