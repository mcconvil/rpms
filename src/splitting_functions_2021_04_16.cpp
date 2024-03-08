// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <R_ext/Utils.h>


using namespace Rcpp;
using namespace arma;
using namespace std;


//##########################################################################################################
//                    External Function Declarations: 
//##########################################################################################################

// function defined in R_exit
void R_CheckUserInterrupt(void);

// functions defined in model_functions 
List survLm_fit(arma::colvec y, arma::mat X, arma::colvec weights); 
List survLm_model(arma::colvec y, arma::mat X, arma::colvec weights, arma::uvec strata, arma::uvec clusters);

//functions defined in variable_tests
int select_var(arma::colvec y, arma::mat mX, arma::mat vars, arma::uvec cat_vec, 
               arma::colvec weights, arma::uvec clusters, arma::uvec des_ind,
               Rcpp::StringVector vnames, arma::uword perm_reps, float pval);

//########################################################################################################


//####################################### get_set_grid  ########################################################################
// makes a grid of {0, 1} for each possible combination splitting the category values into two groups
arma::imat get_set_grid(int C){
  
  int p=(int)pow(2, C-1)-1;
  
  imat set_grid(p, C);
  
  irowvec sets(C); //reps a row of matrix
  sets.zeros();   //change as scroll through rows
  
  for(int j=0; j<p; ++j){
    
    int i=0; //start at first element
    
    // loop through sets until find first 0
    while(i<C){
      
      if(sets(i)==0){
        sets(i)=1;
        i = C;
      }
      else{
        sets(i)=0;
        ++i;
      }
    } //end while-loop
    
    set_grid.row(j)=sets;
  } // end for-loop
  
  return set_grid; 
}
//####################################### End get_set_grid  ########################################################################


//################################################ null_split ##########################################
// [[Rcpp::export]]
List null_split(){
  
  return(
    List::create(
      Named("node") =  NA_INTEGER, 
      Named("cat") =   NA_INTEGER, 
      Named("var") =   NA_STRING,
      Named("xval") =  NA_REAL,
      Named("n") =     NA_INTEGER,
      Named("loss")=   NA_REAL,
      Named("value") = NA_REAL,
      Named("cvar")=   NA_REAL,
      Named("mean")=   NA_REAL
    )
  );
  
} // end function

//############################### end null_splits ######################################################################

//############################### rbind_splits ######################################################################

// [[Rcpp::export]]
List rbind_splits(List split1, List split2) {
  
  // Environment base = Environment("package:base");
  // Function readline = base["readline"];
  // readline("> ");
  
  //------------- handling if one of the splits is the null-split ---------------------------
  if(any(is_na(as<NumericVector>(split1["node"])))){
    if(any(is_na(as<NumericVector>(split2["node"])))){
      List split0=null_split();
      return(split0);
    } // split1 is empty but split2 is not
    else{return(split2);}
  } // if split1 not empty
  else{if(any(is_na(as<NumericVector>(split2["node"])))){return(split1);}}
  //------------------- done handling null ------------------------------------------------- 
  
  arma::vec node = join_vert(as<arma::vec>(split1["node"]), as<arma::vec>(split2["node"]));
  arma::vec cat = join_vert(as<arma::vec>(split1["cat"]), as<arma::vec>(split2["cat"]));
  
  //------------- get var -------------------------
  CharacterVector var1 =  split1["var"];
  CharacterVector var2 =  split2["var"];
  var1 = as< vector< string > >(var1);
  var2 = as< vector< string > >(var2);
  for (int i=0;i<var2.length();i++) {
    var1.push_back(var2(i));
  }
  
  //-------------- get xval ------------------
  vector< vec > xval1=split1["xval"];
  vector< vec > xval2=split2["xval"];
  for (uword i=0;i<xval2.size();i++) {
    xval1.push_back(xval2.at(i));
  }
  
  arma::vec n = join_vert(as<arma::vec>(split1["n"]), as<arma::vec>(split2["n"]));
  
  arma::vec loss = join_vert(as<arma::vec>(split1["loss"]), as<arma::vec>(split2["loss"]));
  
  //---------- get value -----------------------
  vector< vec > val1=split1["value"];
  vector< vec > val2=split2["value"];
  for (uword i=0;i<val2.size();i++) {
    val1.push_back(val2.at(i));
  }
  
  //-------------- get cvar -----------------------------
  vector< vec > cvar1=split1["cvar"];
  vector< vec > cvar2=split2["cvar"];
  for (uword i=0;i<cvar2.size();i++) {
    cvar1.push_back(cvar2.at(i));
  }
  
  arma::vec mean = join_vert(as<arma::vec>(split1["mean"]), as<arma::vec>(split2["mean"]));
  
  return List::create(
    Named("node") =  node, 
    Named("cat") =   cat, 
    Named("var") =   var1,
    Named("xval") =  xval1,
    Named("n") =     n,
    Named("loss")=   loss,
    Named("value") = val1,
    Named("cvar")=   cvar1,
    Named("mean")=   mean
  );
} // end function

//############################### end rbind_splits ######################################################################


//####################################### get_node  ########################################################################
// [[Rcpp::export]]
List get_node(arma::uword node, int cat, std::string vname, arma::colvec y, arma::colvec weights, arma::vec mxval, arma::uvec s, List modfit){
  
  
  if(node==1) { s=find(y<datum::inf); } //all elements in y

  arma::colvec resid=modfit["residuals"];
  double loss = dot(weights(s)%resid, resid); // weighted squared residuals
  //double loss = sum(resid.t()*resid); // squared residuals
  vector <rowvec> xval(1);
  xval.at(0)=mxval.t();
  //cout<<"Got xval"<<endl;
  
  arma::rowvec value1 = modfit["coefficients"];
  vector <rowvec> value;
  value.push_back(value1);
  //cout<<"Got value1"<<endl;
  
  arma::mat spL_covM=modfit["covM"];
  arma::rowvec  cvar1 = spL_covM.diag().t();
  vector <rowvec> cvar;
  cvar.push_back(cvar1);
  //cout<<"Got cvar"<<endl;
  
  double spL_mean = mean(y(s)-resid);
  
  return List::create(
    Named("node") =  node, 
    Named("cat") =   cat, 
    Named("var") =   vname,
    Named("xval") =  xval,
    Named("n") =     s.n_elem,
    Named("loss")=   loss,
    Named("value") = value,
    Named("cvar")=   cvar,
    Named("mean")=   spL_mean
  ); //end create
  
}//end function

//####################################### End get_node  ########################################################################



//############################################################################################################
//
//            Definitions of loss functions one for numeric and one for categorical
//
//##############################################################################################################

// [[Rcpp::export]]
double get_loss(double x, arma::vec &x_val,  arma::colvec &y, arma::mat &mX, arma::colvec &weights) {
  
  uvec s=find(x_val<= x);
  uvec not_s=find(x_val> x);
  
  arma::colvec res_s = survLm_fit(y.elem(s), mX.rows(s), weights.elem(s))["residuals"];
  arma::colvec res_ns = survLm_fit(y.elem(not_s), mX.rows(not_s), weights.elem(not_s))["residuals"];
  
  return(dot(weights.elem(s)%res_s, res_s) + dot(weights.elem(not_s)%res_ns,res_ns));
}


//######################################### small_loss ############################################################################
//------ small_loss function returns the x_value on which to split and the loss associtated 
//after that split where the loss is small

// [[Rcpp::export]]
List small_loss(arma::vec &x_val,  arma::vec &uq_xs, arma::colvec &y, arma::mat &mX, 
                arma::colvec &weights, double Lx, double Ux, double Lloss, double Uloss, arma::uword gridpts) {

  //cout<<"in small_loss "<<endl;
  R_CheckUserInterrupt();
  uvec midxs_I = find(uq_xs>Lx && uq_xs<Ux); //index of unique x vaulues in the middle
  uword mid_n = midxs_I.n_elem; 
  
  // --- middle contains at least 3 values do recursion until no more than three points between  
  if(mid_n > gridpts){

    rowvec loss(gridpts+2), xs(gridpts+2); // vector for loss and x_values checked
    loss.fill(datum::inf);
    loss(0)=Lloss;
    loss(gridpts+1)=Uloss;
    xs(0)=Lx;
    xs(gridpts+1)=Ux;
    
    double delta=(mid_n/(gridpts+1));

    for(uword i=1; i<=gridpts; ++i){
      
      xs(i) = uq_xs(midxs_I((uword)round((float)(i*delta))-1)); //next element -1 because midxs_I begins at 0
      loss(i)=get_loss(xs(i), x_val, y, mX, weights);

    } //end loss for loop

    arma::uword min_indx=loss.index_min();  //get index of min loss
    
    //do recursion
    if(min_indx==0){return(small_loss(x_val,  uq_xs, y, mX, weights,
       xs(min_indx), xs(min_indx+1), loss(min_indx), loss(min_indx+1), gridpts));} else 
         
         if(min_indx==(gridpts+1)){return(small_loss(x_val,  uq_xs, y, mX, weights,
            xs(min_indx-1), xs(min_indx), loss(min_indx-1), loss(min_indx), gridpts));} else{
              double Ux, Lx, Uloss, Lloss;
              Ux=uq_xs(midxs_I((uword)round((float)((min_indx+.5)*delta))-1));  //element midxs_I begins at 0
              Lx=uq_xs(midxs_I((uword)round((float)((min_indx-.5)*delta))-1)); //element midxs_I begins at 0
              Lloss=get_loss(Lx, x_val, y, mX, weights);
              Uloss=get_loss(Ux, x_val, y, mX, weights);  
              return(small_loss(x_val,  uq_xs, y, mX, weights, Lx, Ux, Lloss, Uloss, gridpts));
              
            } //end else min is a middle point
        
  } //end if more than gridpts x values in the middle or delta<2

  // mid_n =< gridpts x-values in the middle, so return smallest loss and split point
  rowvec loss(mid_n+2), xs(mid_n+2); // vectors for loss and x-values
  loss.fill(datum::inf);
  loss(0)=Lloss;
  loss(mid_n+1)=Uloss;
  xs(0)=Lx;
  xs(mid_n+1)=Ux;

  //find loss and xs for each point between Lx and Ux   
  for(uword i=1; i<=mid_n; ++i){
    xs(i)=uq_xs(midxs_I(i-1)); // mdixs_I is index of mid_n unique values between Lx and Ux
    loss(i)=get_loss(xs(i), x_val, y, mX, weights);
  } //end loss for loop

  arma::uword min_indx=loss.index_min();
  uvec min_pt = find(uq_xs==xs(min_indx));

 double xvalue = (uq_xs(min_pt(0)) + uq_xs(min_pt(0)+1))/2;

  
  return List::create(Rcpp::Named("xval") = xvalue,
                    Rcpp::Named("loss")  = loss(min_indx));
  
  // return List::create(Rcpp::Named("xval") = xs(min_indx),
  //                     Rcpp::Named("loss")  = loss(min_indx));
  
}
//******************************* end get_loss numeric *****************************************************

double get_loss_cat(uword xs, arma::imat &sets, arma::vec &cats, arma::vec &x_val, arma::colvec &y, arma::mat &mX,
                    arma::colvec &weights){
 
 R_CheckUserInterrupt();
  
  uvec s, not_s;
  
  uvec lsets=find(sets.row(xs)==1); //find categories on left side
  for(uword c=0; c<lsets.n_elem; ++c){ //loop through categories to get s
    s= join_cols(s, find(x_val==cats(lsets(c))));
  }
  s=vectorise(s);
  
  not_s.reset();
  lsets.reset();
  lsets=find(sets.row(xs)==0); //which categories not in left side group
  for(uword c=0; c<lsets.n_elem; ++c){ //loop through each category not in left group to get not_s
    not_s= join_cols(not_s, find(x_val==cats(lsets(c)))); // which variables have cats
  }
  not_s=vectorise(not_s);
  
  //------------------------------------------------------------------------------------------------------------      
  arma::colvec res_s = survLm_fit(y.elem(s), mX.rows(s), weights.elem(s))["residuals"];
  arma::colvec res_ns = survLm_fit(y.elem(not_s), mX.rows(not_s), weights.elem(not_s))["residuals"];
  double loss=dot(weights.elem(s)%res_s, res_s) + dot(weights.elem(not_s)%res_ns,res_ns);
  
  return(loss);
}


//********************************* small_loss_cat categorical  *****************************************************
// [[Rcpp::export]]
List small_loss_cat(arma::imat &sets, arma::vec &cats, arma::vec &x_val, arma::colvec &y, arma::mat &mX,
                    arma::colvec &weights, std::vector<arma::uvec> cat_ind, 
                    arma::uword Lx, arma::uword Ux, double Lloss, double Uloss, arma::uword gridpts) {

  R_CheckUserInterrupt(); 
  uword mid_n;
  if(Ux==Lx) mid_n= 0; else mid_n=((Ux-Lx)-1); // number of rows between Ux and Lx
  
  if(mX.n_cols > 1) gridpts=sets.n_rows;

  // --- middle contains at least gridpts values, do recursion until no more than gridpts points between ends  
  if(mid_n>gridpts){
    
    rowvec loss((gridpts+2)); // vector for loss value at gridpts points in the middle plus Lx and Ux
    loss.fill(datum::inf);
    loss(0)=Lloss;          // first element of loss is Lloss
    loss(gridpts+1)=Uloss; //last element of loss is Uloss

    urowvec xs(gridpts+2);  // index of rows of sets to consider
    xs(0)=Lx;          // first point is Lx
    xs(gridpts+1)=Ux;  // last point is Ux
    float delta=mid_n/(gridpts+1);  // get difference between each point
    
    //find loss for each midpoint checked
 
    for(uword i=1; i<=gridpts; ++i){
      
      xs(i) = Lx+(uword)round(i*delta);  //next row of sets 
      loss(i)=get_loss_cat(xs(i), sets, cats, x_val, y, mX, weights);

    } //end loss for-loop

   arma::uword min_indx=loss.index_min();  //get index of min loss

    //do recursion.... 
    if(min_indx==0){return(small_loss_cat(sets, cats, x_val, y, mX, weights, cat_ind, 
       xs(min_indx), xs(min_indx+1), loss(min_indx), loss(min_indx+1), gridpts));} else
         
         if(min_indx==(gridpts+1)){return(small_loss_cat(sets, cats, x_val, y, mX, weights, cat_ind, 
            xs(min_indx-1), xs(min_indx), loss(min_indx-1), loss(min_indx), gridpts));} else{
              
              //mid point is the min
              double Ux, Lx, Uloss, Lloss;
              Lx=(uword)round((float)((min_indx-.5)*delta));
              Ux=(uword)round((float)((min_indx+.5)*delta));
              Lloss=get_loss_cat(Lx, sets, cats, x_val, y, mX, weights);
              Uloss=get_loss_cat(Ux, sets, cats, x_val, y, mX, weights); 
              return(small_loss_cat(sets, cats, x_val, y, mX, weights, cat_ind, Lx, Ux, Lloss, Uloss, gridpts));
            } //end else min is a middle point
         
         //else 
           //   return(small_loss_cat(sets, cats, x_val, y, mX, weights, cat_ind, 
            //             xs(min_indx-1), xs(min_indx+1), loss(min_indx-1), loss(min_indx+1), gridpts));
            
  } //end if more than gridpts rows of sets in the middle
  
  // -------- out of recursion ---------------
  // less than gridpts rows of sets in the middle, so return smallest loss and split point
  
  rowvec loss(mid_n+2);  //loss vector length is number of middle points plus Lx an Ux
  loss.fill(datum::inf);
  loss(0)=Lloss;
  loss(mid_n+1)=Uloss;
  urowvec xs(mid_n+2);
  xs(0)=Lx;        //vector of rows of set.table under consideration
  xs(mid_n+1)=Ux;

  //for each row of sets between Lx and Ux get loss and row number   
  for(uword i=1; i<=mid_n; ++i){
    
    xs(i) = Lx+i; //row number
    loss(i)=get_loss_cat(xs(i), sets, cats, x_val, y, mX, weights);

  } //end loss for loop

 arma::uword min_indx=loss.index_min();
  
 // cout<<"xs is "<<xs<<" and min_indx is "<<min_indx<<endl;
//  cout<<"loss is "<<loss<<" at x-value"<<xs(min_indx)<<endl;
  
  return List::create(Rcpp::Named("xval") = xs(min_indx),
                      Rcpp::Named("loss")  = loss(min_indx));
  
 
}  //end middle contains less than gridpoint values;
  

//################################### end small_loss_cat functions  ##################################################################


//####################################################################################################################################
//
//                                                        split_rpms 
//
//  main split function (uses mean squared error for loss function)
//####################################################################################################################################
// [[Rcpp::export]]
List split_rpms(arma::uword node, arma::colvec y, arma::mat mX, arma::mat X, Rcpp::StringVector vnames, arma::uvec cat_vec, 
                arma::colvec weights, arma::uvec strata, arma::uvec clusters, arma::uvec des_ind,  
                arma::uword &bin_size, arma::uword &gridpts, arma::uword &perm_reps, float &pval){
  
  //C++ check for user interupt
  R_CheckUserInterrupt();
  
  arma::uword n = y.n_elem;
  
  if(n<=2*bin_size || n < mX.n_cols) return(null_split());
  
  //cout<<"going to select a variable"<<endl;
  int x = select_var(y, mX, X, cat_vec, weights, clusters, des_ind, vnames, perm_reps, pval);
  if(x ==-1) return(null_split()); //no significant split

  //*********************************** Get split point and loss **************************************** 
  
  arma::vec x_vals=X.col(x); // vector of values to split on
  arma::vec uvar= unique(x_vals); // sorted unique x values
  arma::uword cat=cat_vec.at(x); // indicator = categorical variable 1, or not 0
  std::string var=as<string>(vnames.at(x)); //name of variable 

  //#------------------------ numeric variable  -------------------------------
  if(cat==0) {

    // ------------------ find Lx, Lloss and Ux, Uloss to begin process -----------------
    
    arma::vec svals=sort(X.col(x));//sorted and original x values
    uvec s, not_s;
    
    // -------- Find  Lx and Ux ------------------------------
    s=find(svals>svals(bin_size-1)); //the (bin_size - 1)th largest value because vector starts at 0
    if(s.n_elem < 1 ) return(null_split());
    double Lx=min(svals(s));
    double Lloss=get_loss(Lx, x_vals, y, mX, weights);
    
    s=find(svals<svals(n -bin_size)); //the (n-bin_size)th largest value because vector starts at 0
    if(s.n_elem < 1 ) return(null_split());
    double Ux=max(svals(s));
    double Uloss=get_loss(Ux, x_vals, y, mX, weights);

    //---------------------- Got Lx, Lloss, and Ux, Uloss  -------------------------------------------------------------------------------
    
    //--- get loss and xval
    List sp_point = small_loss(x_vals, uvar, y, mX, weights, Lx, Ux, Lloss, Uloss, gridpts);

    double loss= sp_point["loss"]; 
    if(loss == datum::inf){ return(null_split()); }  //end if no split
    
    arma::rowvec mxval(1); 
    mxval(0)=sp_point["xval"]; //value of x rep as vector of size one

    //----- code to split at .5 if integer
    //   if(all(floor(X.col(x))==ceil(X.col(x)))) mxval.fill(uvar.at(min_idx));
    //  else  mxval.fill((uvar.at(min_idx) + uvar.at(min_idx+1))/2);
    
    //------- the Left node ---------------------------
    arma::uvec sL = find(x_vals<=mxval.at(0));
    
    List modfit = survLm_model(y.elem(sL), mX.rows(sL), weights.elem(sL), strata.elem(sL), clusters.elem(sL));
    
    List SpL=get_node(2*node, cat=0, var, y, weights, mxval, sL, modfit);
    List SpL2=split_rpms(SpL["node"], y(sL), mX.rows(sL), X.rows(sL), vnames, cat_vec,
                         weights(sL), strata(sL), clusters(sL), des_ind, bin_size, gridpts, perm_reps, pval);
    //-------------------------------------------------------------------------------------------------------
    
    //------- the Right node -----------------------------------------------------------------------------
    uvec sR = find(x_vals>mxval.at(0));
    
    modfit = survLm_model(y.elem(sR), mX.rows(sR), weights.elem(sR), strata.elem(sR), clusters.elem(sR));
    
    List SpR=get_node(2*node+1, cat=0, var, y, weights, mxval, sR, modfit);
    List SpR2=split_rpms(SpR["node"], y(sR), mX.rows(sR), X.rows(sR), vnames, cat_vec,
                         weights(sR), strata(sR), clusters(sR), des_ind, bin_size, gridpts, perm_reps, pval);
    
    //---------------------------------------------------------------------------------------------
  
    return rbind_splits(rbind_splits(SpL, SpL2), rbind_splits(SpR, SpR2)); 
    
  } // end numeric
  
  
  //#------------------------ categorical variable  -------------------------------  
  if(cat==1){
    //uvar is "x.cats" unique set of categories
    arma::uword p = uvar.n_elem;  

    if(p<2) return(null_split());  //end if no split
    
    // get vector of uvecs for each category
    vector<uvec> cat_ind(p); //vector to hold indexs for each category
    arma::uvec cat_size(p);
    
    for(uword c=0; c<p; ++c){
      uvec ci = find(X.col(x)==uvar(c)); //find index of each category label 
      cat_ind.at(c)= ci;
      cat_size(c)=ci.n_elem;
    }
    // got the vector of indices for each category
    
    imat power_set;
    
    // ----- check if doing means then do linear order of categories ------
    if(mX.n_cols==1){
      
      power_set=imat(p-1, p);
      power_set.zeros();
      
      vec cat_means(p);
      cat_means.fill(datum::inf);
      
      //for each category find the mean y value
      for(uword i=0; i<p; ++i){
        uvec cs = cat_ind.at(i);
       // changed this to be the survey weighted average
        cat_means(i) = sum(y(cs)%weights(cs))/sum(weights(cs));
      } 
      
      //get matrix of 1-0 whether cat included in left side or not
      // start at smallest and include that in the left then next smallest etc..
      for(uword i=0; i<(p-1); ++i){
        arma::uword cmin = cat_means.index_min();
        power_set.rows(i, p-2).col(cmin).fill(1); //all rows below 
        cat_means(cmin) = datum::inf;  //set min to inf then go find next smallest mean
      }
      
    } //end if doing means
    
    else{  //not doing means
       if(p>14) stop("Categorical variable ", x, " has > than 14 categories.  This would take a long time when fitting models on each node.");
      power_set = get_set_grid(p);
    } 
    
    uvec s, not_s;
    
    // -------- Find size (number of observations) for each rows of power_set-------
    uvec row_size(power_set.n_rows);
    for(uword i=0; i<power_set.n_rows; ++i){
      s=find(power_set.row(i)==1);
      row_size(i)=sum(cat_size(s));
    }
    
    //--- set power_set to only include usable rows ----
    s=find(row_size>bin_size && row_size < n-bin_size);
    if(s.n_elem>0) power_set=power_set.rows(s); else return(null_split()); 
    
    
    // ------------------ find Lx, Lloss and Ux, Uloss to begin process -----------------
    
    // -------- Find  Lx and Lloss ------------------------------
    arma::uword Lx=0; // index of first row in power set
    double Lloss=get_loss_cat(Lx, power_set, uvar, x_vals, y, mX, weights);
    
    //----------------------------------------------------------------- 
    
    //---------- Find Ux and Uloss --------------------------------------
    arma::uword Ux=power_set.n_rows-1;  //index of last row in power set
    double Uloss=get_loss_cat(Ux, power_set, uvar, x_vals, y, mX, weights);
 
    //---------------------- Got Lx, Lloss, and Ux, Uloss  --------------------------------------------------------------------------------  
    
    //--- get loss and xval; xval is index of row of reduced power set that provides the min
    List sp_point = small_loss_cat(power_set, uvar, x_vals, y, mX, weights, cat_ind, Lx, Ux, Lloss, Uloss, gridpts);
    
    double loss= sp_point["loss"];
    if(loss == datum::inf){ return(null_split()); }  //end if no split
    
    
    arma::uword min_idx=sp_point["xval"]; //value of x is row of reduced power set to use
    
    
    //------- the Left node ---------------------------
    
    uvec left_cats = find(power_set.row(min_idx)==1);
    arma::uword ncats_L = left_cats.n_elem;
    vec mxval_L = uvar(left_cats);

    //get index vector of all data with category left side
    uvec sL = cat_ind.at(left_cats(0));
    for(uword i=1; i<ncats_L; ++i) sL= join_vert(sL,cat_ind.at(left_cats(i)));
    
    List modfit = survLm_model(y.elem(sL), mX.rows(sL), weights.elem(sL), strata.elem(sL), clusters.elem(sL));
    
    List SpL=get_node(2*node, cat=1, var, y, weights, mxval_L, sL, modfit);
    
    List SpL2=split_rpms(SpL["node"], y(sL), mX.rows(sL), X.rows(sL), vnames, cat_vec,
                         weights(sL), strata(sL), clusters(sL), des_ind, bin_size, gridpts, perm_reps, pval);
    //-------------------------------------------------------------------------------------------------------
    
    //------- the Right node -----------------------------------------------------------------------------
    
    uvec right_cats = find(power_set.row(min_idx)==0);
    //  if(right_cats.n_elem<1)cout<<"########### right_cats is "<<right_cats<< " "<<endl; 
    arma::uword ncats_R = right_cats.n_elem;
    
    vec mxval_R = uvar(right_cats);

    
    //get index vector of all data with category right side
    uvec sR = cat_ind.at(right_cats(0));
    for(uword i=1; i<ncats_R; ++i) sR= join_vert(sR, cat_ind.at(right_cats(i)));
    
    modfit = survLm_model(y.elem(sR), mX.rows(sR), weights.elem(sR), strata.elem(sR), clusters.elem(sR));
    
    List SpR=get_node(2*node+1, cat=1, var, y, weights, mxval_R, sR, modfit);
    
    List SpR2=split_rpms(SpR["node"], y(sR), mX.rows(sR), X.rows(sR), vnames, cat_vec,
                         weights(sR), strata(sR), clusters(sR), des_ind, bin_size, gridpts, perm_reps, pval);
    
    //---------------------------------------------------------------------------------------------
     // cout<<"leaving rpms_split successfully"<<endl;
    return rbind_splits(rbind_splits(SpL, SpL2), rbind_splits(SpR, SpR2)); 
    
    
  } //end categorical

  return(null_split()); //if neither numerical or categorical 
  
} //end rpms_split

//################################# End rpms_split ##############################################################


//                                                 NEW

//####################################################################################################################################
//
//                                                        random_split 
//
//  Chooses best variable and split point among pV randomly selected splits---- for use with forest 
//####################################################################################################################################
// [[Rcpp::export]]
List random_split2(arma::uword node, arma::colvec y, arma::mat mX, arma::mat X, Rcpp::StringVector vnames, arma::uvec cat_vec, 
                arma::colvec weights, arma::uvec strata, arma::uvec clusters, arma::uvec des_ind,  
                arma::uword &bin_size){
  
  //C++ check for user interupt
  R_CheckUserInterrupt();

  arma::uword n = y.n_elem;
    

// probability of splitting each time is between minp and 100 percent
// if(node>1){
//   double minp=.98;
//   double psplit=minp + (1-minp)*randu(), Runif=randu();
//   if(psplit <Runif) return(null_split());
// }


  // Random chance to stop splitting
  // chance of spliting at each node starting at root
  // 0.98 0.92 0.87 0.82 0.77 0.73 0.68 if .98 multiply by 3
  // 0.99 0.94 0.89 0.85 0.81 0.77 0.73 if .99 muliply by 5
  
  //  double alpha = (int)floor((float)log2(node)), p_split = pow(.99, 5*alpha+1), Runif=randu();
  //if(p_split<Runif) return(null_split());
  
  // check if reached bin-size limit

  if(n<=2*bin_size || n < mX.n_cols) return(null_split());

  //Get a count of the number of unique values each variable has
  uword pV=X.n_cols;

  uvec uvar_count(pV);
  for(uword i=0; i<pV; ++i){
    vec unique_vals=unique(X.col(i));
    uvar_count(i)=unique_vals.n_elem;
  }
  
  
  //Randomly select variable that has more than one value
  uvec Vindx=find(uvar_count>1);  //Vindx is an index of variables with more than one value
  
  if(Vindx.n_elem==0){return(null_split());}


  // p-vectors for randomly selected variables and values and corresponding loss
  uvec x_vec(pV);
  arma::vec loss_vec(pV);
  arma::rowvec xval_vec(pV);  //1x1 matrix for split value
  
  //used for getting ranodm xs and x-values and corresponding losses
  uvec Rindx;  //to get random x
  arma::vec x_vals, svals, uvals;
  arma::uword cat;
  uvec s, not_s, rand_s;

  
  // get pV random splits
  for(uword xindx=0; xindx<pV; ++xindx){
  Rindx=shuffle(Vindx);  
  x_vec(xindx) = Rindx(0); // variable selected randomly   
  x_vals=X.col(x_vec(xindx));
  svals= sort(x_vals); // sorted unique x values
  uvals= unique(x_vals); // sorted unique x values
  cat=cat_vec.at(x_vec(xindx));


  //*********************************** Get split point and loss **************************************** 
  
  //#------------------------ numeric variable  -------------------------------
  if(cat==0) {

    // ------------------ find random split value of x between edges 
   // arma::vec svals=sort(X.col(x));
    
    s=find(svals>svals(bin_size-1) && svals<svals(n-bin_size)); //index of xvalues that are split-able
    if(s.n_elem < 1 ) loss_vec(xindx)=datum::inf;
    else{
    rand_s=shuffle(s); // random index off svals
    
   // get midpoint between of unique values
    uvec indx = find(uvals==svals(rand_s(0)));
    xval_vec(xindx)=(uvals(indx(0))+uvals(indx(0)+1))/2;
    //xval(0)=svals(rand_s(0)); //value of x rep as vector of size one
  
    //--------- Find loss ----------------------------------
    s=find(x_vals<=xval_vec(xindx));
    not_s = find(X.col(x_vec(xindx))>xval_vec(xindx));

    loss_vec(xindx)=get_loss(xval_vec(xindx), x_vals, y, mX, weights);
    
    } //end else 
  } // end numeric
  
  if(cat==1){
    
    //uvar is "x.cats" unique set of categories
    arma::uword p = uvals.n_elem;  
    
    if(p<2) return(null_split());  //end if no split
    
    // get vector of uvecs for each category
    vector<uvec> cat_ind(p); //vector to hold indexs for each category
    arma::uvec cat_size(p);
    
    for(uword c=0; c<p; ++c){
      uvec ci = find(X.col(x_vec(xindx))==uvals(c)); //find index of each category label 
      cat_ind.at(c)= ci;
      cat_size(c)=ci.n_elem;
    }
    // got the vector of indices for each category
    
    imat power_set;
    
    // ----- forests always do means so we do linear order of categories ------

      power_set=imat(p-1, p);
      power_set.zeros();
      
      vec cat_means(p);
      cat_means.fill(datum::inf);
      
      for(uword i=0; i<p; ++i){
        //  uvec cs = find(X.col(x)==uvar(i));
        uvec cs = cat_ind.at(i);
        // cat_means(i) = mean(y(cs));  // changed this to be the survey weighted average
        cat_means(i) = sum(y(cs)%weights(cs))/sum(weights(cs));
      } 
      
      for(uword i=0; i<(p-1); ++i){
        arma::uword cmin = cat_means.index_min();
        power_set.rows(i, p-2).col(cmin).fill(1); //all rows below 
        cat_means(cmin) = datum::inf;  //set min to inf then go find next smallest mean
      }
      
    
    // -------- Find  rows of power_set that are useable------------------------------
    uvec row_size(power_set.n_rows);
    for(uword i=0; i<power_set.n_rows; ++i){
      s=find(power_set.row(i)==1);
      row_size(i)=sum(cat_size(s));
    }
    
    //--- set power_set to only include usable rows ----
    s=find(row_size>bin_size && row_size < n-bin_size);
    if(s.n_elem==0) {loss_vec(xindx)=datum::inf;}
    // ------------- got powwer set --------------------------------------------------
    
    else{
    // ------------------ get split value = xval and loss at xval = loss ---------------------------------------
    
    // get random split
    rand_s=shuffle(s);
    xval_vec(xindx) = rand_s(0); 
    
    s.reset();
    uvec lcats=find(power_set.row(xval_vec(xindx))==1);
    for(uword c=0; c<lcats.n_elem; ++c){ //loop to get s
      s= join_cols(s, find(X.col(x_vec(xindx))==uvals(lcats(c))));
    }
    s=vectorise(s);
    
    uvec rcats=find(power_set.row(xval_vec(xindx))==0);
    for(uword c=0; c<rcats.n_elem; ++c){ //loop to get not_s
      not_s= join_cols(not_s, find(X.col(x_vec(xindx))==uvals(rcats(c))));
    }
    not_s=vectorise(not_s);
    
    // if(s.n_elem < 1 || not_s.n_elem < 1){ return(null_split()); } assume this can't happen  
    
    // get residuals
    arma::colvec res_s = survLm_fit(y.elem(s), mX.rows(s), weights.elem(s))["residuals"];
    arma::colvec res_ns = survLm_fit(y.elem(not_s), mX.rows(not_s), weights.elem(not_s))["residuals"];
    loss_vec(xindx)=dot(weights.elem(s)%res_s, res_s) + dot(weights.elem(not_s)%res_ns,res_ns);
    
    //---------------------- got xval and loss ------------------------------------------- 
    } //end else
  } // end categorical
        
  } //end i loop for number of randoms to try

  // cout<<"x_vec"<<endl<<x_vec<<endl;
  // cout<<"xval_vec"<<endl<<xval_vec<<endl;
  // cout<<"loss_vec"<<endl<<loss_vec<<endl;
  
  uword mindx = index_min(loss_vec);
  uword x = x_vec(mindx);
  std::string var=as<string>(vnames.at(x));
  

  arma::vec xval(1);
  xval(0)= xval_vec(mindx);

  
  double loss = loss_vec(mindx);
  if(loss == datum::inf){ return(null_split()); }
  
// cout<<"mindx"<<endl<<mindx<<endl;
// cout<<"x"<<endl<<x<<endl;
// cout<<"var"<<endl<<var<<endl;


// List L1;
// 
// return(L1);

//return(null_split());
 //================================================= Get nodes =======================================================

 
 x_vals=X.col(x);
 svals= sort(x_vals); // sorted unique x values
 uvals= unique(x_vals); // sorted unique x values
 cat=cat_vec.at(x);
 var=as<string>(vnames.at(x));
 
 
 //#------------------------ numeric variable  -------------------------------
 if(cat==0) {
   
   //--------- Find loss ----------------------------------
   s=find(x_vals<=xval(0));
   not_s = find(X.col(x)>xval(0));
   if(s.n_elem < 1 || not_s.n_elem < 1){ return(null_split()); }  
   
   arma::colvec res_s = survLm_fit(y.elem(s), mX.rows(s), weights.elem(s))["residuals"];
   arma::colvec res_ns = survLm_fit(y.elem(not_s), mX.rows(not_s), weights.elem(not_s))["residuals"];
   double loss=get_loss(xval(0), x_vals, y, mX, weights);
   if(loss == datum::inf){ return(null_split()); }  //end if no split
   
  
   //================================================= Get nodes =======================================================
   //------- the Left node ---------------------------
   List modfit = survLm_model(y.elem(s), mX.rows(s), weights.elem(s), strata.elem(s), clusters.elem(s));
   
   List SpL=get_node(2*node, cat=0, var, y, weights, xval, s, modfit);
   List SpL2=random_split2(SpL["node"], y(s), mX.rows(s), X.rows(s), vnames, cat_vec,
                          weights(s), strata(s), clusters(s), des_ind, bin_size);
   //-------------------------------------------------------------------------------------------------------
   
   //------- the Right node -----------------------------------------------------------------------------
   // sR is not_s
   modfit = survLm_model(y.elem(not_s), mX.rows(not_s), weights.elem(not_s), strata.elem(not_s), clusters.elem(not_s));
   
   List SpR=get_node(2*node+1, cat=0, var, y, weights, xval, not_s, modfit);
   List SpR2=random_split2(SpR["node"], y(not_s), mX.rows(not_s), X.rows(not_s), vnames, cat_vec,
                          weights(not_s), strata(not_s), clusters(not_s), des_ind, bin_size);
   //---------------------------------------------------------------------------------------------
   
   return rbind_splits(rbind_splits(SpL, SpL2), rbind_splits(SpR, SpR2)); 
   
 } // end numeric
 
 //------------------------------------- categorical variable  -------------------------------  
 if(cat==1){
   //uvar is "x.cats" unique set of categories
   arma::uword p = uvals.n_elem;  
   
   if(p<2) return(null_split());  //end if no split
   
   // get vector of uvecs for each category
   vector<uvec> cat_ind(p); //vector to hold indexs for each category
   arma::uvec cat_size(p);
   
   for(uword c=0; c<p; ++c){
     uvec ci = find(X.col(x)==uvals(c)); //find index of each category label 
     cat_ind.at(c)= ci;
     cat_size(c)=ci.n_elem;
   }
   // got the vector of indices for each category
   
   imat power_set;
   
   // ----- check if doing means then do linear order of categories ------
   if(mX.n_cols==1){
     
     power_set=imat(p-1, p);
     power_set.zeros();
     
     vec cat_means(p);
     cat_means.fill(datum::inf);
     
     for(uword i=0; i<p; ++i){
       //  uvec cs = find(X.col(x)==uvar(i));
       uvec cs = cat_ind.at(i);
       // cat_means(i) = mean(y(cs));  // changed this to be the survey weighted average
       cat_means(i) = sum(y(cs)%weights(cs))/sum(weights(cs));
     } 
     
     for(uword i=0; i<(p-1); ++i){
       arma::uword cmin = cat_means.index_min();
       power_set.rows(i, p-2).col(cmin).fill(1); //all rows below 
       cat_means(cmin) = datum::inf;  //set min to inf then go find next smallest mean
     }
     
   } //end if doing means
   
   else{  //not doing means
     // if(p>12) warning("Categorical variable ", x, " has > than 12 categories.  This may take a long time when fitting models on each node.");
     power_set = get_set_grid(p);
   } 
   //cout<<power_set;
   
   
   uvec s, not_s, rand_s;
   
   // -------- Find  rows of power_set that are useable------------------------------
   uvec row_size(power_set.n_rows);
   for(uword i=0; i<power_set.n_rows; ++i){
     s=find(power_set.row(i)==1);
     row_size(i)=sum(cat_size(s));
   }
   
   //--- set power_set to only include usable rows ----
   s=find(row_size>bin_size && row_size < n-bin_size);
   if(s.n_elem==0) return(null_split()); 
   // ------------- got powwer set --------------------------------------------------
   
   // ------------------ get split value = xval and loss at xval = loss ---------------------------------------
   
   // get random split
   
   s.reset();
   uvec lcats=find(power_set.row(xval(0))==1);
   for(uword c=0; c<lcats.n_elem; ++c){ //loop to get s
     s= join_cols(s, find(X.col(x)==uvals(lcats(c))));
   }
   s=vectorise(s);
   
   uvec rcats=find(power_set.row(xval(0))==0);
   for(uword c=0; c<rcats.n_elem; ++c){ //loop to get not_s
     not_s= join_cols(not_s, find(X.col(x)==uvals(rcats(c))));
   }
   not_s=vectorise(not_s);
   
   if(s.n_elem < 1 || not_s.n_elem < 1){ return(null_split()); }  
   
   // get residuals
   arma::colvec res_s = survLm_fit(y.elem(s), mX.rows(s), weights.elem(s))["residuals"];
   arma::colvec res_ns = survLm_fit(y.elem(not_s), mX.rows(not_s), weights.elem(not_s))["residuals"];
   //    double loss=dot(weights.elem(s)%res_s, res_s) + dot(weights.elem(not_s)%res_ns,res_ns);
   
   
   //------- the Left node ---------------------------
   
   // s= sL   uvec left_cats = find(power_set.row(xval)==1);
   arma::uword ncats_L = lcats.n_elem;
   vec mxval_L = uvals(lcats);
   
   //get index vector of all data with category left side
   
   List modfit = survLm_model(y.elem(s), mX.rows(s), weights.elem(s), strata.elem(s), clusters.elem(s));
   List SpL=get_node(2*node, cat=1, var, y, weights, mxval_L, s, modfit);
   List SpL2=random_split2(SpL["node"], y(s), mX.rows(s), X.rows(s), vnames, cat_vec,
                          weights(s), strata(s), clusters(s), des_ind, bin_size);
   //-------------------------------------------------------------------------------------------------------
   
   
   //------- the Right node -----------------------------------------------------------------------------
   
   
   //  if(right_cats.n_elem<1)cout<<"########### right_cats is "<<right_cats<< " "<<endl; 
   arma::uword ncats_R = rcats.n_elem;
   vec mxval_R = uvals(rcats);
   
   modfit = survLm_model(y.elem(not_s), mX.rows(not_s), weights.elem(not_s), strata.elem(not_s), clusters.elem(not_s));
   List SpR=get_node(2*node+1, cat=1, var, y, weights, mxval_R, not_s, modfit);
   List SpR2=random_split2(SpR["node"], y(not_s), mX.rows(not_s), X.rows(not_s), vnames, cat_vec,
                          weights(not_s), strata(not_s), clusters(not_s), des_ind, bin_size);
   
   //---------------------------------------------------------------------------------------------
   
   return rbind_splits(rbind_splits(SpL, SpL2), rbind_splits(SpR, SpR2)); 
   
 } //end categorical
 
 return(null_split()); //if neither numerical or categorical 
 
 
  
} //end random_split

//################################# End random_split ##############################################################



//         OLD ONE

//####################################################################################################################################
//
//                                                        random_split 
//
//  produces random variable and split point ---- for use with forest 
//####################################################################################################################################
// [[Rcpp::export]]
List random_split(arma::uword node, arma::colvec y, arma::mat mX, arma::mat X, Rcpp::StringVector vnames, arma::uvec cat_vec, 
                  arma::colvec weights, arma::uvec strata, arma::uvec clusters, arma::uvec des_ind,  
                  arma::uword &bin_size, arma::uvec randv){
  
  //C++ check for user interupt
  R_CheckUserInterrupt();
  
  arma::uword n = y.n_elem;
  
  
  // probabilit of splitting each time is between minp and 100 percent
  if(node>1){
    double minp=.98;
    double psplit=minp + (1-minp)*randu(), Runif=randu();
    if(psplit <Runif) return(null_split());
  }
  
  
  // Random chance to stop splitting
  // chance of spliting at each node starting at root
  // 0.98 0.92 0.87 0.82 0.77 0.73 0.68 if .98 multiply by 3
  // 0.99 0.94 0.89 0.85 0.81 0.77 0.73 if .99 muliply by 5
  
  //  double alpha = (int)floor((float)log2(node)), p_split = pow(.99, 5*alpha+1), Runif=randu();
  //if(p_split<Runif) return(null_split());
  
  // check if reached bin-size limit
  
  if(n<=2*bin_size || n < mX.n_cols) return(null_split());
  
  //cout<<"randv is "<<randv<<endl;
  
  //Get a count of the number of unique values each variable has
  uword pV=randv.n_elem;
  
  uvec uvar_count(pV);
  for(uword i=0; i<pV; ++i){
    vec unique_vals=unique(X.col(randv.at(i)));
    uvar_count(i)=unique_vals.n_elem;
  }
  
  
  //Randomly select variable that has more than one value
  uvec Vindx=find(uvar_count>1);  //Vindx is an index of variables with more than one value
  
  if(Vindx.n_elem==0){return(null_split());}
  uvec Rindx=shuffle(Vindx);
  uword x = randv.at(Rindx(0)); // variable selected randomly 
  
  arma::vec x_vals=X.col(x);
  arma::vec svals= sort(x_vals); // sorted unique x values
  arma::vec uvals= unique(x_vals); // sorted unique x values
  arma::uword cat=cat_vec.at(x);
  std::string var=as<string>(vnames.at(x));
  
  
  //*********************************** Get split point and loss **************************************** 
  
  //#------------------------ numeric variable  -------------------------------
  if(cat==0) {
    
    // ------------------ find random split value of x between edges 
    // arma::vec svals=sort(X.col(x));
    uvec s, not_s, rand_s;
    
    s=find(svals>svals(bin_size-1) && svals<svals(n-bin_size)); //index of xvalues that are split-able
    if(s.n_elem < 1 ) return(null_split());
    
    arma::rowvec xval(1);  //1x1 matrix for split value
    rand_s=shuffle(s); // random index off svals
    
    // cout<<" first value " << svals(rand_s(0))<< " second value " << svals(rand_s(0)+1) << endl;
    
    uvec indx = find(uvals==svals(rand_s(0)));
    
    xval(0)=(uvals(indx(0))+uvals(indx(0)+1))/2;
    //xval(0)=svals(rand_s(0)); //value of x rep as vector of size one
    
    //--------- Find loss ----------------------------------
    s=find(x_vals<=xval(0));
    not_s = find(X.col(x)>xval(0));
    if(s.n_elem < 1 || not_s.n_elem < 1){ return(null_split()); }  
    
    arma::colvec res_s = survLm_fit(y.elem(s), mX.rows(s), weights.elem(s))["residuals"];
    arma::colvec res_ns = survLm_fit(y.elem(not_s), mX.rows(not_s), weights.elem(not_s))["residuals"];
    double loss=get_loss(xval(0), x_vals, y, mX, weights);
    if(loss == datum::inf){ return(null_split()); }  //end if no split
    
    
    //================================================= Get nodes =======================================================
    //------- the Left node ---------------------------
    List modfit = survLm_model(y.elem(s), mX.rows(s), weights.elem(s), strata.elem(s), clusters.elem(s));
    
    List SpL=get_node(2*node, cat=0, var, y, weights, xval, s, modfit);
    List SpL2=random_split(SpL["node"], y(s), mX.rows(s), X.rows(s), vnames, cat_vec,
                           weights(s), strata(s), clusters(s), des_ind, bin_size, randv);
    //-------------------------------------------------------------------------------------------------------
    
    //------- the Right node -----------------------------------------------------------------------------
    // sR is not_s
    modfit = survLm_model(y.elem(not_s), mX.rows(not_s), weights.elem(not_s), strata.elem(not_s), clusters.elem(not_s));
    
    List SpR=get_node(2*node+1, cat=0, var, y, weights, xval, not_s, modfit);
    List SpR2=random_split(SpR["node"], y(not_s), mX.rows(not_s), X.rows(not_s), vnames, cat_vec,
                           weights(not_s), strata(not_s), clusters(not_s), des_ind, bin_size, randv);
    //---------------------------------------------------------------------------------------------
    
    return rbind_splits(rbind_splits(SpL, SpL2), rbind_splits(SpR, SpR2)); 
    
  } // end numeric
  //================================================= got nodes =======================================================
  
  //#------------------------ categorical variable  -------------------------------  
  if(cat==1){
    
    //uvar is "x.cats" unique set of categories
    arma::uword p = uvals.n_elem;  
    
    if(p<2) return(null_split());  //end if no split
    
    // get vector of uvecs for each category
    vector<uvec> cat_ind(p); //vector to hold indexs for each category
    arma::uvec cat_size(p);
    
    for(uword c=0; c<p; ++c){
      uvec ci = find(X.col(x)==uvals(c)); //find index of each category label 
      cat_ind.at(c)= ci;
      cat_size(c)=ci.n_elem;
    }
    // got the vector of indices for each category
    
    imat power_set;
    
    // ----- check if doing means then do linear order of categories ------
    if(mX.n_cols==1){
      
      power_set=imat(p-1, p);
      power_set.zeros();
      
      vec cat_means(p);
      cat_means.fill(datum::inf);
      
      for(uword i=0; i<p; ++i){
        //  uvec cs = find(X.col(x)==uvar(i));
        uvec cs = cat_ind.at(i);
        // cat_means(i) = mean(y(cs));  // changed this to be the survey weighted average
        cat_means(i) = sum(y(cs)%weights(cs))/sum(weights(cs));
      } 
      
      for(uword i=0; i<(p-1); ++i){
        arma::uword cmin = cat_means.index_min();
        power_set.rows(i, p-2).col(cmin).fill(1); //all rows below 
        cat_means(cmin) = datum::inf;  //set min to inf then go find next smallest mean
      }
      
    } //end if doing means
    
    else{  //not doing means
      // if(p>12) warning("Categorical variable ", x, " has > than 12 categories.  This may take a long time when fitting models on each node.");
      power_set = get_set_grid(p);
    } 
    //cout<<power_set;
    
    
    uvec s, not_s, rand_s;
    
    // -------- Find  rows of power_set that are useable------------------------------
    uvec row_size(power_set.n_rows);
    for(uword i=0; i<power_set.n_rows; ++i){
      s=find(power_set.row(i)==1);
      row_size(i)=sum(cat_size(s));
    }
    
    //--- set power_set to only include usable rows ----
    s=find(row_size>bin_size && row_size < n-bin_size);
    if(s.n_elem==0) return(null_split()); 
    // ------------- got powwer set --------------------------------------------------
    
    // ------------------ get split value = xval and loss at xval = loss ---------------------------------------
    
    // get random split
    rand_s=shuffle(s);
    uword xval = rand_s(0); 
    
    s.reset();
    uvec lcats=find(power_set.row(xval)==1);
    for(uword c=0; c<lcats.n_elem; ++c){ //loop to get s
      s= join_cols(s, find(X.col(x)==uvals(lcats(c))));
    }
    s=vectorise(s);
    
    uvec rcats=find(power_set.row(xval)==0);
    for(uword c=0; c<rcats.n_elem; ++c){ //loop to get not_s
      not_s= join_cols(not_s, find(X.col(x)==uvals(rcats(c))));
    }
    not_s=vectorise(not_s);
    
    if(s.n_elem < 1 || not_s.n_elem < 1){ return(null_split()); }  
    
    // get residuals
    arma::colvec res_s = survLm_fit(y.elem(s), mX.rows(s), weights.elem(s))["residuals"];
    arma::colvec res_ns = survLm_fit(y.elem(not_s), mX.rows(not_s), weights.elem(not_s))["residuals"];
    //    double loss=dot(weights.elem(s)%res_s, res_s) + dot(weights.elem(not_s)%res_ns,res_ns);
    
    //---------------------- got xval and loss ------------------------------------------- 
    
    
    //------- the Left node ---------------------------
    
    // s= sL   uvec left_cats = find(power_set.row(xval)==1);
    arma::uword ncats_L = lcats.n_elem;
    vec mxval_L = uvals(lcats);
    
    //get index vector of all data with category left side
    
    List modfit = survLm_model(y.elem(s), mX.rows(s), weights.elem(s), strata.elem(s), clusters.elem(s));
    List SpL=get_node(2*node, cat=1, var, y, weights, mxval_L, s, modfit);
    List SpL2=random_split(SpL["node"], y(s), mX.rows(s), X.rows(s), vnames, cat_vec,
                           weights(s), strata(s), clusters(s), des_ind, bin_size, randv);
    //-------------------------------------------------------------------------------------------------------
    
    
    //------- the Right node -----------------------------------------------------------------------------
    
    
    //  if(right_cats.n_elem<1)cout<<"########### right_cats is "<<right_cats<< " "<<endl; 
    arma::uword ncats_R = rcats.n_elem;
    vec mxval_R = uvals(rcats);
    
    modfit = survLm_model(y.elem(not_s), mX.rows(not_s), weights.elem(not_s), strata.elem(not_s), clusters.elem(not_s));
    List SpR=get_node(2*node+1, cat=1, var, y, weights, mxval_R, not_s, modfit);
    List SpR2=random_split(SpR["node"], y(not_s), mX.rows(not_s), X.rows(not_s), vnames, cat_vec,
                           weights(not_s), strata(not_s), clusters(not_s), des_ind, bin_size, randv);
    
    //---------------------------------------------------------------------------------------------
    
    return rbind_splits(rbind_splits(SpL, SpL2), rbind_splits(SpR, SpR2)); 
    
  } //end categorical
  
  return(null_split()); //if neither numerical or categorical 
  
} //end rpms_split

//################################# End random_split ##############################################################




//                             End of File
//###############################################################################################################
