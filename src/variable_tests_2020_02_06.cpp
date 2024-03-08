// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include<Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

List survLm_fit(arma::colvec y, arma::mat X, arma::colvec weights); //declaration: defined in model_functions

//int select_var(arma::colvec y, arma::mat mX, arma::mat vars, arma::uvec cat_vec, arma::colvec weights, arma::uvec clusters, arma::uvec des_ind,
//               Rcpp::StringVector vnames, arma::uword perm_reps, float pval);


List get_clus_ids(arma::uvec clus){
  
  arma::uvec clus_id= unique(clus); //ids of clusters
  uword C=clus_id.n_elem; //C is number of clusters
  //if(C<1) cout<<"############ C is "<<C<<endl;
  vector<uvec> cind(C);
  
  //find cluster effect for each cluster 
  for(uword c=0; c<C; ++c){
    //find elements in cluster c 
    uvec ci = find(clus==clus_id(c));
    //if(ci.n_elem<1)cout<<"########### ci has "<<ci.n_elem<< " elements"<<endl;
    cind.at(c) = ci;
    //  clus_size(c)=ci.n_elem;  //uncomment if cluster sizes wanted later
  }
  
  return List::create(
    Named("C")=C,
    Named("clusindx") = cind
  );
  
}//end get_clus_ids

arma::vec get_clus_effect(arma::vec rij, uword C, vector<uvec> cind){
  
//  arma::ivec clus_size(C); //vector to hold cluster sizes
  arma::vec clus_eff(C); // vector to hold cluster effects
 // if(C < 1) cout<<"C is "<<C<<endl;
  //find cluster effect for each cluster 
  for(uword c=0; c<C; ++c){
    clus_eff(c)= mean(rij(cind.at(c)));
    //  clus_size(c)=ci.n_elem;  //uncomment if cluster sizes wanted later
  }
  // cout<<"----------- The cind vector of vectors ----------" <<endl<<endl;
  // cout<<cind.at(0) <<endl << cind.at(1) <<endl << cind.at(2) <<endl;
  
  return clus_eff;
}

//#################### clus_perm ######################################
// Returns an n x M matrix containing M permutations of the weighted residuals 
//
// res: n-vec of weighted residuals; 
// clus_ids: C-vec of unique cluster ids
// effs: C-vec of estimated cluster effects;  clus: n-vec of cluster lables
// M: integer of number of permuted values to return
// 
//
//============================================================================
// [[Rcpp::export]] 
arma::mat clus_perm(arma::vec res, arma::uword C, std::vector <arma::uvec> clus_indx, arma::vec effs, arma::uword M) {
  
  int  n = res.n_elem;
  arma::mat r_vals(n, M); //matrix to hold M reordered weighted residuals
  
  for(uword m=0; m<M; ++m){
    arma::vec reffs=shuffle(effs);   // shuffled cluster effects
    
    arma::vec nr=res;
    
    for(uword i=0; i<C; ++i){
      
      // uvec j=find(clus==clus_ids(i)); //index of each observation in cluster i
      uvec j = clus_indx.at(i);
      
  //    if(j.n_elem<1) cout<<"j is "<<j<<endl;
      nr(j) -= effs(i); //subtract average old effect
      nr(j) += reffs(i); //add average reshuffled effect
      
      uvec sj = shuffle(j); //shuffle labels j within cluster
      
      nr(j)=nr(sj); // new value is shuffled values
      
    } //end mult-cluster loop
    
    r_vals.col(m)=nr;  //mth column gets permuted value of weighted residual
    
  } //end M perms loop
  
  return(r_vals);
  
} // end clus_perm

// --------------------------- end clus_perm ---------------------------------

//********************************* unclustered perm ***************************
// if design does not have clusters use faster unclustered permuation
// ****************************************************************************
//
// [[Rcpp::export]] 
arma::mat perm(arma::vec res, arma::uword M) {
  
  int  n = res.n_elem;
  
  //considering weighted residuals rij
  
  arma::mat r_vals(n, M+1); //matrix to hold M reordered weighted residuals
  
  r_vals.col(0)=res; // first column is original order of weighted residuals
  
  // begin permutaions
  
  //perm loop
  for(uword m=1; m<=M; ++m){ // m starts at 1 because first element is original
    
    r_vals.col(m)=shuffle(res);  //mth column gets permuted value of weighted residual
    
  } //end M perms loop
  
  return(r_vals);
  
}
//------------------------ end perm -----------------------------------------


//********************* peak functions for num and cat *************************
//Returns the highest value of the cumulative score function ordered
//by the numeric variable to find model variable instability


double peak_num(arma::colvec score, arma::vec var){
  
  uvec sindx = sort_index(var);
  
  double peak=max(abs(cumsum(score(sindx))));
  
  return peak;
}

// -------------------- end peak_num ---------------------------------------

// ################ peak for categorical variables  ######################
// Returns the highest value of the cumulative score function ordered
// by the categories to find model variable instability
// #####################################################################3


double peak_cat(arma::colvec score, arma::vec var){
  
  arma::vec cat_id= unique(var); // get the unique categories
  if(cat_id.n_elem <= 1) return (double)0;  //if only one category then can't split
  
  uword C=cat_id.n_elem;
  
  arma::vec cat_sums(C); // make vector to hold category sums
  
  for(uword i=0; i<C; ++i){
    arma::uvec s = find(var == cat_id(i));
    cat_sums(i)=sum(score(s));
  }
  
  uvec pos = find(cat_sums >= 0);
 // if(pos.n_elem<1) cout<<"pos is "<< pos << endl;
  uvec neg = find(cat_sums < 0);
//  if(neg.n_elem<1) cout<<"neg is "<<neg << endl;
  double peak, a, b;
  
  a=fabs((float)sum(cat_sums(neg)));
  b=sum(cat_sums(pos));
  
  if(a > b) peak=a; else peak=b;
  
  return peak;
}
// -------------------- end peak_cat ---------------------------------------
//*********************** end peak functions ********************************




//###############################################################################################################################
//          New get pvec  returns 2 p-vectors: p-valuse and cooresponding max_peak for each variable 
//
//          replaces  old  get_pvec and eliminates var_test
//################################################################################################################################

List get_pvec(arma::colvec scores, arma::mat mX, arma::mat vars, arma::uvec cat_vec,  
              arma::uvec clusters, arma::uvec des_ind, arma::uword perm_reps, float pval){
  

  arma::uword pX = mX.n_cols, pV = vars.n_cols, n=scores.n_elem, M=100;

  arma::uword alpha=(uword)ceil(pval*perm_reps); //find number of failures until dismiss
  perm_reps=(uword)ceil(perm_reps/M); //find number of reps of M to do

  arma::mat pvals(pV, pX), peaks(pV, pX); // matrices that store p-vals and peaks for each test
  pvals.fill(0);
  peaks.fill(0);

  arma::vec pvec(pV), maxpeak(pV); //these are two vectors that are returned
  
  List  clus_id = get_clus_ids(clusters);
  
  for(arma::uword j=0; j<pX; ++j){  //do for each variable in the model matrix
    
    arma::vec clus_eff;
    
    arma::colvec chi = scores;
    chi %= mX.col(j);         //multiply res'*X_j

//cout<<chi(0)<<" and "<<chi(1)<<endl;     
    double a;  //original peak
    
    if(des_ind[2] == 1){  // if clustered (vectors begin at 0)
   
      clus_eff = get_clus_effect(chi, clus_id["C"], clus_id["clusindx"]); // get clus effects
    } // end if clusterd
    
//############### Permutation Loop over perm_reps X 100 #############
    //------ get permuted values of residuals ---------------------    
   // for(arma::uword pr_indx=0; pr_indx<perm_reps; ++pr_indx){

   //---- compare each variable to obsereved peak for each permutation 

   for(uword i=0; i<pV; ++i){
    arma::uword pr_indx=0;
    while(pr_indx<perm_reps and pvals(i,j)<alpha){
      arma::mat pres(n, M); // to store M permuted values of the original residuals
      
      if(des_ind[2] == 1){ // if clustered (vectors begin at 0)
        pres = clus_perm(chi,  clus_id["C"], clus_id["clusindx"], clus_eff, M);
      }
      else{  // if not clustered 
        pres = perm(chi, M);
      } // end cluster else 
      //----------------------------------------------------
 
        // ------- categorical variable --------------
        if(cat_vec[i] == 1){
          
          a = peak_cat(chi, vars.col(i));
          peaks(i, j) = a;
          double peak;
          arma::uword m=0;
          while(m <M and  pvals(i,j)<alpha){
         // for(arma::uword m=0; m<M; ++m){
            
            peak = peak_cat(pres.col(m), vars.col(i));
            
            if(peak > a) pvals(i, j) += 1; //higher than observed peak
         //maybe stupid   if(peak > peaks(i,j)) peaks(i, j) = peak; //find max peak
            ++m;
          } //end m loop
          ++pr_indx;
        }//end if cat
        else{
          //--------- numeric variable -----------------------
          a = peak_num(chi, vars.col(i));
      //    cout << "a is "<<a<<endl;
          peaks(i, j) = a;
          arma::uword m=0;
          double peak;
          while(m <M and  pvals(i,j)<alpha){
          //for(arma::uword m=0; m<M; ++m){
            
            peak = peak_num(pres.col(m), vars.col(i));
          //  cout<<"peak is "<<peak<<endl;
            if(peak >= a){ pvals(i, j) += 1;}  //higher or equal to peak on observed data
          // might be stupid  if(peak > peaks(i,j)){ peaks(i, j) = peak; }  //find max peak
            ++m;
          } //end m loop
          ++pr_indx;
        } //end else numerical
      } // end while  loop over perm_reps 
    } //end i = 0..pV loop over variables
  } // end j loop over mX variables
  pvals = (.5 + pvals)/((double)(perm_reps*M+1));
  
  // for each i return lowest pval and corresponding peak
  // maybe return full vector and add this code-snippet to split to collect importance measure

  for(arma::uword i=0; i<pV; ++i){
    double minpval = min(pvals.row(i));
    arma::uvec s = find(pvals.row(i)==minpval);
   // if(s.n_elem<1)cout<<"########### ci has "<<s.n_elem<< " elements"<<endl;
    arma::rowvec potpeaks = peaks.row(i); // peaks for ith variable
    double mpeak = max(potpeaks(s)); // maximum of all potential peaks
    pvec(i) = minpval;
    maxpeak(i) = mpeak;
  }

  return List::create(
    Named("pvec") = pvec,
    Named("peaks") = maxpeak
  );
  
} // end get_pvec function
//#########################################################################################################################################################


//###########3#################################### select_var funciton #############################################
//
//   returns variable index or -1 if no variable is found to be significant
//
//#####################################################################################################################
int select_var(arma::colvec y, arma::mat mX, arma::mat vars, arma::uvec cat_vec, arma::colvec weights, arma::uvec clusters, arma::uvec des_ind,
               Rcpp::StringVector vnames, arma::uword perm_reps, float pval){
  
  //cout<<"in select var "<<endl;
  //--- check for variables with only 1 unique value ----------
  uword pV=vars.n_cols;
  uvec uvar_count(pV);
  for(uword i=0; i<pV; ++i){
    vec unique_vals=unique(vars.col(i));
    uvar_count(i)=unique_vals.n_elem;
  }
  
  uvec Vindx=find(uvar_count>1);  //Vindx is an index of variables with more than one value
  if(Vindx.n_elem==0){ return(-1); }
 //don't like this sloppy reusing variable names
 // vars = vars.cols(Vindx);  //only consider variables with at least two unique values
  //------------- removed vars with only 1 unique value------
  

  //----- get residuals ---- 
  arma::colvec res = survLm_fit(y, mX, weights)["residuals"];
  res%=weights; //considering weighted residuals rij
  
  
  //cout<<"going to get_pvec"<<endl;
  List pvec = get_pvec(res, mX, vars.cols(Vindx), cat_vec,  clusters, des_ind, perm_reps, pval);
  arma::vec pvals = pvec["pvec"];
  if(pvals.min() <= pval){
 //   cout<<"pvals is "<<pvals<<endl;
   //uvec s = find(pvals == pvals.min());
   uvec s = find(pvals <= pval);
    arma::vec peaks = pvec["peaks"];

    
    if(s.n_elem == 1) return(Vindx(s(0)));
    else{
      uvec a = shuffle(find(peaks.elem(s)==max(peaks(s)))); //
      return(Vindx(s(a(0)))); // return random max if more than 1 minimum
    }//end else more than 1 min
    
  }//end if there are significant variables
  
  return(-1); //  -1 no significant variables
}//end select_var
//####################################### end var_select #################################################################
//############################################################ End var select file ##########################################################




