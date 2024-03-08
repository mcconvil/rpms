// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include<Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;


//############################## BREAK UP INTO SEPARATE FILE  var select ######################################################################

// [[Rcpp::export]] 
List survLm_fit(arma::colvec y, arma::mat X, arma::colvec weights){

  //weights = sqrt(weights);
  arma::mat wX = X;
  wX.each_col() %= sqrt(weights);  // multiply each column by sqrt(weights)
  
  arma::mat iA = arma::trans(wX)*wX;

  // check for colinearitiy in X before inverting
  if(iA.is_sympd()){ 
    
    iA = inv_sympd(iA);
    
    vec wy=y;
    wy %=sqrt(weights);
    
    arma::colvec coef = iA*(arma::trans(wX)*wy);
    
    // residuals-- these are unweighted
    arma::colvec resid = y - X*coef;   
    
    return List::create(Rcpp::Named("coefficients") = coef,
                        Rcpp::Named("residuals")    = resid);
    
  } // end positve definite if singular
  else if(arma::det(iA) != 0){
    
    iA = inv(iA);
    
    vec wy=y;
    wy %=sqrt(weights);
    
    arma::colvec coef = iA*(arma::trans(wX)*wy);

    // residuals-- these are unweighted
    arma::colvec resid = y - X*coef;   
    
    return List::create(Rcpp::Named("coefficients") = coef,
                        Rcpp::Named("residuals")    = resid);
  
  } // end not singular
  
  else{
    
    arma::colvec coef(X.n_cols, fill::zeros);
    arma::colvec resid(y.n_elem); resid.fill(datum::inf);
    
    return List::create(Rcpp::Named("coefficients") = coef,
                        Rcpp::Named("residuals")    = resid); 
  } // end singular
  
  
  
} //end survLm_fit


//**************** getting design based estimates of std.error for coefficient estimates ***
arma::mat survLM_covM(arma::colvec resid, arma::mat X, 
                      arma::colvec weights, arma::uvec strata, arma::uvec clusters) {
  
  size_t n = X.n_rows, p = X.n_cols; 
  
  arma::colvec sweights = sqrt(weights);
  X.each_col() %= sweights;  // multiply each column by sqrt(weights)
  
  arma::mat iA = arma::trans(X)*X;
  
  if(arma::det(iA)!=0){
    if(iA.is_sympd()) iA = inv_sympd(iA); else  iA = inv(iA);
    
    arma::mat D = X;        //has weight sqrt(Wt) since X has sqrt(Wt)
    D.each_col() %= resid;  //is sqrt(Wt)*X*resid
    
    arma::mat G(p, p, fill::zeros);
    
    arma::uvec h_labs=unique(strata);
    
    uword H = h_labs.n_elem;
    
    for(uword h=0; h<H; ++h){
      arma::mat G_h(p, p, fill::zeros);
      uvec uh = find(strata==h_labs[h]);
      
      //arma::uvec ids = unique(clusters.elem(find(strata==h_labs[h])));  //why am i finding this again?
      arma::uvec ids = unique(clusters.elem(uh));
      uword nh = ids.n_elem;
      
      for(uword i=0; i<nh; ++i ){
        uvec s=find(clusters.elem(uh) == ids[i]);  //in strata h and cluster i 
        arma::rowvec dhi = sum(D.rows(s));
        
        G_h += (trans(dhi)*dhi);
        
      } //end in cluster loop
      
      if(nh <= 1) {G+=datum::inf;}  // single psu
      else G += (nh/(nh-1))*G_h; // 
      
    }
    
    if(n <= p) {G=datum::inf;}  // single psu
    else G=((n-1)/(n-p))*G; // 
    G=((n-1)/(n-p))*G;
    
    arma::mat covM = iA*G*iA;
    return covM;
  }
  else{  // matrix is singular
    
    arma::mat covM(p, p);
    covM.fill(datum::inf);
    return covM;
  
  } // end det==0

  
} //End get_CovM



//' Fit a linear model using data collected from a complex sample
//' 
//' @param y A vector of values
//' @param X The design matrix of the linear model  
//' @param weights A vector of sample weights for each observation 
//' @param strata A vector of strata labels
//' @param clusters A vector of cluster labels
//' 
//' @return list containing coefficients, covariance matrix and the residuals
//' 
//' @keywords internal
//' 
//' @export
// [[Rcpp::export]] 
List survLm_model(arma::colvec y, arma::mat X, arma::colvec weights, 
                  arma::uvec strata, arma::uvec clusters) {

 // cout<<"in survLm_model"<<endl;
  List Lnmod = survLm_fit(y, X, weights);
//  cout<<"did the fit"<<endl;
  arma::mat covM = survLM_covM(Lnmod["residuals"], X,
                              weights, strata, clusters);
  
  return List::create(Rcpp::Named("coefficients") = Lnmod["coefficients"],
                      Rcpp::Named("covM")         = covM,
                      Rcpp::Named("residuals")    = Lnmod["residuals"]);
  
}

//########## end model functions ######################################################################################
