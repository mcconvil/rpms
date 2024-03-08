-----------   Version 0_4_1  ------------------

Allows rpms_forest() an option to run parallel on non-Windows systems

------------  Version 0_4_0  ------------------

Provides function rpms_forest() which fits a random forest model to survey data

Uses loss that is not weighted by proportion of sample
  - leads to more interpretable splits
  
Improved version of qtree() function  
  - allows option to show sample size n

-------- Version 0_3_0 --------------

  The whole recusive partitioning algorithm is now done in C++
  
    -fixes rounding errors resulting from moving between R and C++
  
  Missing values for categorical variables used for recursive partitioning are treated as their own category 
  
  More options added to the qtree() function
  
  Provides a full year of Consumer Expenditure (CE) dataset with more interesting variables



--------- Version 0_2_1 ---------------

  Fixes a number of bugs.
  
  Makes splitting on categories faster when e_equ=y~1
  
  Updates variable selection procedure with new perm reps
  
  New version of qtree with more options



---------- Version 0_2_0 ----------
  
  Provide a new function end_nodes()
  
  Fixed bugs in examples using end_nodes()
  
  Included vignette "rpms"
  
  Changed variable selection procedure in case of ties in p-value evaluation
  
  Added labels and caption options to qtree() function
  

---------- Version 0_1_0 -----------

  provides main rpms function with 
    
    -two methods: print() and predict()

  other functions: in_node(), qtree(), node_plot()

   