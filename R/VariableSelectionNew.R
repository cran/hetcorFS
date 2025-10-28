
#######################################################
###                                                 ###
###               Selection of Variables            ###
###                                                 ###
#######################################################

#' Selection of Variables
#'
#' Filtering feature selection based on the heterogeneous correlation matrix.
#'
#' @param data A data frame. Values of type 'numeric' or 'integer' are treated as numerical.
#' @param alpha significance level to be used for testing, default = 0.05.
#' @param missing pairwise complete by default, set to TRUE for complete deletion
#' @param pv_adj correction method for p-value "none" by default, for options see p.adjust
#' @param smooth.tol the minimum acceptable eigenvalue for the smoothing default 10^-12
#' @param method algorithm used, c (cell-wise) by default,the other option is r (row-wise)
#'
#' @return An list of elements:
#'    \item{rearranged.data.set}{A data frame with the original variables with numerical features first }
#'    \item{selected.variables}{A data frame of the selected variables}
#'    \item{variable.indices}{The indices of the selected variables from the original data frame }
#'    \item{original.corr.matrix}{The \eqn{p} by \eqn{p} extended correlation matrix of all the inputted variables}
#'    \item{corr.matrix}{The \eqn{d} by \eqn{d} extended correlation matrix of the selected variables}
#'    \item{original.p.value.matrix}{The \eqn{p} by \eqn{p} p-value matrix of all the inputted variables}
#'    \item{p.value.matrix}{The \eqn{d} by \eqn{d} p-value matrix of the selected variables}
#'
#' @examples
#' data(Employee_Satisfaction_Index)
#' feature_selection(Employee_Satisfaction_Index)


UFS <- function(data=NULL,alpha=0.05,missing=FALSE,pv_adj='none',smooth.tol=10^-12,method='c'){
  if(missing){
    rows_with_na <- which(apply(data, 1, function(row) any(is.na(row))))
    data=data[-rows_with_na,]
  }
  if(method=='c'){method='TPM'}else if(method=='r'){method='GHSVK'}else{
    warning('the method specified does not exsist, using default cell-wise c')
    method='TPM' }
  the_data=data
  the_data_unchanged <- the_data
  data_names <- colnames(the_data)
  m_and_p<- HCPM(the_data)
  corr_mat=m_and_p$cor_mat
  p_mat=m_and_p$p_value
  #rownames(corr_mat)=data_names
  #colnames(corr_mat)=data_names
  #rownames(p_mat)=data_names
  #colnames(p_mat)=data_names
  #### smoothing the corr_mat
  corr_mat=psych::cor.smooth(corr_mat,eig.tol=smooth.tol)
  orig_corr_mat <- corr_mat #extended correlation matrix, left as is
  if(method=='TPM'){
  ##### Adjusting the p-values
  if(pv_adj!="none"){
    p_mat=matrix(p.adjust(p_mat,method = pv_adj),nrow(corr_mat),nrow(corr_mat))
  }

  orig_p_mat <-  p_mat#p-value matrix, left as is


  main=core_var_sel(the_data=the_data,corr_mat=corr_mat,p_mat=p_mat,alpha=alpha)
  }else{
    orig_p_mat <-  p_mat#p-value matrix, left as is
    main=core_var_selH(the_data=the_data,corr_mat=corr_mat,p_mat=p_mat,alpha=alpha,pv_adj=pv_adj)

  }

  output_list <- list('rearranged.data.set' = the_data_unchanged,'selected.features' = main$selected.variables,'feature.indices' = main$variable.indices,'original.corr.matrix'=orig_corr_mat,'corr.matrix' = main$corr.matrix,'original.p.value.matrix' = orig_p_mat,'p.value.matrix' =main$p.value.matrix) #creates a list of outputs
   return(output_list) #return the outputs
}




core_var_sel<-function(the_data,corr_mat,p_mat,alpha){
  data_names=colnames(the_data)
  or_corr_mat=corr_mat #store original correlation matrix
  p = nrow(corr_mat) #find the number of rows of the correlation matrix (same as number of columns), should match p-value matrix
  diag(p_mat) = 0.1 + alpha #make diagonals a function of alpha so function does not drop all columns
  while (any(p_mat<=alpha) & p>2){ #make a condition so the process will repeat
    corr_vector <- c() #create an empty vector
    for (i in 1:(p-1)){
      for (j in (i+1):p){ #loop over the number of rows and columns (lower-triangular)
        if (p_mat[i,j] <= alpha){ #check if the specific p-value is <=0.05
          corr_vector=c(corr_vector,abs(corr_mat[i,j])) #if so, append absolute value of respective correlation coefficient to empty vector
        }
      }
    }
    corr_vector <- sort(corr_vector,decreasing = T) #sort the correlation vector in decreasing order
    indices <- as.data.frame(which(abs(corr_mat)==corr_vector[1],arr.ind = T)) #find the indices of the highest correlation coefficient from the correlation matrix
    sum_vector <- c() #create an empty vector
    for (k in 1:length(indices$row)){ #loop over the length of the indices dataframe
      sum_vector=c(sum_vector,sum(abs(corr_mat[indices$row[k],]))) #sum the absolute value of correlation coefficients in the row of the stores indices, append the sum
      #to empty vector
    }
    indices <- cbind(indices,sum_vector) #column bind the row sums vector to the index data frame
    indices_arranged <- arrange(indices,desc(sum_vector)) #arrange the data frame based on a descending order for the sum vector
    desired_index <- indices_arranged$row[1] #the desired index will be stored as the row corresponding to the highest row sum
    corr_mat <- corr_mat[-c(desired_index),] #remove that indexed row from the correlation matrix
    corr_mat <- corr_mat[,-c(desired_index)] #remove that indexed column from the correlation matrix
    p_mat <- p_mat[-c(desired_index),] #remove that indexed row from the p-value matrix
    p_mat <- p_mat[,-c(desired_index)] #remove that indexed column from the p-value matrix
    p = nrow(corr_mat) #save the number of rows (same as columns) of the correlation matrix (same as p-value matrix)
  }
  matrix_names <- colnames(corr_mat) #save the column names of the reduced correlation matrix as a vector
  final_indices <- NULL #create an empty vector
  for (m in 1:length(matrix_names)){ #loop over the length of the matrix names vector
    final_indices = append(final_indices,which(data_names==matrix_names[m])) #append the indices of matrix names from the original dataframe
  }
# var_vector <- c()
 # Mylevels=c("character","factor","integer","numeric") #creates a vector of class types
  if (any(p_mat<=alpha) & p==2){ #check to see if two variables are left
    sumcor=apply(abs(or_corr_mat),2,sum)##row sum of original correlation
    selected = which.max(sumcor[final_indices]) #selects the index of the maximum correlation
    final_indices = which(data_names==names(selected)) #finds original index

    }
  updated_data <- the_data[,final_indices] #subset the original dataframe using the newly-obtained indices

  diag(p_mat) = 0 #make diagonals of p-value matrix 0 again
  output_list <- list('selected.variables' = updated_data,'variable.indices' = final_indices,'corr.matrix' = corr_mat,'p.value.matrix' = p_mat) #creates a list of outputs


}


core_var_selH<-function(the_data,corr_mat,p_mat,alpha,pv_adj){
  data_names=colnames(the_data)
  or_corr_mat=corr_mat #store original correlation matrix
  p = nrow(corr_mat) #find the number of rows of the correlation matrix (same as number of columns), should match p-value matrix
  diag(corr_mat)=0
  O_p_mat <-  p_mat#p-value matrix, left as is
  orig_p_mat <-  p_mat#p-value matrix to change

  meanCor=apply(abs(corr_mat),1,mean)
  meanCorO=meanCor
  corLevels=NULL
  diag(orig_p_mat) = 0.1 + alpha #make diagonals a function of alpha so function does not drop all columns


  index=which.max(meanCor)
   if(pv_adj!="none"){
      p_mat[index,-index]=p.adjust(orig_p_mat[index,-index],method = pv_adj)

     while (any(p_mat[index,]<=alpha)&p>2 ){ #make a condition so the process will repeat
        corLevels=c(corLevels,meanCor[index])
        #if(any(p_mat[index,]<=alpha)){
        desired_index <- index
        corr_mat <- corr_mat[-c(desired_index),] #remove that indexed row from the correlation matrix
        corr_mat <- corr_mat[,-c(desired_index)] #remove that indexed column from the correlation matrix
        meanCor=apply(abs(corr_mat),1,mean)
        index=which.max(meanCor)
        orig_p_mat <- orig_p_mat[-c(desired_index),] #remove that indexed row from the p-value matrix
        orig_p_mat<- orig_p_mat[,-c(desired_index)] #remove that indexed column from the p-value matrix
        p_mat=orig_p_mat
        p_mat[index,-index]=p.adjust(orig_p_mat[index,-index],method = pv_adj)
        p = nrow(corr_mat) #save the number of rows (same as columns) of the correlation matrix (same as p-value matrix)

    # }
      }
   }else{
     p_mat=orig_p_mat
     while (any(p_mat[index,]<=alpha)&p>2 ){ #make a condition so the process will repeat
       corLevels=c(corLevels,meanCor[index])
       #if(any(p_mat[index,]<=alpha)){
       desired_index <- index
       corr_mat <- corr_mat[-c(desired_index),] #remove that indexed row from the correlation matrix
       corr_mat <- corr_mat[,-c(desired_index)] #remove that indexed column from the correlation matrix
       meanCor=apply(abs(corr_mat),1,mean)
       index=which.max(meanCor)
       orig_p_mat <- orig_p_mat[-c(desired_index),] #remove that indexed row from the p-value matrix
       orig_p_mat<- orig_p_mat[,-c(desired_index)] #remove that indexed column from the p-value matrix
       p_mat <- orig_p_mat #remove that indexed row from the p-value matrix
        p = nrow(corr_mat) #save the number of rows (same as columns) of the correlation matrix (same as p-value matrix)

    # }
   }
   }
  matrix_names <- colnames(corr_mat) #save the column names of the reduced correlation matrix as a vector
  final_indices <- NULL #create an empty vector
  for (m in 1:length(matrix_names)){ #loop over the length of the matrix names vector
    final_indices = append(final_indices,which(data_names==matrix_names[m])) #append the indices of matrix names from the original dataframe
  }
  if (any(orig_p_mat<=alpha) & p==2){ #check to see if two variables are left
    selected = which.max(meanCorO[final_indices]) #selects the index of the maximum correlation
    final_indices = which(data_names==names(selected)) #finds original index
  }



  updated_data <- the_data[,final_indices] #subset the original dataframe using the newly-obtained indices

  # var_vector <- c()
  # Mylevels=c("character","factor","integer","numeric") #creates a vector of class types
  diag(p_mat) = 0 #make diagonals of p-value matrix 0 again
  output_list <- list('selected.variables' = updated_data,'variable.indices' = final_indices,'corr.matrix' = corr_mat,'p.value.matrix' = p_mat,average_corr=meanCorO, corLevels= corLevels) #creates a list of outputs


}



#######################################################
###                                                 ###
###           Extended Correlation Matrix           ###
###                                                 ###
#######################################################

#' Extended Correlation Matrix
#'
#' Extends the traditional correlation matrix (between numerical data) to also
#' include binary and ordinal categorical data.
#'
#' @param data A data frame. Values of type 'numeric' or 'integer' are treated as numerical.
#'
#' @return A list with with elements:
#'    \item{cor_mat}{An \eqn{p} by \eqn{p} matrix of correlation coefficients}
#'    \item{p_value}{An \eqn{p} by \eqn{p} matrix of p-values}
#'
#' @examples
#'
#' data(Employee_Satisfaction_Index)
#' HCPM(Employee_Satisfaction_Index)
#'
HCPM <- function(data=NULL){
 class=sapply(data.frame(data),class)
 num=data[,which(class=="integer"|class=="numeric")]
 if(any(class=="integer"|class=="numeric")){
 ord=data[,-which(class=="integer"|class=="numeric")]
 }else{ord=data}
  if(!is.null(num)){
    num=as.matrix(num)
    pcon=ncol(num)
    n=nrow(num)
  }else{
    pcon=0
  }
  if(!is.null(ord)){
    ord=as.matrix(ord)
    pord=ncol(ord)
    n=nrow(ord)
  #  ord=apply(ord, 2,as.factor)
  }else{pord=0}

  p=pord+pcon
  cor_mat=matrix(0,p,p)
  diag(cor_mat)=1
  p_value=matrix(0,p,p)
  result=0
  if(pord>1){
    for(i in 1:(pord-1)){
      for(j in (i+1):(pord)){
        result=try(polycor::polychor(ord[,i],ord[,j],std.err = T),silent=TRUE)
        test=class(result)
        if(test != "try-error"){
        cor_mat[(pcon+i),(pcon+j)]=result$rho
        p_value[(pcon+i),(pcon+j)]= 1-pchisq((result$rho^2/result$var),1)
        }else{p_value[(pcon+i),(pcon+j)]=0.99}
      }
    }
  }
  result=0
  if(pord!=0&pcon!=0){
    for(i in 1:pcon){
      for(j in 1:pord){
        result=polycor::polyserial(num[,i], ord[,j],std.err = T)
        cor_mat[i,(pcon+j)]=result$rho
        p_value[i,(pcon+j)]= 1-pchisq((result$rho^2/result$var),1)
      }
    }
  }
  if(pcon>1){
    cor_mat[1:pcon,1:pcon]=cor(num,use='pairwise.complete.obs')
    ptest <- cor.test_extended(num[,1:pcon]) #create pearson p value matrix

    out=polycor::hetcor(num[,1:pcon],use='pairwise.complete.obs')
    if(any(which(out$tests>0.05))){
      HC <- out$correlations       # Heterogeneous Correlation matrix
      SE <- out$std.errors         # Standard errors
      WT <- (HC/SE)^2              # Wald's Statistic
      Wp <- 1 - pchisq(WT, df = 1) # p-values from Wald's Test
      ptest[which(out$tests>0.05)]=Wp[which(out$tests>0.05)]
    }
    p_value[1:pcon,1:pcon]=ptest
  }

  cn_num=colnames(num)
  cn_ord=colnames(ord)
  if(is.null(cn_num)&pcon>0){cn_num="n1"}
  if(is.null(cn_ord)&pord>0){cn_ord="o1"}
  col_name=c(cn_num,cn_ord)

  cor_mat[lower.tri(cor_mat)] <- t(cor_mat)[lower.tri(cor_mat)] #since lower-triangle was focus, we must match upper to lower
  diag(cor_mat)=1 #diagonals will always be 0 for the correlation coefficients
  p_value[lower.tri(p_value)] <- t(p_value)[lower.tri(p_value)] #match upper to lower for p-value matrix as well
  diag(p_value)=0 #diagonals will always be 0 for the p-values
  colnames(cor_mat)=col_name
  rownames(cor_mat)=col_name
  colnames(p_value)=col_name
  rownames(p_value)=col_name

  return(list(cor_mat=cor_mat,p_value=p_value))
}


cor.test_extended <- function(num){
  num = as.data.frame(num) #ensure num is a dataframe, especially if it is 1 column
  p = ncol(num) #save the number of columns in the numerical data
  p.mat <- matrix(0,nrow=p,ncol=p) #create a p-value matrix of 0s sized based on p
  if (p==1){
    result <- cor.test(num[,p],num[,p]) #run the cor.test function on each pair of variables
    p_value <- result$p.value #pull the p-value from the result
    p.mat[p,p] <- p_value #stores p-value in its respective position in the created matrix
  }
  else if (p>1){
    for (i in 1:(p-1)){
      for (j in (i+1):p){ #loop over the number of rows and columns in the created matrices (lower-triangular)
        result <- cor.test(num[,i],num[,j]) #run the cor.test function on each pair of variables
        p_value <- result$p.value #pull the p-value from the result
        p.mat[i,j] <- p_value #stores p-value in its respective position in the created matrix
      }
    }
  }
  p.mat[lower.tri(p.mat)] <- t(p.mat)[lower.tri(p.mat)] #since lower-triangle was focus, we must match upper to lower
  diag(p.mat)=0 #diagonals will always be 0 for the p-values
  return(p.mat) #return the p-value matrix
}




#######################################################
###                                                 ###
###                   Tuning Alpha                  ###
###                                                 ###
#######################################################

#' Tuning Alpha
#'
#' Creates a barplot of all variables to see which variables are selected at
#' different values of alpha.
#'
#' @param data A data frame. Values of type 'numeric' or 'integer' are treated as numerical.
#' @param grid.alpha A vector of alpha values to be plotted, default = seq(0.01,0.99,by=0.01).
#' @param missing pairwise complete by default, set to TRUE for complete deletion
#' @param pv_adj correction method for p-value "none" by default, for options see p.adjust
#' @param smooth.tol the minimum acceptable eigenvalue for the smoothing default 10^-12
#' @param method algorithm used, c (cell-wise) by default,the other option is r (row-wise)
#'
#' @return Displays a barplot depicting which variables are selected at each value
#'          of alpha (multiplied by 100) and a list with elements:
#'    \item{survivors}{Vector depicting how many alphas a variable is selected under}
#'    \item{data_names}{Vector depicting the corresponding names of the variables}
#' @examples
#' data(Employee_Satisfaction_Index)
#' tuning_alpha(Employee_Satisfaction_Index)


FS_barplot <- function(data = NULL, grid.alpha = seq(0.01,0.99,by=0.01),missing=FALSE,pv_adj='none',smooth.tol=10^-12,method="c"){
  if(missing){
    rows_with_na <- which(apply(data, 1, function(row) any(is.na(row))))
    data=data[-rows_with_na,]
  }
  if(method=='c'){method='TPM'}else if(method=='r'){method='GHSVK'}else{
    warning('the method specified does not exsist, using default cell-wise c')
    method='TPM' }
  the_data=data
  the_data_unchanged <- the_data
  data_names <- colnames(the_data)
  m_and_p<- HCPM(the_data)
  corr_mat=m_and_p$cor_mat
  p_mat=m_and_p$p_value
  rownames(p_mat)=data_names
  colnames(p_mat)=data_names
  #### smoothing the corr_mat
  corr_mat=psych::cor.smooth(corr_mat,eig.tol=smooth.tol)
  orig_corr_mat <- corr_mat#HCPM(num = num,bin = bin,ord_cat = ord_cat) #extended correlation matrix, left as is

  if(method=='TPM'){
    ##### Adjusting the p-values
    if(pv_adj!="none"){
      p_mat=matrix(p.adjust(p_mat,method = pv_adj),nrow(corr_mat),nrow(corr_mat))
    }

    orig_p_mat <-  p_mat#p-value matrix, left as is


  orig_p_mat <-  p_mat#p_value_matrix(num = num,bin = bin,ord_cat = ord_cat) #p-value matrix, left as is



  p = ncol(the_data) #save the number of columns in the dataset
  survivors <- rep(0,p) #create an empty vector as long as the dataset
  for (i in grid.alpha){ #loop over the values in the grid.alpha vector
    indices<- core_var_sel(the_data=the_data,corr_mat=corr_mat,p_mat=p_mat,alpha=i)$variable.indices
     survivors[indices] <- survivors[indices] + 1 #increase the counter for the selected variables
  
  }
  }else{
    orig_p_mat <-  p_mat#p-value matrix, left as is
    p = ncol(the_data) #save the number of columns in the dataset
    survivors <- rep(0,p) #create an empty vector as long as the dataset
    for (i in grid.alpha){ #loop over the values in the grid.alpha vector
      indices<- core_var_selH(the_data=the_data,corr_mat=corr_mat,p_mat=p_mat,alpha=i,pv_adj = pv_adj)$variable.indices
      survivors[indices] <- survivors[indices] + 1 #increase the counter for the selected variables
    }

  }
  grid=c(0,grid.alpha)
  pp= barplot(height=grid[survivors+1]*100, names.arg=data_names,las=3,yaxt='n',ylab='Percent level of significance',ylim=c(0,max(grid[survivors+1]*100)+5))#expression(paste(italic(alpha) %.% 100," (%)",sep=""))
  # barplot(height=survivors/lenght(grid.alpha),names.arg=data_names, yaxt='n') #create a barplot of the counts
  #axis(side=2,at=grid*100,las=1) #create a y axis for the alpha values
  axis(side=2,las=1) #create a y axis for the alpha values
  text(pp,grid[survivors+1]*100+3,labels=round(rank(-survivors,ties.method = 'min')))

  # title(xlab = "Variables",ylab = expression(paste(italic(alpha) %.% 100," (%)",sep=""))) #labels
  return(list(survivors,data_names)) #returns the survivors vector along with plotting the barplot
}



RedRate<-function(data_red){
  p=ncol(as.matrix(data_red))
  if(p==1){RedRate=0
  }else{
    cor=HCPM(data_red)
    RedRate=1/(p*(p-1))*sum(abs(cor$cor_mat[upper.tri(cor$cor_mat)]))
  }
  return(RedRate)
}

JaccardRate<-function(data,data_red,k=6){
  gower.distR <- daisy(data_red, metric = c("gower"))
  gower.distF <- daisy(data, metric = c("gower"))
  GDR <- as.matrix(gower.distR)
  GDF <- as.matrix(gower.distF)
  n=nrow(data)
  JR=0
  for(i in 1:n){
    RN<-order(GDR[i,])[2:(k+1)]
    FN<-order(GDF[i,])[2:(k+1)]
    Un=length(union(RN,FN))
    In=length(intersect(RN,FN))
    JR=JR+In/Un
  }
  JR=JR/n
  return(JR)
}
