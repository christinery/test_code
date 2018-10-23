################################################################################
## This file include some elementary computation for array mentioned in the ####
## file "mainQRnoweight" for paper "Copula-based M-estimate for FVCM" ##########
################################################################################


################################################################################
##This function is used to compute the difference for each element of a array###
##minus a same matrix, and get the mean of these difference's square ###########
################################################################################
matrixdiff = function(M1, M2){
 ##M1 is array,M2 is matrix,the first and two dimension of M1 is the same with dim(M2)
 matrixnum = dim(M1)[3]
 ##comdiffM store all difference square value with replacement sample times B
 ##number of "matrixnum" matrix
 comdiffM = array(0, c(dim(M2)[1], dim(M2)[2], matrixnum))
 ##"comdiffsum" used to store sum of difference among these num. B matrix
 comdiffsum = matrix(0, dim(M2)[1], dim(M2)[2])
 for (k in 1:matrixnum) {
   comdiffM[,,k] = (M1[,,k] - M2)^2
   comdiffsum = comdiffsum + comdiffM[,,k]
 }
 comdiffmean = comdiffsum/(matrixnum-1)
 return(comdiffmean) ##dim(M2)[1] * dim(M2)[2] matrix
}
###############################################################################
########################## end ################################################
###############################################################################


#############################################################################
##This function is used to compute the mean for a array by sum each element##
##or sqrt element (matrix) of array##########################################
#############################################################################
matrixsqmean = function(A,I){
  ##A is a array
  ##I takes value 1 or 0, 0 denote take sqrt firstly then compute mean, 1 denote
  ##compute the mean for a array by sum each element of array
  index = dim(A)[3]
  ASUM = matrix(0, dim(A)[1], dim(A)[2])
  if(I==0){
  for (l in 1:index) { ASUM = ASUM + sqrt(A[,,l]) }
  return(ASUM/index)
  } else{
    for (l in 1:index) {
      ASUM = ASUM + (A[,,l])
    }
    return(ASUM/index) ##retrun "dim(A)[1] * dim(A)[2]" matrix
  }
}
###############################################################################
########################## end ################################################
###############################################################################


#################################################################################
##This function is used to find the empirical coverage probability for function##
##coefficient discreted by grid points###########################################
#################################################################################
matrixecp = function(A1, A2, A0){
  ##A1, A2 are array and with the same dimension, A0 is matrix with the same dimension
  ##for the first and two element of dim(A1)
  index1 = dim(A1)[1]
  index2 = dim(A1)[2]
  index3 = dim(A1)[3]
  if (prod(dim(A1)==dim(A2))==0)
    stop("The dimension of the array must be the same!")
    ##compute the left and right point for interval
    intervall = A1 - 1.96*sqrt(A2)
    intervalr = A1 + 1.96*sqrt(A2)
    A = array(A0, c(index1, index2, index3))
    ##the difference between left or right of interval and A
    intervalld = intervall - A
    intervalrd = intervalr - A
    ##product of difference, return number "index3" array, matrix dimension is "index1 * index2" 
    intervalpro = intervalld * intervalrd
    ##identify the positive or negative of each element of the array
    matrixI = (intervalpro<=0)
    ##element of array is added, return one "index1 * index2" matrix
    TOT = matrix(0, index1, index2)
    for (k in 1:index3) {
      TOT = TOT + matrixI[,,k]
    }
    return(TOT/index3)
}
###############################################################################
########################## end ################################################
###############################################################################


#############################################################################
##This function is used to compute the row mean for each element (matrix)####
##of array.##################################################################
#############################################################################
arraymatrixmean = function(A){
  ##A is a array
  index = dim(A)[3]
  COMMEAN = NULL
  for (l in 1:index) {
    COMMEAN = rbind(COMMEAN, apply(A[,,l], 1, mean))
  }
  return(COMMEAN)
}
###############################################################################
########################## end ################################################
###############################################################################














