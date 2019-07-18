computeDist <- function(ds1, ds2, region_bounds, sigma, bootstrap_n, n_background_1, n_background_2, maxval, lut){ #####Call this. 
  ###ds1, ds2 is data
  ###boot_strap is how subsamples to be passed to mmd
  ###n_background_1, n_background_2 how much contrast should be added
  #######normalise data
  
  if (is.null(ds1) | is.null(ds2)) {
    print("Error: Null arguments")
    return(NA)
  }
  
  if ( ((length(ds1) < 2) & (length(ds2) < 2)) & ((n_background_1<2 )|(n_background_2<2)) ) {
    return(NA)
  }
  
  #### create joint, augument with noise
  ds1_augumented <- createJoint(ds1, n_background_1, region_bounds)
  ds2_augumented <- createJoint(ds2, n_background_2, region_bounds)
  
  if (bootstrap_n < 2) {
    
    sample1 <- ds1_augumented
    sample2 <- ds2_augumented
    
  } else {
    
    bootstrap_n1 <- bootstrap_n
    bootstrap_n2 <- bootstrap_n
    
    sample1 <- ds1_augumented[sample(nrow(ds1_augumented), bootstrap_n1, replace = TRUE), ]
    sample2 <- ds2_augumented[sample(nrow(ds2_augumented), bootstrap_n2, replace = TRUE), ]
  
  }
  
  m <- NROW(sample1)
  n <- NROW(sample2)
  
  positive1 <- subset(sample1, obs_type==1)
  positive2 <- subset(sample2, obs_type==1)
  
  negative1 <- subset(sample1, obs_type==0)
  negative2 <- subset(sample2, obs_type==0)
  
  kxx_p <- kernelSumSymmetric(positive1, maxval, lut)
  kyy_p <- kernelSumSymmetric(positive2, maxval, lut)
  kxy_p <- kernelSum(positive1, positive2, maxval, lut)
  
  kxx_n <- NROW(negative1)^2
  kyy_n <- NROW(negative2)^2
  kxy_n <- NROW(negative1)*NROW(negative2)
  
  result <- (1/(m*m))*(kxx_p+kxx_n)+(1/(n*n))*(kyy_p+kyy_n)-(2/(m*n))*(kxy_p+kxy_n)
  
  return(sqrt(result))
}

computeDist_stoch <- function(ds1, ds2, region_bounds, sigma, bootstrap_n, n_background_1, n_background_2, maxval, lut){ #####Call this. 
  ###ds1, ds2 is data
  ###boot_strap is how subsamples to be passed to mmd
  ###n_background_1, n_background_2 how much contrast should be added
  #######normalise data
  
  if (is.null(ds1) | is.null(ds2)) {
    print("Error: Null arguments")
    return(NA)
  }
  
  if ( ((length(ds1) < 2) & (length(ds2) < 2)) & ((n_background_1<2 )|(n_background_2<2)) ) {
    return(NA)
  }
  
  #### create joint, augument with noise
  ds1_augumented <- createJoint(ds1, n_background_1, region_bounds)
  ds2_augumented <- createJoint(ds2, n_background_2, region_bounds)
  
  if (bootstrap_n < 2) {
    
    sample1 <- ds1_augumented
    sample2 <- ds2_augumented
    
  } else {
    
    bootstrap_n1 <- bootstrap_n
    bootstrap_n2 <- bootstrap_n
    
    sample1 <- ds1_augumented[sample(nrow(ds1_augumented), bootstrap_n1, replace = TRUE), ]
    sample2 <- ds2_augumented[sample(nrow(ds2_augumented), bootstrap_n2, replace = TRUE), ]
    
  }
  
  m <- NROW(sample1)
  n <- NROW(sample2)
  
  positive1 <- subset(sample1, obs_type==1)
  positive2 <- subset(sample2, obs_type==1)
  
  negative1 <- subset(sample1, obs_type==0)
  negative2 <- subset(sample2, obs_type==0)
  
  kxx_p <- kernelSumSymmetric(positive1, maxval, lut)
  kyy_p <- kernelSumSymmetric(positive2, maxval, lut)
  kxy_p <- kernelSum(positive1, positive2, maxval, lut)
  
  kxx_n <- kernelSumSymmetric(negative1, maxval, lut) 
  kyy_n <- kernelSumSymmetric(negative2, maxval, lut)
  kxy_n <- kernelSum(negative1, negative2, maxval, lut)
  
  result <- (1/(m*m))*(kxx_p+kxx_n)+(1/(n*n))*(kyy_p+kyy_n)-(2/(m*n))*(kxy_p+kxy_n)
  
  return(sqrt(result))
}

createJoint <- function(ds, n_background, region_bounds){
  if (NROW(ds)>0) {
    ds_fg <- data.frame(positions=ds, obs_type=1)
  } else {
    ds_fg <- data.frame()
  }
  
  if (n_background > 0) {
    rand_noise <- runif(n_background, min = region_bounds[1], max = region_bounds[2])
    ds_bg <- data.frame(positions=rand_noise, obs_type=0)
    ds_augumented<- rbind(ds_fg, ds_bg)
  } else {
    ds_augumented<- ds_fg
  }

  return(ds_augumented)
}

kernelSum <- function(joint_ds1, joint_ds2, maxVal, lut){
  a1 <- as.integer(joint_ds1[["positions"]])
  a2 <- as.integer(joint_ds1[["obs_type"]])
  
  b1 <- as.integer(joint_ds2[["positions"]]) 
  b2 <- as.integer(joint_ds2[["obs_type"]])
  
  maxVali <- as.integer(maxVal)
  lutd <- as.double(lut)
  
  if(!(is.double(lutd)))
    stop("Double object")
  
  if(!(is.integer(a2) && 
       is.integer(a1) && 
       is.integer(b1) && 
       is.integer(b2) && 
       is.integer(maxVali)))
    stop("Integer object")
  
  a<-.Call("kernel_sum", a1, a2, b1, b2, maxVali, lutd)
  return(a)
}

kernelSumSymmetric <- function(joint_ds1, maxVal, lut){
  a1 <- as.integer(joint_ds1[["positions"]])
  a2 <- as.integer(joint_ds1[["obs_type"]])

  maxVali <- as.integer(maxVal)
  lutd <- as.double(lut)
  
  if(!(is.double(lutd)))
    stop("Double object")
  
  if(!(is.integer(a2) && 
       is.integer(a1) && 
       is.integer(maxVali)))
    stop("Integer object")
  
  a<-.Call("kernel_sum_symmetric", a1, a2, maxVali, lutd)
  return(a)
}

buildMMDLUT <- function(maxVal,sigma){
  
  imax <- as.integer(maxVal)
  dsigma <- as.double(sigma)

  if(!(is.integer(imax)))
    stop("Integer object")
  
  if(!(is.double(dsigma)))
    stop("Double object")
  
  lut <- .Call("kernel_lut", imax, dsigma) 
  return(lut)
}

isParallel <- function(){
  result <- .Call("is_parallel") 
  return(result)
}