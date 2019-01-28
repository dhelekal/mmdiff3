computeDist <- function(ds1, ds2, bootstrap_n, n_background){
  #######normalise data
  ds1_r <-rescale(ds1)
  ds2_r <-rescale(ds2)
  
  #### estimate sigma here
  sigma_est <- 1.0
  #### create joint, augument with noise
  ds1_augumented <- ds_fg
  ds2_augumented <- ds_bg
  
  sample1 <- ds1_augumented[sample(nrow(df), bootstrap_n, replace = TRUE), ]
  sample2 <- ds2_augumented[sample(nrow(df), bootstrap_n, replace = TRUE), ]
  
  result <- runMMD(sample1, sample2, sigma_est)
  
  return(result)
}

createJoint <- function(ds, n_background){
  rand_noise <- runif(n_background, min = 0, max = 1)
  
  ds_fg <- data.frame(positions=ds, obs_type=1)
  ds_bg <- data.frame(positions=rand_noise, obs_type=0)

  ds_augumented<- merge(ds_fg, ds_bg)
  return(ds_augumented)
}



runMMD <- function(joint_ds1, joint_ds2, sigma){
  
  a1 = as.double(joint_ds1[,1])
  a2 = as.integer(joint_ds1[,2])
  
  b1 = as.double(joint_ds2[,1]) 
  b2 = as.integer(joint_ds2[,2])
  
  d_sigma = as.double(sigma)

  if(!(is.double(a1) || is.double(b1) || is.double(d_sigma)))
    stop("Double object")
  
  if(!(is.integer(a2) || is.integer(b2)))
    stop("Integer object")
  
  a<-.Call("jmmd", a1, a2, b1, b2, d_sigma)
  return(a)
}

rescale <- function(x){
  return((x-mean(x))/(max(x)-min(x)))
}