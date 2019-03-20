computeDist <- function(ds1, ds2, region_bounds, sigma, bootstrap_n, n_background_1, n_background_2){ #####Call this. 
  ###ds1, ds2 is data
  ###boot_strap is how subsamples to be passed to mmd
  ###n_background_1, n_background_2 how much contrast should be added
  #######normalise data
  
  if (is.null(ds1) | is.null(ds2)) {
    print("Error: Null arguments")
    return(NA)
  }
  
  if (length(ds1) < 2 | length(ds2)< 2) {
    return(NA)
  }
  
  if (bootstrap_n < 2) {
    bootstrap_n1 = length(ds1)
    bootstrap_n2 = length(ds2)
  } else {
    bootstrap_n1 = bootstrap_n
    bootstrap_n2 = bootstrap_n
  }
  
  #### create joint, augument with noise
  ds1_augumented <- createJoint(ds1, n_background_1, region_bounds)
  ds2_augumented <- createJoint(ds2, n_background_2, region_bounds)
  
  sample1 <- ds1_augumented[sample(nrow(ds1_augumented), bootstrap_n1, replace = TRUE), ]
  sample2 <- ds2_augumented[sample(nrow(ds2_augumented), bootstrap_n2, replace = TRUE), ]
  
  result <- runMMD(sample1, sample2, region_bounds, sigma)
  
  return(result)
}

createJoint <- function(ds, n_background, region_bounds){
  rand_noise <- floor(runif(n_background, region_bounds[1], region_bounds[2]))

  ds_fg <- data.frame(positions=ds, obs_type=1)
  ds_bg <- data.frame(positions=rand_noise, obs_type=0)

  ds_augumented<- rbind(ds_fg, ds_bg)
  return(ds_augumented)
}

runMMD <- function(joint_ds1, joint_ds2, region_bounds, sigma){
  a1 = as.integer(joint_ds1[[1]])
  a2 = as.integer(joint_ds1[[2]])
  
  b1 = as.integer(joint_ds2[[1]]) 
  b2 = as.integer(joint_ds2[[2]])
  
  minb = as.integer(region_bounds[1]) 
  maxb = as.integer(region_bounds[2])#region_bounds[2])
  
  d_sigma = as.double(sigma)

  if(!(is.double(d_sigma)))
    stop("Double object")
  
  if(!(is.integer(a2) && 
       is.integer(a1) && 
       is.integer(b1) && 
       is.integer(b2) && 
       is.integer(minb) && 
       is.integer(maxb)))
    stop("Integer object")
  
  a<-.Call("jmmd", a1, a2, b1, b2, minb, maxb, d_sigma)
  return(a)
}