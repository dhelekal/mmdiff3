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