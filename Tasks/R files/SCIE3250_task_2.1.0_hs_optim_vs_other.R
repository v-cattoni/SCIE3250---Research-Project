#TASK_2.0
#Estimating parameters and calculating sum of squares.


####READING_IN_FILES####
load("/Users/vcatt31/Desktop/Uni/Summer 2022/SCIE3250/Tasks/Data files/RAMCore[asmt][v4.495].rdata")


####CONSTANTS####
#sorts through data removing any empty columns
non_det_cols = c()
are_there_only_nas <- rep(NA, length(Bio[,1]))
for(i in 1:length(Bio[1,])){
  for(j in 1:length(Bio[,i])){
    are_there_only_nas[j] <- is.na(Bio[,i][j])
  }
  if(FALSE %in% are_there_only_nas){
    non_det_cols = append(non_det_cols, i)
  }
}

#removes any colums with trouble process (negs)
non_det_cols = non_det_cols[c(-20, -22, -27, -33, -35,
                              -38, -51, -60, -62, -70,
                              -75, -88, -108, -109, -116,
                              -122, -124, -135, -139, -141,
                              -144, -146, -177, -182, -207,
                              -213, -228, -237, -240, -241,
                              -257, -263, -264, -272, -277,
                              -284, -289, -299, -303, -304,
                              -315, -321, -325, -327, -330,
                              -336, -337, -339, -340, -347,
                              -352, -356, -362, -363, -365,
                              -366, -367, -387, -388, -398,
                              -413, -420, -423, -438)]
years <- seq(1950, 2020, 1) #year of each row of data
accuracy_of_k <- 100 #for estimating b_h parameters
accuracy_of_r <- 0.001
num_rows <- length(years) #num rows to iterate through all data
num_cols <- length(non_det_cols) #num cols to iterate through all data

####INITIALISING####
bio_mat <- matrix(NA, num_rows, num_cols) #will store biomass(MT)
catch_mat <- matrix(NA, num_rows, num_cols) #will store catch(MT)
residuals <- matrix(NA, num_cols, 2) #will store rss for 'best-fit' of b_h, h_s, and 
num_h_s_optim <- 0 
num_h_s_other <- 0
num_h_s_same <- 0


####FUNCTIONS####

h_s <- function(r, k, X){
  h <- c()
  for(i in X){
    if(i < k/r){
      h <- append(h, r*i)
    } else {
      h <- append(h, k)
    }
  }
  return(h)
}


rss <- function(X, Y){
  rs <- 0
  for(i in 1:length(X)){
    rs <- rs + (X[i] - Y[i])**2
  }
  return(rs)
}

rss_hs = function(p, X, Y){
  r <- p[1]
  k <- p[2]
  return(rss(Y, h_s(r, k, X)))
}

for(i in non_det_cols){
  for(j in 1:num_rows){
    if(FALSE %in% is.na(Bio[,i][j]) & FALSE %in% is.na(MCatch[,i][j])){
      bio_mat[j, match(i, non_det_cols)] <- c(Bio[,i][j])
      catch_mat[j, match(i, non_det_cols)] <- c(MCatch[,i][j])
    }
  }
}

h_s_parameter_est <- function(X, Y){
  r_est = sum(X*Y)/sum(X**2)
  k_est = sum(Y)/length(Y)
  for(i in 1:length(X)){
    parms = alpha(X, Y, r_est, k_est)
    r_est = parms[1]
    k_est = parms[2]
  }
  
  return(c(r_est, k_est))
}
alpha = function(X, Y, r_est, k_est){
  XL = c()
  YL = c()
  XU = c()
  YU = c()
  for(i in 1:length(X)){
    if(X[i] < k_est/r_est){
      XL = append(XL, X[i])
      YL = append(YL, Y[i])
    } else{
      XU = append(XU, X[i])
      YU = append(YU, Y[i])
    }
  }
  r_est = sum(XL*YL)/sum(XL**2)
  k_est = sum(YU)/length(YU)
  return(c(r_est, k_est))
}



####MAIN####
for(i in 1:num_cols){
  bios <- bio_mat[, i][!is.na(bio_mat[,i])]
  catchs <- catch_mat[, i][!is.na(catch_mat[,i])]
  X <- bios[1:(length(bios)-1)] - catchs[1:(length(bios)-1)]
  Y <- bios[2:length(bios)] 
  
  h_s_param_est = h_s_parameter_est(X, Y)
  hs_opt_r <- h_s_param_est[1]
  hs_opt_k <- h_s_param_est[2]
  
  hs_p_guess <- c(mean(Y/X), max(Y))
  hs_opt <- optim(hs_p_guess, rss_hs, X = X, Y = Y)
  hs_optim_r <- hs_opt[[1]][1]
  hs_optim_k <- hs_opt[[1]][2]
  
  
  residuals[i,1] <- rss(Y, h_s(hs_opt_r, hs_opt_k, X))
  residuals[i,2] <- rss(Y, h_s(hs_optim_r, hs_optim_k, X))

  # jpeg(file = as.character(i))
  # plot(X, Y, main = non_det_cols[i])
  # lines(X, h_s(hs_opt_r, hs_opt_k, X), col = 'blue', lwd = 2, type = 'p')
  # lines(X, h_s(hs_optim_r, hs_optim_k, X), col = 'green', lwd = 2, type = 'p')
  # dev.off()
  
  if(1.01*residuals[i, 1] > (residuals[i, 2]) & 1.01*(residuals[i, 2]) > residuals[i, 1]){
    num_h_s_same = num_h_s_same + 1
  } else if(1.01*residuals[i, 1] < residuals[i, 2]){
    num_h_s_other = num_h_s_other + 1
  } else if (residuals[i, 1] > 1.01*residuals[i, 2]){
    num_h_s_optim = num_h_s_optim + 1
  }
}
print(num_h_s_same)
print(num_h_s_other)
print(num_h_s_optim)
