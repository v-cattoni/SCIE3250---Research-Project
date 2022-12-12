#TASK_2.0
#Estimates parameters for BH and HS for each col using optim. 
#Plots BH vs HS
#Calculates RSS for each col vs BH and HS
#Counts the number of cols which fit BH or HS better

####READING_IN_FILES####
load("/Users/vcatt31/Desktop/Uni/Summer 2022/SCIE3250/Tasks/Data files/RAMCore[asmt][v4.495].rdata")
non_det_cols = read.csv("non_det_cols.csv", header = FALSE)
non_det_cols = non_det_cols[,1]

####CONSTANTS####
years <- seq(1950, 2020, 1) #year of each row of data
num_rows <- length(years) #num rows to iterate through all data
num_cols <- length(non_det_cols) #num cols to iterate through all data

####INITIALISING####
bio_mat <- matrix(NA, num_rows, num_cols) #will store biomass(MT)
catch_mat <- matrix(NA, num_rows, num_cols) #will store catch(MT)
residuals <- matrix(NA, num_cols, 2) #will store rss for 'best-fit' of b_h, h_s, and 
num_b_h <- 0 
num_h_s <- 0

####FUNCTIONS####
#Beverton-Holt
b_h <- function(r, k, X){
  return((r*X)/(1+(r-1)*X/k))
}

#Hockey-Stick
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

#Residual Sum of Squares
rss <- function(X, Y){
  rs <- 0
  for(i in 1:length(X)){
    rs <- rs + (X[i] - Y[i])**2
  }
  return(rs)
}

#Residual Sum of Squares (Beverton-Holt) 
#This is the function which will optimised in terms of r and k
rss_bh = function(p, X, Y){
  r <- p[1]
  k <- p[2]
  return(rss(Y, b_h(r, k, X)))
}

#Residual Sum of Squares (Hockey-Stick)
#This is the function which will optimised in terms of r and k
rss_hs = function(p, X, Y){
  r <- p[1]
  k <- p[2]
  return(rss(Y, h_s(r, k, X)))
}

####SORTING_DATA####
#Creates a matrix of Biomass data and a matrix of Catch data with only
#overlapping entries.
for(i in non_det_cols){
  for(j in 1:num_rows){
    if(FALSE %in% is.na(Bio[,i][j]) & FALSE %in% is.na(MCatch[,i][j])){
      bio_mat[j, match(i, non_det_cols)] <- c(Bio[,i][j])
      catch_mat[j, match(i, non_det_cols)] <- c(MCatch[,i][j])
    }
  }
}

####MAIN####
for(i in 1:length(non_det_cols)){
  bios <- bio_mat[, i][!is.na(bio_mat[,i])]
  catchs <- catch_mat[, i][!is.na(catch_mat[,i])]
  X <- bios[1:(length(bios)-1)] - catchs[1:(length(bios)-1)]
  Y <- bios[2:length(bios)] 
  
  bh_p_guess <- c(mean(Y/X), max(Y))
  hs_p_guess <- c(mean(Y/X), max(Y))
  
  bh_opt <- optim(bh_p_guess, rss_bh, X = X, Y = Y)
  hs_opt <- optim(hs_p_guess, rss_hs, X = X, Y = Y)
  
  bh_opt_r <- bh_opt[[1]][1]
  bh_opt_k <- bh_opt[[1]][2]
  hs_opt_r <- hs_opt[[1]][1]
  hs_opt_k <- hs_opt[[1]][2]
  
  residuals[i,1] <- rss(Y, b_h(bh_opt_r, bh_opt_k, X))
  residuals[i,2] <- rss(Y, h_s(hs_opt_r, hs_opt_k, X))

  # jpeg(file = as.character(i))
  # plot(X, Y,
  #      main = non_det_cols[i],
  #      xlim = c(0, max(X)),
  #      ylim = c(0, max(c(max(Y),
  #                        max(b_h(bh_opt_r, bh_opt_k, X)),
  #                        max(h_s(hs_opt_r, hs_opt_k, X))))))
  # lines(X, b_h(bh_opt_r, bh_opt_k, X), col = 'blue', lwd = 2, type = 'p')
  # lines(X, h_s(hs_opt_r, hs_opt_k, X), col = 'green', lwd = 2, type = 'p')
  # dev.off()
 
  if(residuals[i, 1] == min(residuals[i,])){
    num_b_h = num_b_h + 1   #num_b_h = number of times BH has lower SSE
  } else {
    num_h_s = num_h_s + 1   #num_h_s = number of times HS has lower SSE
  }
}




