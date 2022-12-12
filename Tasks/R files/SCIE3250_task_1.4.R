#TASK_1.5
#Plot Bio_{t}-Mcatch_{t} vs Bio_{t+1}. 

####READING_IN_FILES####
load("/Users/vcatt31/Desktop/Uni/Summer 2022/SCIE3250/Tasks/Data files/RAMCore[asmt][v4.495].rdata")

####CONSTANTS####
years <- seq(1950, 2020, 1)
num_rows <- length(years)

####REMOVING_BAD_DATA####
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
non_det_cols = non_det_cols[c(-116, -182, -315)]

####CONSTANTS####
num_cols <- length(non_det_cols)

####INITIALISING####
bio_mat <- matrix(NA, num_rows, num_cols)
catch_mat <- matrix(NA, num_rows, num_cols)

####SORTING_DATA####
for(i in non_det_cols){
  for(j in 1:num_rows){
    if(FALSE %in% is.na(Bio[,i][j]) & FALSE %in% is.na(MCatch[,i][j])){
      bio_mat[j, match(i, non_det_cols)] <- c(Bio[,i][j])
      catch_mat[j, match(i, non_det_cols)] <- c(MCatch[,i][j])
    }
  }
}

####MAIN####
for(i in 1:num_cols){
  bios = bio_mat[, i][!is.na(bio_mat[,i])]
  catchs = catch_mat[, i][!is.na(catch_mat[,i])]
  r_t = bios[1:(length(bios)-1)] - catchs[1:(length(bios)-1)]
  r_t1 = bios[2:length(bios)] 
  # jpeg(file = as.character(i))
  plot(r_t, r_t1, 
       main = non_det_cols[i],
       xlim = c(0, max(r_t)),
       ylim = c(0, max(r_t1)))
  # dev.off()
}


