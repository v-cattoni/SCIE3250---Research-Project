#TASK_1.1
#Create biomass time series plots for all species with data

####READING_IN_FILES####
load("/Users/vcatt31/Desktop/Uni/Summer 2022/SCIE3250/Tasks/Data files/RAMCore[asmt][v4.495].rdata")

####CONSTANTS####
years <- seq(1950, 2020, 1)

####INITIALISING####
are_there_only_nas <- rep(NA, length(Bio[,1]))

####MAIN####
#goes through each column in 'Bio.csv' and determines if it is empty
#(this is to avoid an error that occurs when plotting empty vectors)
for(i in 1:length(Bio[1,])){
  for(j in 1:length(Bio[,i])){
    are_there_only_nas[j] <- is.na(Bio[,i][j])
  }
  #plots the non-empty columns of 'Bio.csv'
  if(FALSE %in% are_there_only_nas){
    # jpeg(file = as.character(i)) #saves each plot as the column in 'Bio.csv'
    plot(years, Bio[,i], main = i) #sets title to be the column in 'Bio.csv'
    # dev.off()
  }
}
