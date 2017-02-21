#!/usr/bin/env Rscript

#readLines(con = file("stdin"))

#options(rgl.useNULL = TRUE) 

source("~/Box Sync/Research/Time-Series/Code/TS_Functions.R")
#source("~/Dropbox/Time-Series/Code/TS_Functions.R")
# source("C:/Users/XeroXtancy/Dropbox/Research/Time-Series/Code/TS_Functions.R")

##############################################################################################
# LOAD IN DATA
##############################################################################################

# files <- c("50words", "Adiac", "ArrowHead", "Beef","BeetleFly","BirdChicken","Car","CBF","CinC_ECG_torso","Coffee","DiatomSizeReduction",
#   "Earthquakes","ECGFiveDays","FaceFour","FISH","Gun_Point","Ham","Haptics","Herring","InlineSkate",
# "Lighting2","Lighting7","MALLAT","Meat","OliveOil","Plane","ShapeletSim","Symbols","ToeSegmentation1",
# "ToeSegmentation2","Trace","Wine","WordsSynonyms","Worms")
# 
# 
# for (k in files) {

fname = "ECGFiveDays"
path = paste("/Users/johnkimnguyen/Box Sync/Research/Time-Series/Code/UCR_To_Do_w40_a5/", fname,"/", sep="")
#path = paste("~/Dropbox/Research/Time-Series/Code/UCR_To_Do_w40_a5/", fname,"/", sep="")
#path = paste("~/Dropbox/Time-Series/Code/UCR_To_Do_w40_a5/", fname,"/", sep="")
#path = paste("~/Dropbox/Research/Time-Series/Code/UCR_To_Do_w40_a5/", fname,"/", sep="")
#path = paste("C:/Users/XeroXtancy/Dropbox/Research/Time-Series/Code/Unused_Data/", fname,"/", sep="")

setwd(path)

data <- read.table(paste(fname,'_TRAIN', sep=""), sep=',')
data <- as.data.frame(t(data)) # Transpose the data since # seq is num of rows
data <-  data[-1,] # Drop the first row since we dont need class label
data <- data.frame(TIME = 1:nrow(data), data)

data_ts <- read.zoo(data) # Convert data to zoo time-series
colnames(data_ts) <- gsub("V","X",colnames(data_ts))

# Plot the time-series | Be sure to rename the columns
# g <-autoplot(data_ts[,c(1:10)]) +
#   theme_bw() +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),legend.position="none",
#         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),plot.background=element_blank())
# g

##############################################################################################
# SET PARAMETERS FOR SAX
##############################################################################################

w = 50 # Number of breakpoints = length of sequence
alpha = 20 # Size of alphabet

##############################################################################################
# GET AND SAVE DISTANCE MATRIX
##############################################################################################
# 
# name <- list_of_SAX(data_ts, L, alpha)
# lapply(name, write, "list_of_SAX_w40_a5.txt", append=TRUE, ncolumns=1000)
# 
# #Rprof("file.out")
# 
# Get the Distance Matrix and Save them
# EUCLDIST <- eucldist(data_ts)
# write.csv(EUCLDIST, file='EUCLDIST.csv')
# }
# DTWDIST <- dtwdist(data_ts)
# write.csv(DTWDIST, file='DTWDIST.csv')
# MINDIST <- mindist(data_ts, w, alpha)
# write.csv(MINDIST, file='MINDIST.csv')
#EDITDIST <- edit(data_ts, alpha)
# write.csv(EDITDIST, file='EDITDIST.csv')
# QGRAMDIST <- qgram_dist(data_ts, alpha)
# write.csv(QGRAMDIST, file='QGRAM.csv')
# APXDIST <- apxdist(data_ts, w, alpha)
# write.csv(APXDIST, file='APXDIST.csv')


# EUCL2 = euclidean(data_ts, 2, alpha)
# write.csv(EUCL2, file='EUCL-2.csv')
# EUCL3 = euclidean(data_ts, 3, alpha)
# write.csv(EUCL3, file='EUCL-3.csv')
# 
# KLMEAN2 = kldiv(data_ts, 2, alpha)
# write.csv(KLMEAN2, file='KLMEAN-2.csv')
# KLMEAN3 = kldiv(data_ts, 3, alpha)
# write.csv(KLMEAN3, file='KLMEAN-3.csv')
KSDIST-2 <- ksdist(data_ts, 2, alpha)
write.csv(KSDIST-2, file='KSDIST-2.csv')
KSDIST-3 <- ksdist(data_ts, 3, alpha)
write.csv(KSDIST, file='KSDIST-3.csv')
JSDIST-2 <- jsdiv(data_ts, 2, alpha)
write.csv(JSDIST, file='JSDIST-2.csv')
JSDIST-3 <- jsdiv(data_ts, 3, alpha)
write.csv(JSDIST, file='JSDIST-3.csv')

#Rprof(NULL)
#summaryRprof("file.out")

# USE THIS IF THE COMPUTATION WAS ALREADY SAVED!
APXDIST = read.table("APXDIST.csv", sep=",", header=TRUE)[-1]
MINDIST = read.table("MINDIST.csv", sep=",", header=TRUE)[-1]
DTWDIST = read.table("DTWDIST.csv", sep=",", header=TRUE)[-1]
#EDITDIST = read.table("EDITDIST.csv", sep=",", header=TRUE)[-1]
# #QGRAMDIST = read.table("QGRAM.csv", sep=",", header=TRUE)[-1]
#EUCL2 = read.table('EUCL-2.csv', sep=",", header=TRUE)[-1]
#EUCL3 = read.table('EUCL-3.csv', sep=",", header=TRUE)[-1]
KLMEAN2 = read.table('KLMEAN-2.csv', sep=",", header=TRUE)[-1]
KLMEAN3 = read.table('KLMEAN-3.csv', sep=",", header=TRUE)[-1]
# KLMEAN4 = read.table('KLMEAN-4.csv', sep=",", header=TRUE)[-1]
JSDIST = read.table('JSDIST-2.csv', sep=",", header=TRUE)[-1]
KSDIST = read.table('KSDIST.csv', sep=",", header=TRUE)[-1]


##############################################################################################
# COMPUTE F-MEAUSRES WITH VARYING CLUSTER
##############################################################################################

min_fmeas_pam <- list()
#edit_fmeas_pam <- list()
#qgram_fmeas_pam <- list()
apx_fmeas_pam <- list()
# eucl2_fmeas_pam <- list()
# eucl3_fmeas_pam <- list()
klmean2_fmeas_pam <- list()
klmean3_fmeas_pam <- list()
#klmean4_fmeas_pam <- list()

# mindist_clusters <- hclust(as.dist(MINDIST))
# apxdist_clusters <- hclust(as.dist(APXDIST))
# dtwclusters <- hclust(as.dist(DTWDIST))
# KLmean_clusters <- hclust(as.dist(KLMEAN2))
# 
# par(mfrow = c(2,2), mai = c(0.1, 0.1, 0.5, 0.5))
# plot(dtwclusters, main='DTWDIST', axes=FALSE)
# plot(mindist_clusters, main='MINDIST', axes=FALSE)
# plot(apxdist_clusters, main='APXDIST', axes=FALSE)
# plot(KLmean_clusters, main='KLAVG2', axes=FALSE)

for (k in 2:4) {
  # Cluster with PAM
  dtw_clust <- pam(DTWDIST, k)
  min_clust <- pam(MINDIST, k)
  #edit_clust <- pam(EDITDIST, k)
  #qgram_clust <- pam(QGRAMDIST, k)
  apx_clust <- pam(APXDIST, k)
  #eucl2_clust <- pam(EUCL2, k)
  #eucl3_clust <- pam(EUCL3, k)
  klmean2_clust <- pam(KLMEAN2, k)
  klmean3_clust <- pam(KLMEAN3, k)
  #klmean4_clust <- pam(KLMEAN4, k)

  # Get the F-measures
  min_fmeas_pam[[k]] <- fmeasure(min_clust$clustering, dtw_clust$clustering)
  #edit_fmeas_pam[[k]] <- fmeasure(edit_clust$clustering, dtw_clust$clustering)
  #qgram_fmeas_pam[[k]] <- fmeasure(qgram_clust$clustering, dtw_clust$clustering)
  apx_fmeas_pam[[k]] <- fmeasure(apx_clust$clustering, dtw_clust$clustering)
  #eucl2_fmeas_pam[[k]] <- fmeasure(eucl2_clust$clustering, dtw_clust$clustering)
  #eucl3_fmeas_pam[[k]] <- fmeasure(eucl3_clust$clustering, dtw_clust$clustering)
  klmean2_fmeas_pam[[k]] <- fmeasure(klmean2_clust$clustering, dtw_clust$clustering)
  klmean3_fmeas_pam[[k]] <- fmeasure(klmean3_clust$clustering, dtw_clust$clustering)
  #klmean4_fmeas_pam[[k]] <- fmeasure(klmean4_clust$clustering, dtw_clust$clustering)
}

df <- data.frame(unlist(min_fmeas_pam),unlist(apx_fmeas_pam),#unlist(edit_fmeas_pam),
                 #unlist(eucl2_fmeas_pam), unlist(eucl3_fmeas_pam),
                 unlist(klmean2_fmeas_pam), unlist(klmean3_fmeas_pam))#, unlist(klmean4_fmeas_pam))

colnames(df) <- c("mindist","apxdist", "klmean2", "klmean3")#, "klmean4")#, "cosine2", "cosine3")
df<- data.frame(Num_Cluster = 1:nrow(df), df)
df.m <- melt(df, id.vars = "Num_Cluster")
colMeans(df)

ggplot(df.m) + geom_line(aes(x = Num_Cluster, y = value, colour = variable))
# #df["klmean4"]
# #colMeans(df)
# write.csv(t(df), file='cluster_size_result.csv')
pamk(DTWDIST)$nc
# 
# 
# 
# min_jaccard_pam <- list()
# edit_jaccard_pam <- list()
# apx_jaccard_pam <- list()
# eucl2_jaccard_pam <- list()
# eucl3_jaccard_pam <- list()
# klmean2_jaccard_pam <- list()
# klmean3_jaccard_pam <- list()
# 
# for (k in 2:16) {
#   # Cluster with PAM
#   dtw_clust <- pam(DTWDIST, k)
#   min_clust <- pam(MINDIST, k)
#   edit_clust <- pam(EDITDIST, k)
#   apx_clust <- pam(APXDIST, k)
#   eucl2_clust <- pam(EUCL2, k)
#   eucl3_clust <- pam(EUCL3, k)
#   klmean2_clust <- pam(KLMEAN2, k)
#   klmean3_clust <- pam(KLMEAN3, k)
# 
#   # Get the F-measures
#   min_jaccard_pam[[k]] <- extCriteria(min_clust$clustering, dtw_clust$clustering, "Rand")
#   edit_jaccard_pam[[k]] <- extCriteria(edit_clust$clustering, dtw_clust$clustering, "Rand")
#   apx_jaccard_pam[[k]] <- extCriteria(apx_clust$clustering, dtw_clust$clustering, "Rand")
# 
#   eucl2_jaccard_pam[[k]] <- extCriteria(eucl2_clust$clustering, dtw_clust$clustering, "Rand")
#   eucl3_jaccard_pam[[k]] <- extCriteria(eucl3_clust$clustering, dtw_clust$clustering, "Rand")
# 
#   klmean2_jaccard_pam[[k]] <- extCriteria(klmean2_clust$clustering, dtw_clust$clustering, "Rand")
#   klmean3_jaccard_pam[[k]] <- extCriteria(klmean3_clust$clustering, dtw_clust$clustering, "Rand")
# }
# df <- data.frame(unlist(min_jaccard_pam),unlist(apx_jaccard_pam),unlist(edit_jaccard_pam),
#                  unlist(eucl2_jaccard_pam), unlist(eucl3_jaccard_pam),
#                  unlist(klmean2_jaccard_pam), unlist(klmean3_jaccard_pam))
# colnames(df) <- c("mindist","apxdist", "eed", "eucl2","eucl3", "klmean2", "klmean3")#, "cosine2", "cosine3")
# df<- data.frame(Num_Cluster = 1:nrow(df), df)
# #df.m <- melt(df, id.vars = "Num_Cluster")
# #ggplot(df.m) + geom_line(aes(x = Num_Cluster, y = value, colour = variable))
# colMeans(df)
# write.csv(t(df), file='cluster_size_result_rand.csv')



##############################################################################################
# COMPUTE F-MEAUSRES WITH VARYING SEQUENCE LENGTH
##############################################################################################
# 
# min_fmeas_pam <- list()
# edit_fmeas_pam <- list()
# apx_fmeas_pam <- list()
# eucl2_fmeas_pam <- list()
# eucl3_fmeas_pam <- list()
# klmean2_fmeas_pam <- list()
# klmean3_fmeas_pam <- list()
# 
# alpha = 25
# 
# path = "/Users/johnkimnguyen/Dropbox/Research/Time-Series/Code/Varied_Sequence_Size/BirdChicken/"
# setwd(path)
# 
# DTWDIST = read.table("DTWDIST.csv", sep=",", header=TRUE)[-1]
# 
# for (w in seq(50,300,10)) {
# 
#   # fname <- list_of_SAX(data_ts, L, alpha)
#   # lapply(fname, write, paste("list_of_SAX_a", alpha, "_w", w, ".txt", sep=""), append=TRUE, ncolumns=1000)
# # 
# #   MINDIST <- mindist(data_ts, w, alpha)
# #   EDITDIST <- edit(data_ts, alpha)
# #   APXDIST <- apxdist(data_ts, w, alpha)
# # 
# #   write.csv(MINDIST, file=paste('MINDIST', w, ".csv", sep=""))
# #   write.csv(APXDIST, file=paste('APXDIST', w, ".csv", sep=""))
# #   write.csv(EDITDIST, file=paste('EDITDIST', w, ".csv", sep=""))
# 
#   DTWDIST = read.table("DTWDIST.csv", sep=",", header=TRUE)[-1]
#   APXDIST = read.table(paste('APXDIST', w, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   MINDIST = read.table(paste('MINDIST', w, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   EDITDIST = read.table(paste('EDITDIST', w, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   EUCL2 = read.table(paste('EUCL-2w', w, ".csv"), sep=",", header=TRUE)[-1]
#   EUCL3 = read.table(paste('EUCL-2w', w, ".csv"), sep=",", header=TRUE)[-1]
#   KLMEAN2 = read.table(paste('KLMEAN-2w', w, ".csv"), sep=",", header=TRUE)[-1]
#   KLMEAN3 = read.table(paste('KLMEAN-3w', w, ".csv"), sep=",", header=TRUE)[-1]
# 
#   dtw_clust <- pam(DTWDIST, 7)
#   min_clust <- pam(MINDIST, 7)
#   edit_clust <- pam(EDITDIST, 7)
#   apx_clust <- pam(APXDIST, 7)
#   eucl2_clust <- pam(EUCL2, 7)
#   eucl3_clust <- pam(EUCL3, 7)
#   klmean2_clust <- pam(KLMEAN2, 7)
#   klmean3_clust <- pam(KLMEAN3, 7)
# 
#   # Get the F-measures
#   min_fmeas_pam[[w]] <- fmeasure(min_clust$clustering, dtw_clust$clustering)
#   edit_fmeas_pam[[w]] <- fmeasure(edit_clust$clustering, dtw_clust$clustering)
#   apx_fmeas_pam[[w]] <- fmeasure(apx_clust$clustering, dtw_clust$clustering)
#   eucl2_fmeas_pam[[w]] <- fmeasure(eucl2_clust$clustering, dtw_clust$clustering)
#   eucl3_fmeas_pam[[w]] <- fmeasure(eucl3_clust$clustering, dtw_clust$clustering)
#   klmean2_fmeas_pam[[w]] <- fmeasure(klmean2_clust$clustering, dtw_clust$clustering)
#   klmean3_fmeas_pam[[w]] <- fmeasure(klmean3_clust$clustering, dtw_clust$clustering)
# 
# }
# 
# df <- data.frame(unlist(min_fmeas_pam),unlist(apx_fmeas_pam),unlist(edit_fmeas_pam),
#                  unlist(eucl2_fmeas_pam), unlist(eucl3_fmeas_pam),
#                  unlist(klmean2_fmeas_pam), unlist(klmean3_fmeas_pam))#,unlist(cosine2_fmeas_pam), unlist(cosine3_fmeas_pam))
# colnames(df) <- c("mindist","apxdist", "eed", "eucl2","eucl3", "klmean2", "klmean3")#, "cosine2", "cosine3")
# df<- data.frame(SequenceSize = seq(50,300,10), df)
# 
# df.m <- melt(df, id.vars = "SequenceSize")
# ggplot(df.m) + geom_line(aes(x = SequenceSize, y = value, colour = variable))
# colMeans(df)
