source("~/Dropbox/Research/Time-Series/Code/TS_Functions.R")
#source("C:/Users/XeroXtancy/Dropbox/Research/Time-Series/Code/TS_Functions.R")

##############################################################################################
# LOAD IN DATA
##############################################################################################

fname = "ShapeletSim"

path = paste("/Users/johnkimnguyen/Dropbox/Research/Time-Series/Code/Varied_Sequence_Size/", fname,"/", sep="")

setwd(path)

data <- read.table(paste(fname,'_TRAIN', sep=""), sep=',')
data <- as.data.frame(t(data)) # Transpose the data since # seq is num of rows
data <-  data[-1,] # Drop the first row since we dont need class label
data <- data.frame(TIME = 1:nrow(data), data)
data_ts <- read.zoo(data) # Convert data to zoo time-series

DTWDIST = read.table("DTWDIST.csv", sep=",", header=TRUE)[-1]
alpha = 8


##############################################################################################
# COMPUTE F-MEASURES WITH VARYING SEQUENCE LENGTH
##############################################################################################

min_fmeas_pam <- list()
edit_fmeas_pam <- list()
apx_fmeas_pam <- list()
#eucl2_fmeas_pam <- list()
#eucl3_fmeas_pam <- list()
klmean2_fmeas_pam <- list()
klmean3_fmeas_pam <- list()

k = 2

#### RESUME AT 140
for (w in seq(10,150,10)) {
# fname <- list_of_SAX(data_ts, L, alpha)
# lapply(fname, write, paste("list_of_SAX_a", alpha, "_w", w, ".txt", sep=""), append=TRUE, ncolumns=1000)
#
# MINDIST <- mindist(data_ts, w, alpha)
# EDITDIST <- edit(data_ts, alpha)
# APXDIST <- apxdist(data_ts, w, alpha)
#   write.csv(MINDIST, file=paste('MINDIST', w, ".csv", sep=""))
#   write.csv(APXDIST, file=paste('APXDIST', w, ".csv", sep=""))
#   write.csv(EDITDIST, file=paste('EDITDIST', w, ".csv", sep=""))
# print(w)
# }

  DTWDIST = read.table("DTWDIST.csv", sep=",", header=TRUE)[-1]
  APXDIST = read.table(paste('APXDIST', w, ".csv", sep=""), sep=",", header=TRUE)[-1]
  MINDIST = read.table(paste('MINDIST', w, ".csv", sep=""), sep=",", header=TRUE)[-1]
  EDITDIST = read.table(paste('EDITDIST', w, ".csv", sep=""), sep=",", header=TRUE)[-1]
#  EUCL2 = read.table(paste('EUCL-2w', w, ".csv", sep=""), sep=",", header=TRUE)[-1]
#  EUCL3 = read.table(paste('EUCL-2w', w, ".csv", sep=""), sep=",", header=TRUE)[-1]
  KLMEAN2 = read.table(paste('KLMEAN-2w', w, ".csv", sep=""), sep=",", header=TRUE)[-1]
  KLMEAN3 = read.table(paste('KLMEAN-3w', w, ".csv", sep=""), sep=",", header=TRUE)[-1]

  dtw_clust <- pam(DTWDIST, k)
  min_clust <- pam(MINDIST, k)
  edit_clust <- pam(EDITDIST, k)
  apx_clust <- pam(APXDIST, k)
  #eucl2_clust <- pam(EUCL2, k)
  #eucl3_clust <- pam(EUCL3, k)
  klmean2_clust <- pam(KLMEAN2, k)
  klmean3_clust <- pam(KLMEAN3, k)

  # Get the F-measures
  min_fmeas_pam[[w]] <- fmeasure(min_clust$clustering, dtw_clust$clustering)
  edit_fmeas_pam[[w]] <- fmeasure(edit_clust$clustering, dtw_clust$clustering)
  apx_fmeas_pam[[w]] <- fmeasure(apx_clust$clustering, dtw_clust$clustering)
  #eucl2_fmeas_pam[[w]] <- fmeasure(eucl2_clust$clustering, dtw_clust$clustering)
  #eucl3_fmeas_pam[[w]] <- fmeasure(eucl3_clust$clustering, dtw_clust$clustering)
  klmean2_fmeas_pam[[w]] <- fmeasure(klmean2_clust$clustering, dtw_clust$clustering)
  klmean3_fmeas_pam[[w]] <- fmeasure(klmean3_clust$clustering, dtw_clust$clustering)

}

df <- data.frame(unlist(min_fmeas_pam),unlist(apx_fmeas_pam),
                 unlist(klmean2_fmeas_pam), unlist(klmean3_fmeas_pam))#,unlist(cosine2_fmeas_pam), unlist(cosine3_fmeas_pam))
colnames(df) <- c("mindist","apxdist", "klavg2", "klavg3")#, "cosine2", "cosine3")
df <- data.frame(SequenceSize = seq(10,150,10), df)

df.m <- melt(df, id.vars = "SequenceSize")
ggplot(df.m) + geom_smooth(aes(x = SequenceSize, y = value, linetype = variable), size=.75, color="black", se=F) +
  ylab("F-Measure") + xlab("Sequence Size") + ggtitle(fname) +
  scale_linetype_manual(values=c("solid", "dashed", "dotdash", "dotted")) +
  theme_bw() +
  theme(panel.border = element_blank(), legend.key = element_blank(),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  legend.key.width = unit(2, "line"),
  axis.line = element_line(colour = "black")) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(10))
colMeans(df)
