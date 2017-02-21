source("~/Dropbox/Research/Time-Series/Code/TS_Functions.R")
#source("C:/Users/XeroXtancy/Dropbox/Research/Time-Series/Code/TS_Functions.R")

##############################################################################################
# LOAD IN DATA
##############################################################################################

fname = "ToeSegmentation2"

path = paste("/Users/johnkimnguyen/Dropbox/Research/Time-Series/Code/Varied_Alpha/", fname,"/", sep="")
#path = paste("C:/Users/XeroXtancy/Dropbox/Research/Time-Series/Code/UCR_To_Do_w40_a5/", fname,"/", sep="")

setwd(path)

data <- read.table(paste(fname,'_TRAIN', sep=""), sep=',')
data <- as.data.frame(t(data)) # Transpose the data since # seq is num of rows
data <-  data[-1,] # Drop the first row since we dont need class label
data <- data.frame(TIME = 1:nrow(data), data)
data_ts <- read.zoo(data) # Convert data to zoo time-series

DTWDIST = read.table("DTWDIST.csv", sep=",", header=TRUE)[-1]
w=20

##############################################################################################
# COMPUTE F-MEAUSRES WITH VARYING ALPHABET SIZE
##############################################################################################

min_fmeas_pam <- list()
edit_fmeas_pam <- list()
apx_fmeas_pam <- list()
# eucl2_fmeas_pam <- list()
# eucl3_fmeas_pam <- list()
klmean2_fmeas_pam <- list()
klmean3_fmeas_pam <- list()

k = 2

for (alpha in 2:20) {

# fname <- list_of_SAX(data_ts, L, alpha)
# lapply(fname, write, paste("list_of_SAX_a", alpha, "_w20", ".txt", sep=""), append=TRUE, ncolumns=1000)

#   print(alpha)
#   MINDIST <- mindist(data_ts, w, alpha)
#   EDITDIST <- edit(data_ts, alpha)
#   APXDIST <- apxdist(data_ts, w, alpha)
#   write.csv(MINDIST, file=paste('MINDIST', alpha, ".csv", sep=""))
#   write.csv(APXDIST, file=paste('APXDIST', alpha, ".csv", sep=""))
#   write.csv(EDITDIST, file=paste('EDITDIST', alpha, ".csv", sep=""))
# }

DTWDIST = read.table("DTWDIST.csv", sep=",", header=TRUE)[-1]
APXDIST = read.table(paste('APXDIST', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
MINDIST = read.table(paste('MINDIST', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
EDITDIST = read.table(paste('EDITDIST', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#EUCL2 = read.table(paste('EUCL-2a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#EUCL3 = read.table(paste('EUCL-2a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
KLMEAN2 = read.table(paste('KLMEAN-2a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
KLMEAN3 = read.table(paste('KLMEAN-3a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]

  dtw_clust <- pam(DTWDIST, k)
  min_clust <- pam(MINDIST, k)
  edit_clust <- pam(EDITDIST, k)
  apx_clust <- pam(APXDIST, k)
 #eucl2_clust <- pam(EUCL2, 3)
  #eucl3_clust <- pam(EUCL3, 3)
  klmean2_clust <- pam(KLMEAN2, k)
  klmean3_clust <- pam(KLMEAN3, k)

#Get the F-measures
min_fmeas_pam[[alpha]] <- fmeasure(min_clust$clustering, dtw_clust$clustering)
edit_fmeas_pam[[alpha]] <- fmeasure(edit_clust$clustering, dtw_clust$clustering)
apx_fmeas_pam[[alpha]] <- fmeasure(apx_clust$clustering, dtw_clust$clustering)
#eucl2_fmeas_pam[[alpha]] <- fmeasure(eucl2_clust$clustering, dtw_clust$clustering)
#eucl3_fmeas_pam[[alpha]] <- fmeasure(eucl3_clust$clustering, dtw_clust$clustering)
klmean2_fmeas_pam[[alpha]] <- fmeasure(klmean2_clust$clustering, dtw_clust$clustering)
klmean3_fmeas_pam[[alpha]] <- fmeasure(klmean3_clust$clustering, dtw_clust$clustering)
}

df <- data.frame(unlist(min_fmeas_pam),unlist(apx_fmeas_pam),unlist(edit_fmeas_pam),
                 #unlist(eucl2_fmeas_pam), unlist(eucl3_fmeas_pam),
                 unlist(klmean2_fmeas_pam), unlist(klmean3_fmeas_pam))#,unlist(cosine2_fmeas_pam), unlist(cosine3_fmeas_pam))
colnames(df) <- c("mindist","apxdist", "eed", "klavg2", "klavg3")#, "cosine2", "cosine3")
df<- data.frame(AlphabetSize = 1:nrow(df), df)
df['AlphabetSize'] = df['AlphabetSize'] + 1
# write.csv(t(df), file='varied_alpha_fmeas_RESULT.csv')

df.m <- melt(df, id.vars = "AlphabetSize")
df.m <- subset(df.m, variable != "eed" & variable != "eucl2" & variable != "eucl3")
fmeas_plot <- ggplot(df.m) + 
  #geom_point(aes(x = AlphabetSize, y = value, colour = variable)) +
  geom_smooth(aes(x = AlphabetSize, y = value, linetype = variable), color="black", size=.75, se = F) +
  ylab("F-Measure") + xlab("Alphabet Size") + ggtitle(fname) + 
  theme_bw() + 
  scale_linetype_manual(values=c("solid", "dashed", "dotdash", "dotted"))+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.width = unit(2, "line"),
        legend.key = element_blank(),
        legend.title=element_blank()) +
  scale_x_continuous(breaks=number_ticks(10)) +
  scale_y_continuous(breaks=number_ticks(10)) 
fmeas_plot
colMeans(df)



##############################################################################################
# COMPUTE F-RECALL WITH VARYING ALPHABET SIZE
##############################################################################################
# 
# min_prec_pam <- list()
# edit_prec_pam <- list()
# apx_prec_pam <- list()
# eucl2_prec_pam <- list()
# eucl3_prec_pam <- list()
# klmean2_prec_pam <- list()
# klmean3_prec_pam <- list()
# 
# for (alpha in 2:26) {
# 
#   
#   DTWDIST = read.table("DTWDIST.csv", sep=",", header=TRUE)[-1]
#   APXDIST = read.table(paste('APXDIST', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   MINDIST = read.table(paste('MINDIST', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   EDITDIST = read.table(paste('EDITDIST', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   EUCL2 = read.table(paste('EUCL-2a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   EUCL3 = read.table(paste('EUCL-2a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   KLMEAN2 = read.table(paste('KLMEAN-2a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   KLMEAN3 = read.table(paste('KLMEAN-3a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   
#   dtw_clust <- pam(DTWDIST, 3)
#   min_clust <- pam(MINDIST, 3)
#   edit_clust <- pam(EDITDIST, 3)
#   apx_clust <- pam(APXDIST, 3)
#   eucl2_clust <- pam(EUCL2, 3)
#   eucl3_clust <- pam(EUCL3, 3)
#   klmean2_clust <- pam(KLMEAN2, 3)
#   klmean3_clust <- pam(KLMEAN3, 3)
#   
#   #Get the F-measures
#   min_prec_pam[[alpha]] <- extCriteria(min_clust$clustering, dtw_clust$clustering, "Recall")
#   edit_prec_pam[[alpha]] <- extCriteria(edit_clust$clustering, dtw_clust$clustering, "Recall")
#   apx_prec_pam[[alpha]] <- extCriteria(apx_clust$clustering, dtw_clust$clustering, "Recall")
#   eucl2_prec_pam[[alpha]] <- extCriteria(eucl2_clust$clustering, dtw_clust$clustering, "Recall")
#   eucl3_prec_pam[[alpha]] <- extCriteria(eucl3_clust$clustering, dtw_clust$clustering, "Recall")
#   klmean2_prec_pam[[alpha]] <- extCriteria(klmean2_clust$clustering, dtw_clust$clustering, "Recall")
#   klmean3_prec_pam[[alpha]] <- extCriteria(klmean3_clust$clustering, dtw_clust$clustering, "Recall")
# }
# 
# df_recall <- data.frame(unlist(min_prec_pam),unlist(apx_prec_pam),unlist(edit_prec_pam),
#                       unlist(eucl2_prec_pam), unlist(eucl3_prec_pam),
#                       unlist(klmean2_prec_pam), unlist(klmean3_prec_pam))#,unlist(cosine2_fmeas_pam), unlist(cosine3_fmeas_pam))
# colnames(df_recall) <- c("mindist","apxdist", "eed", "eucl2","eucl3", "klmean2", "klmean3")#, "cosine2", "cosine3")
# df_recall<- data.frame(AlphabetSize = 1:nrow(df_recall), df_recall)
# df_recall['AlphabetSize'] = df_recall['AlphabetSize'] + 1
# 
# df_recall.m <- melt(df_recall, id.vars = "AlphabetSize")
# df_recall.m <- subset(df_recall.m, variable != "eed" & variable != "eucl2" & variable != "eucl3")
# ggplot(df_recall.m) + geom_line(aes(x = AlphabetSize, y = value, colour = variable)) + ylab("Recall") + xlab("Alphabet Size") + ggtitle("CinC_ECG_torso (w = 50)") + 
# colMeans(df_recall)

# write.csv(t(df_recall), file='varied_alpha_recall_RESULT.csv')


##############################################################################################
# COMPUTE F-RECALL WITH VARYING ALPHABET SIZE
##############################################################################################

# min_prec_pam <- list()
# edit_prec_pam <- list()
# apx_prec_pam <- list()
# eucl2_prec_pam <- list()
# eucl3_prec_pam <- list()
# klmean2_prec_pam <- list()
# klmean3_prec_pam <- list()
# 
# for (alpha in 2:26) {
#   
#   
#   DTWDIST = read.table("DTWDIST.csv", sep=",", header=TRUE)[-1]
#   APXDIST = read.table(paste('APXDIST', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   MINDIST = read.table(paste('MINDIST', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   EDITDIST = read.table(paste('EDITDIST', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   EUCL2 = read.table(paste('EUCL-2a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   EUCL3 = read.table(paste('EUCL-2a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   KLMEAN2 = read.table(paste('KLMEAN-2a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   KLMEAN3 = read.table(paste('KLMEAN-3a', alpha, ".csv", sep=""), sep=",", header=TRUE)[-1]
#   
#   dtw_clust <- pam(DTWDIST, 3)
#   min_clust <- pam(MINDIST, 3)
#   edit_clust <- pam(EDITDIST, 3)
#   apx_clust <- pam(APXDIST, 3)
#   eucl2_clust <- pam(EUCL2, 3)
#   eucl3_clust <- pam(EUCL3, 3)
#   klmean2_clust <- pam(KLMEAN2, 3)
#   klmean3_clust <- pam(KLMEAN3, 3)
#   
#   #Get the F-measures
#   min_prec_pam[[alpha]] <- extCriteria(min_clust$clustering, dtw_clust$clustering, "Precision")
#   edit_prec_pam[[alpha]] <- extCriteria(edit_clust$clustering, dtw_clust$clustering, "Precision")
#   apx_prec_pam[[alpha]] <- extCriteria(apx_clust$clustering, dtw_clust$clustering, "Precision")
#   eucl2_prec_pam[[alpha]] <- extCriteria(eucl2_clust$clustering, dtw_clust$clustering, "Precision")
#   eucl3_prec_pam[[alpha]] <- extCriteria(eucl3_clust$clustering, dtw_clust$clustering, "Precision")
#   klmean2_prec_pam[[alpha]] <- extCriteria(klmean2_clust$clustering, dtw_clust$clustering, "Precision")
#   klmean3_prec_pam[[alpha]] <- extCriteria(klmean3_clust$clustering, dtw_clust$clustering, "Precision")
# }
# 
# df_prec <- data.frame(unlist(min_prec_pam),unlist(apx_prec_pam),unlist(edit_prec_pam),
#                  unlist(eucl2_prec_pam), unlist(eucl3_prec_pam),
#                  unlist(klmean2_prec_pam), unlist(klmean3_prec_pam))#,unlist(cosine2_fmeas_pam), unlist(cosine3_fmeas_pam))
# colnames(df_prec) <- c("mindist","apxdist", "eed", "eucl2","eucl3", "klmean2", "klmean3")#, "cosine2", "cosine3")
# df_prec<- data.frame(AlphabetSize = 1:nrow(df_prec), df_prec)
# df_prec['AlphabetSize'] = df_prec['AlphabetSize'] + 1
# 
# df_prec.m <- melt(df_prec, id.vars = "AlphabetSize")
# df_prec.m <- subset(df_prec.m, variable != "eed" & variable != "eucl2" & variable != "eucl3")
# ggplot(df_prec.m) + geom_line(aes(x = AlphabetSize, y = value, colour = variable)) + ylab("Precision") + xlab("Alphabet Size") + ggtitle("CinC_ECG_torso (w = 50)")
# colMeans(df_prec)
# 
# write.csv(t(df_prec), file='varied_alpha_precision_RESULT.csv')
