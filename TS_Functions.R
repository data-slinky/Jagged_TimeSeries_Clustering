##############################################################################################
# Loading the necessary packages
##############################################################################################

installIfNeeded = function(cliblist){
  libsNeeded = cliblist
  libsNeeded = libsNeeded[!(libsNeeded %in% installed.packages()[,"Package"])]
  if(length(libsNeeded)>0) install.packages(libsNeeded)
}

packages_needed <- c("fpc","plyr", "devtools", "data.table", "stringi","entropy", "rgl", "biogram","ggplot2", "stringdist",
                     "TSclust", "zoo", "ape", "reshape2", "seewave", "clusterCrit", "fpc", "doParallel", "foreach", "lsa", "jsd")

installIfNeeded(packages_needed)

for (i in packages_needed){
  if( !is.element(i, .packages(all.available = TRUE)) ) {
    install.packages(i)
  }
  library(i,character.only = TRUE)
}

registerDoParallel()

##############################################################################################
# FUNCTIONS
##############################################################################################

eucldist <- function(time_series) {
  size <- length(colnames(time_series))
  dist = data.table(matrix(0, nrow = size, ncol = size))
  colnames(dist) <- colnames(time_series) # Assign the column names to the distance matrix
  dist <- foreach(i = 1:size, .combine='rbind') %:%
    foreach(j = 1:size, .combine='c', .packages=c("TSclust","seewave")) %dopar% {
      diss.EUCL(time_series[,i], time_series[,j])
    }
  return (dist)
}

dtwdist <- function(time_series) {
  size <- length(colnames(time_series))
  dist = data.table(matrix(0, nrow = size, ncol = size))
  colnames(dist) <- colnames(time_series) # Assign the column names to the distance matrix
  dist <- foreach(i = 1:size, .combine='rbind') %:%
    foreach(j = 1:size, .combine='c', .packages=c("TSclust","seewave")) %dopar% {
      diss.DTWARP(time_series[,i], time_series[,j], plot=FALSE)
    }
  return (dist)
}

mindist <- function(time_series, w, alpha){
  size <- length(colnames(time_series))
  dist = df = data.table(matrix(0, nrow = size, ncol = size))
  colnames(dist) <- colnames(time_series) # Assign the column names to the distance matrix
  dist <- foreach(i = 1:size, .combine='rbind') %:%
    foreach(j = 1:size, .combine='c', .packages=c("TSclust","seewave")) %dopar% {
      diss.MINDIST.SAX(time_series[,i], time_series[,j], w, alpha, plot=FALSE)
  }
  return (dist)
}

apxdist <- function(time_series, w, alpha){
  size <- length(colnames(time_series))
  dist = data.table(matrix(0, nrow = size, ncol = size))
  colnames(dist) <- colnames(time_series) # Assign the column names to the distance matrix
  breakpoint = qnorm(0:w/w)
  breakpoint[w] = max(time_series)
  breakpoint[1] = min(time_series)
  sequences <- apply(data_ts,2,SAX,alpha,w,"gaussian",collapse="")
  dist <- foreach(i = 1:size, .combine='rbind') %:%
    foreach(j = 1:size, .combine='c', .packages=c("TSclust","seewave")) %dopar% {
      vi = match(unlist(strsplit(sequences[i], split="")),letters)
      vj = match(unlist(strsplit(sequences[j], split="")),letters)
      d = 0
      for (k in 1:nchar(sequences[i])){
        d = d + abs((breakpoint[vi[k]]+breakpoint[vi[k]+1]) - (breakpoint[vj[k]]+breakpoint[vj[k]+1]))/2
      }
      (sqrt(size/w)) * sqrt(d)
    }
  return (dist)
}

edit <- function(time_series, alpha){
  size <- length(colnames(time_series))
  dist = df = data.table(matrix(0, nrow = size, ncol = size))
  colnames(dist) <- colnames(time_series) # Assign the column names to the distance matrix
  sequences <- apply(data_ts,2,SAX,alpha,w,"gaussian",collapse="")
  dist <- foreach(i = 1:length(sequences), .combine='rbind') %:% 
    foreach(j = 1:length(sequences), .combine='c', .packages=c("stringdist","seewave")) %dopar% {
      if (sequences[i]==sequences[j]) {0}
      else {stringdist(sequences[i], sequences[j], "osa") + 0.25*((nchar(sequences[i]) + nchar(sequences[j])) + 2* min(table(strsplit(sequences[i], "")[[1]]), table(strsplit(sequences[j], "")[[1]])))}
    }
  return(dist)
}

qgram_dist <- function(time_series, alpha){
  size <- length(colnames(time_series))
  dist = df = data.table(matrix(0, nrow = size, ncol = size))
  colnames(dist) <- colnames(time_series) # Assign the column names to the distance matrix
  sequences <- apply(data_ts,2,SAX,alpha,w,"gaussian",collapse="")
  dist <- foreach(i = 1:length(sequences), .combine='rbind') %:% 
    foreach(j = 1:length(sequences), .combine='c', .packages=c("stringdist","seewave")) %dopar% {
      stringdist(sequences[i], sequences[j], method="hamming")
    }
  return(dist)
}

# Count of the L-tuples
Lcount <- function(string, L, alpha){
  alphabet <- unlist( strsplit("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ", "") )
  x1 <- strsplit(string, '')[[1]] # Split the string into a vector values
  trigrams <- count_ngrams(x1, L, alphabet[1:alpha]) # Get the ngrams frequencies
  trigrams <- as.matrix(trigrams) # Convert triplet matrix to regular matrix
  colnames(trigrams) <- gsub("[^a-zA-Z]", "", colnames(trigrams)) # Clean up the column names
  return(trigrams)
}

# Compute the cosine distance with SVD
cosinedist <- function(time_series, L, alpha, dim=20) {
  size <- length(colnames(time_series))
  dist = data.table(matrix(0, nrow = size, ncol = size))
  colnames(dist) <- colnames(time_series) # Assign the column names to the distance matrix
  sequences <- apply(data_ts,2,SAX,alpha,w,"gaussian",collapse="")
  dist <- foreach (i = 1:size) %:%
    foreach (j = 1:size) %dopar% {
      1 - sum(count1*count2)/sqrt(sum(count1^2)*sum(count2^2))
    }
  return(dist)
}

# Compute the cosine distance with SVD
cosine_test <- function(time_series, L, alpha, dim=20) {
  symbols = list()
  sequences <- apply(time_series,2,SAX, alpha, w, "gaussian",collapse="")
  symbols <- lapply(sequences, Lcount, L=L, alpha=alpha)
  return(symbols)
}

#count_matrix <- as.data.table(cosine_test(data_ts, 5, 5))

reduce <- function(A,dim) {
  #Calculates the SVD
  sing <- svd(A)
  #Approximate each result of SVD with the given dimension
  u<-as.matrix(sing$u[, 1:dim])
  v<-as.matrix(sing$v[, 1:dim])
  d<-as.matrix(diag(sing$d)[1:dim, 1:dim])
  return(u%*%d%*%t(v))
}


# With the scaling
freq_distr <- function(string, L, alpha){
  alphabet <- unlist( strsplit("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ", "") )
  x1 <- strsplit(string, '')[[1]] # Split the string into a vector values
  trigrams <- count_ngrams(x1, L, alphabet[1:alpha]) # Get the ngrams frequencies
  trigrams <- as.matrix(trigrams) # Convert triplet matrix to regular matrix
  trigrams <- trigrams + 0.5*length(alphabet)^-L # Add the scaling counts to avoid KL = infinity
  colnames(trigrams) <- gsub("[^a-zA-Z]", "", colnames(trigrams)) # Clean up the column names
  distribution <- trigrams/((nchar(string)-L+1)+0.5) # Compute the distribution
  return(distribution)
}

kldiv <- function(time_series, L, alpha) {
  size <- length(colnames(time_series))
  dist <- data.table(matrix(0, nrow = size, ncol = size))
  colnames(dist) <- colnames(time_series) # Assign the column names to the distance matrix
  sequences <- apply(data_ts,2,SAX,alpha,w,"gaussian",collapse="")
  dist <- foreach(i = 1:length(sequences), .combine='rbind') %:% 
    foreach(j = 1:length(sequences),.combine='c', .export=c("freq_distr"), .packages=c("biogram","seewave")) %dopar% {
      kl.dist(freq_distr(sequences[i],L,alpha),freq_distr(sequences[j],L,alpha))$D
  }
  return(dist)
}

ksdist <- function(time_series, L, alpha) {
  size <- length(colnames(time_series))
  dist <- data.table(matrix(0, nrow = size, ncol = size))
  colnames(dist) <- colnames(time_series) # Assign the column names to the distance matrix
  sequences <- apply(data_ts,2,SAX,alpha,w,"gaussian",collapse="")
  dist <- foreach(i = 1:length(sequences), .combine='rbind') %:% 
    foreach(j = 1:length(sequences),.combine='c', .export=c("freq_distr"), .packages=c("biogram","seewave")) %dopar% {
      ks.dist(freq_distr(sequences[i],L,alpha),freq_distr(sequences[j],L,alpha))$D
    }
  return(dist)
}

jsdiv <- function(time_series, L, alpha) {
  size <- length(colnames(time_series))
  dist <- data.table(matrix(0, nrow = size, ncol = size))
  colnames(dist) <- colnames(time_series) # Assign the column names to the distance matrix
  sequences <- apply(data_ts,2,SAX,alpha,w,"gaussian",collapse="")
  dist <- foreach(i = 1:length(sequences), .combine='rbind') %:% 
    foreach(j = 1:length(sequences),.combine='c', .export=c("freq_distr"), .packages=c("biogram","jsd")) %dopar% {
      JSD(freq_distr(sequences[i],L,alpha),freq_distr(sequences[j],L,alpha))$D
    }
  return(dist)
}

fmeasure <- function(predict, actual){
  precision <- as.numeric(extCriteria(predict, actual, "Precision"))
  recall <- as.numeric(extCriteria(predict, actual, "Recall"))
  fmeasure <- 2 * precision * recall / (precision + recall)
 return(fmeasure)
}

ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

list_of_SAX <- function(time_series, L, alpha) {
  symbols = list()
  size <- length(colnames(time_series))
  for (i in 1:size){
    symbols[[i]] <- SAX(time_series[,i], alphabet=alpha, PAA=w, breakpoints= "gaussian", collapse=" ")
  }
  return(symbols)
}

euclidean <- function(time_series, L, alpha) {
  size <- length(colnames(time_series))
  dist = data.table(matrix(0, nrow = size, ncol = size))
  colnames(dist) <- colnames(time_series) # Assign the column names to the distance matrix
  for (i in 1:size){
    si = SAX(time_series[,i], alphabet=alpha, PAA=w, breakpoints= "gaussian", collapse="")
    count1 = Lcount(si,L,alpha)
    for (j in 1:size){
      sj = SAX(time_series[,j], alphabet=alpha, PAA=w, breakpoints= "gaussian", collapse="")
      count2 = Lcount(sj,L,alpha)
      dist[i,j] = sqrt(sum((count1 - count2) ^ 2))
    }
  }
  return(dist)
}

number_ticks <- function(n) {function(limits) pretty(limits, n)}