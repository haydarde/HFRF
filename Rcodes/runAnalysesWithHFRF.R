  rm(list = ls(all.names = TRUE))
  rootPath <- "~/Documents/makaleler/Furkan_II/toGitHub/"
  setwd(rootPath)
  library(caret)
  library(kernlab)
  library(cluster)
  
  source(paste0(rootPath,"Rcodes/HFRF.R"))
  

  d <- 1
  contVars <- array() 
  data <- read.csv(paste0(rootPath,"/datasets/dataset", d,".csv"))
  
  if ( d == 1 ){
    contVars <- c(1:6)
    ranges <- apply(data,2,range)
    allRanges <- ranges[2,]-ranges[1,]
    allRanges[9] <- allRanges[9]+1
    load(file = paste0(rootPath,"distanceMatrices/genDist1.RData"))
    genDist <- genDist1
    c <- 2
    m <- 3
    p <- 1.1
    svmHypers <- matrix(c(0.02575229, 19.12452,
                          0.1575426, 20.16324),
                        nrow = 2, ncol =2, byrow = T)
    cat("Running Dataset", d,"\n","with m = ", m," and p = ", p, "\n")
  } else if ( d == 2 ){
    contVars <- c(1:5)
    ranges <- apply(data,2,range)
    allRanges <- ranges[2,]-ranges[1,]
    allRanges[8:10] <- allRanges[c(8:10)]+1
    load(file = paste0(rootPath,"distanceMatrices/genDist2.RData"))
    genDist <- genDist2
    c <- 2
    m <- 1.5
    p <- 1.1
    svmHypers <- matrix(c(0.297187, 3.289575,
                          0.04284169, 12.79675),
                        nrow = 2, ncol =2, byrow = T) 
    cat("Running Dataset", d,"\n","with m = ", m," and p = ", p, "\n")
  } else if( d == 3 ){
    contVars <- c(1:5)
    ranges <- apply(data,2,range)
    allRanges <- ranges[2,]-ranges[1,]
    load(file = paste0(rootPath,"distanceMatrices/genDist3.RData"))
    genDist <- genDist3
    c <- 5
    m <- 1.5
    p <- 1.1
    svmHypers <- matrix(c(0.1876616, 104.4838,
                          0.01801197, 164.3519,
                          0.1216698, 1.467394,
                          0.1329445, 3.166073,
                          0.006121353, 443.7304),
                        nrow = 5, ncol =2, byrow = T) 
    cat("Running Dataset", d,"\n","with m = ", m," and p = ", p, "\n")
  } else if( d == 4 ){
    contVars <- c(1:6)
    ranges <- apply(data,2,range)
    allRanges <- ranges[2,]-ranges[1,]
    allRanges[c(8,9)] <- allRanges[c(8,9)]+1
    load(file = paste0(rootPath,"distanceMatrices/genDist4.RData"))
    genDist <- genDist4
    c <- 3 
    m <- 1.5
    p <- 1.1
    svmHypers <- matrix(c(0.1676458, 6.194395,
                          0.04083332, 20.88579,
                          0.1541289, 16.73755),
                        nrow = 5, ncol =2, byrow = T)
    cat("Running Dataset", d,"\n","with m = ", m," and p = ", p, "\n")
  }  else if ( d == 5 ){
    contVars <- c(1:8)
    ranges <- apply(data,2,range)
    allRanges <- ranges[2,]-ranges[1,]
    allRanges[10:12] <- allRanges[10:12]+1
    load(file = paste0(rootPath,"distanceMatrices/genDist5.RData"))
    genDist <- genDist5
    c <- 6 
    m <- 3
    p <- 1.25
    svmHypers <- matrix(c(0.01424481, 736.5986,
                          0.01256321, 2.717128,
                          0.01570372, 382.1433,
                          0.0581623, 533.1814,
                          0.1173727, 8.38066,
                          0.1416243, 10.76749),
                        nrow = 6, ncol =2, byrow = T)
    cat("Running Dataset", d,"\n","with m = ", m," and p = ", p, "\n")
  } else if ( d == 6 ){
    contVars <- c(1:5)
    ranges <- apply(data,2,range)
    allRanges <- ranges[2,]-ranges[1,]
    allRanges[7] <- allRanges[7]+1
    load(file = paste0(rootPath,"distanceMatrices/genDist6.RData"))
    genDist <- genDist6
    c <- 3 
    m <- 3
    p <- 1.1
    svmHypers <- matrix(c( 0.02209247, 387.4695,
                           0.234909, 43.10924,
                           0.7577936, 84.65255),
                        nrow = 3, ncol =2, byrow = T)
    cat("Running Dataset", d,"\n","with m = ", m," and p = ", p, "\n")
  }
  
  set.seed(13579)
  data_part <- createDataPartition(y = data$lprice, p = 0.8, list = F)
  train <- data[data_part,]
  test <- data[-data_part,]
  
  train_x <- train[, -c(3)]
  train_y <- train[, 3]
  train_xy <- train
  test_x <- test[, -c(3)] 
  test_y <- test[, 3]
  test_xy <- test
  
  specResTest <- matrix(NA,20,4)
  HFRFfit1 <- HFRF(formula = as.formula(lprice ~ .), m = m, c = c, p = p, 
                       train = train, test = test, train_x = train_x, test_x = test_x, train_y = train_y, test_y = test_y,
                       train_control = train_control, contVars = contVars,
                       ranges = allRanges, genDist = genDist, svmHypers = svmHypers)
  specResTest[1,] <- HFRFfit1$GOF_test
  nClust <- HFRFfit1$nClust
  print(HFRFfit1$GOF_test)
  
  for (tt in 2:20){
    HFRFfit <- HFRF(formula = as.formula(lprice ~ .), m = m, c = nClust, p = p, train = train, test = test, train_x = train_x, test_x = test_x, train_y = train_y, test_y = test_y,
                        train_control = train_control, contVars = contVars,
                        ranges = allRanges, genDist = genDist, svmHypers = NULL)
    specResTest[tt,] <- HFRFfit$GOF_test
    print(HFRFfit$GOF_test)
  }
  