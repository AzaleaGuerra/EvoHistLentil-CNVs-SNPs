
dp <- read.delim("../data/xhmm/LensXhmm.filtered_centered.RD.txt")
dp <- na.omit(dp)

meta <- read.delim("../meta/LDP_meta.txt")

#keep only the depth of CNV regions present in more than 20% of the samples
dp.f <- numeric(0)
for (i in 2:ncol(dp)) {
  if (length(dp[dp[,i]>80 | dp[,i]< -80, i]) > 5) {
    dp.f <- cbind(dp.f, dp[,i])
  }
}
dp.f <- as.data.frame(dp.f)
dp.f <- cbind(Matrix=dp$Matrix, dp.f)

dp.meta <- numeric(0)
for (i in 1:nrow(meta)) {
  a <- dp.f[dp.f$Matrix == as.character(meta[i,1]),]
  b <- cbind(a, meta[i,c(8:9)])
  dp.meta <- rbind(dp.meta, b)
}

#Vst for Pop 1 and 2
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==1 | dp.meta$GenGroup==2, i])
  b <- (var(dp.meta[dp.meta$GenGroup==1, i])*length(dp.meta[dp.meta$GenGroup==1, i]) + 
        var(dp.meta[dp.meta$GenGroup==2, i])*length(dp.meta[dp.meta$GenGroup==2, i])) / 
        (length(dp.meta[dp.meta$GenGroup==1, i]) + length(dp.meta[dp.meta$GenGroup==2, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v12 <- mean(v$V2, na.rm = T)

#Vst Pop 1 and 3
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==1 | dp.meta$GenGroup==3, i])
  b <- (var(dp.meta[dp.meta$GenGroup==1, i])*length(dp.meta[dp.meta$GenGroup==1, i]) + 
          var(dp.meta[dp.meta$GenGroup==3, i])*length(dp.meta[dp.meta$GenGroup==3, i])) / 
    (length(dp.meta[dp.meta$GenGroup==1, i]) + length(dp.meta[dp.meta$GenGroup==3, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v13 <- mean(v$V2, na.rm = T)

#Vst Pop 1 and 4
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==1 | dp.meta$GenGroup==4, i])
  b <- (var(dp.meta[dp.meta$GenGroup==1, i])*length(dp.meta[dp.meta$GenGroup==1, i]) + 
          var(dp.meta[dp.meta$GenGroup==4, i])*length(dp.meta[dp.meta$GenGroup==4, i])) / 
    (length(dp.meta[dp.meta$GenGroup==1, i]) + length(dp.meta[dp.meta$GenGroup==4, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v14 <- mean(v$V2, na.rm = TRUE)

#Vst Pop 1 and 5
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==1 | dp.meta$GenGroup==5, i])
  b <- (var(dp.meta[dp.meta$GenGroup==1, i])*length(dp.meta[dp.meta$GenGroup==1, i]) + 
          var(dp.meta[dp.meta$GenGroup==5, i])*length(dp.meta[dp.meta$GenGroup==5, i])) / 
    (length(dp.meta[dp.meta$GenGroup==1, i]) + length(dp.meta[dp.meta$GenGroup==5, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v15 <- mean(v$V2, na.rm = T)

#Vst Pop 1 and 6
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==1 | dp.meta$GenGroup==6, i])
  b <- (var(dp.meta[dp.meta$GenGroup==1, i])*length(dp.meta[dp.meta$GenGroup==1, i]) + 
          var(dp.meta[dp.meta$GenGroup==6, i])*length(dp.meta[dp.meta$GenGroup==6, i])) / 
    (length(dp.meta[dp.meta$GenGroup==1, i]) + length(dp.meta[dp.meta$GenGroup==6, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v16 <- mean(v$V2, na.rm = T)

#Vst Pop 1 and 7
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==1 | dp.meta$GenGroup==7, i])
  b <- (var(dp.meta[dp.meta$GenGroup==1, i])*length(dp.meta[dp.meta$GenGroup==1, i]) + 
          var(dp.meta[dp.meta$GenGroup==7, i])*length(dp.meta[dp.meta$GenGroup==7, i])) / 
    (length(dp.meta[dp.meta$GenGroup==1, i]) + length(dp.meta[dp.meta$GenGroup==7, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v17 <- mean(v$V2, na.rm = T)

#Vst Pop 1 and 8
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==1 | dp.meta$GenGroup==8, i])
  b <- (var(dp.meta[dp.meta$GenGroup==1, i])*length(dp.meta[dp.meta$GenGroup==1, i]) + 
          var(dp.meta[dp.meta$GenGroup==8, i])*length(dp.meta[dp.meta$GenGroup==8, i])) / 
    (length(dp.meta[dp.meta$GenGroup==1, i]) + length(dp.meta[dp.meta$GenGroup==8, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v18 <- mean(v$V2, na.rm = T)

#Vst Pop 2 and 3
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==2 | dp.meta$GenGroup==3, i])
  b <- (var(dp.meta[dp.meta$GenGroup==2, i])*length(dp.meta[dp.meta$GenGroup==2, i]) + 
          var(dp.meta[dp.meta$GenGroup==3, i])*length(dp.meta[dp.meta$GenGroup==3, i])) / 
    (length(dp.meta[dp.meta$GenGroup==2, i]) + length(dp.meta[dp.meta$GenGroup==3, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v23 <- mean(v$V2, na.rm = T)

#Vst Pop 2 and 4
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==2 | dp.meta$GenGroup==4, i])
  b <- (var(dp.meta[dp.meta$GenGroup==2, i])*length(dp.meta[dp.meta$GenGroup==2, i]) + 
          var(dp.meta[dp.meta$GenGroup==4, i])*length(dp.meta[dp.meta$GenGroup==4, i])) / 
    (length(dp.meta[dp.meta$GenGroup==2, i]) + length(dp.meta[dp.meta$GenGroup==4, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v24 <- mean(v$V2, na.rm = T)

#Vst Pop 2 and 5
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==2 | dp.meta$GenGroup==5, i])
  b <- (var(dp.meta[dp.meta$GenGroup==2, i])*length(dp.meta[dp.meta$GenGroup==2, i]) + 
          var(dp.meta[dp.meta$GenGroup==5, i])*length(dp.meta[dp.meta$GenGroup==5, i])) / 
    (length(dp.meta[dp.meta$GenGroup==2, i]) + length(dp.meta[dp.meta$GenGroup==5, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v25 <- mean(v$V2, na.rm = T)

#Vst Pop 2 and 6
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==2 | dp.meta$GenGroup==6, i])
  b <- (var(dp.meta[dp.meta$GenGroup==2, i])*length(dp.meta[dp.meta$GenGroup==2, i]) + 
          var(dp.meta[dp.meta$GenGroup==6, i])*length(dp.meta[dp.meta$GenGroup==6, i])) / 
    (length(dp.meta[dp.meta$GenGroup==2, i]) + length(dp.meta[dp.meta$GenGroup==6, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v26 <- mean(v$V2, na.rm = T)

#Vst Pop 2 and 7
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==2 | dp.meta$GenGroup==7, i])
  b <- (var(dp.meta[dp.meta$GenGroup==2, i])*length(dp.meta[dp.meta$GenGroup==2, i]) + 
          var(dp.meta[dp.meta$GenGroup==7, i])*length(dp.meta[dp.meta$GenGroup==7, i])) / 
    (length(dp.meta[dp.meta$GenGroup==2, i]) + length(dp.meta[dp.meta$GenGroup==7, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v27 <- mean(v$V2, na.rm = T)

#Vst Pop 2 and 8
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==2 | dp.meta$GenGroup==8, i])
  b <- (var(dp.meta[dp.meta$GenGroup==2, i])*length(dp.meta[dp.meta$GenGroup==2, i]) + 
          var(dp.meta[dp.meta$GenGroup==8, i])*length(dp.meta[dp.meta$GenGroup==8, i])) / 
    (length(dp.meta[dp.meta$GenGroup==2, i]) + length(dp.meta[dp.meta$GenGroup==8, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v28 <- mean(v$V2, na.rm = T)

#Vst Pop 3 and 4
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==3 | dp.meta$GenGroup==4, i])
  b <- (var(dp.meta[dp.meta$GenGroup==3, i])*length(dp.meta[dp.meta$GenGroup==3, i]) + 
          var(dp.meta[dp.meta$GenGroup==4, i])*length(dp.meta[dp.meta$GenGroup==4, i])) / 
    (length(dp.meta[dp.meta$GenGroup==3, i]) + length(dp.meta[dp.meta$GenGroup==4, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v34 <- mean(v$V2, na.rm = T)

#Vst Pop 3 and 5
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==3 | dp.meta$GenGroup==5, i])
  b <- (var(dp.meta[dp.meta$GenGroup==3, i])*length(dp.meta[dp.meta$GenGroup==3, i]) + 
          var(dp.meta[dp.meta$GenGroup==5, i])*length(dp.meta[dp.meta$GenGroup==5, i])) / 
    (length(dp.meta[dp.meta$GenGroup==3, i]) + length(dp.meta[dp.meta$GenGroup==5, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v35 <- mean(v$V2, na.rm = T)

#Vst Pop 3 and 6
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==3 | dp.meta$GenGroup==6, i])
  b <- (var(dp.meta[dp.meta$GenGroup==3, i])*length(dp.meta[dp.meta$GenGroup==3, i]) + 
          var(dp.meta[dp.meta$GenGroup==6, i])*length(dp.meta[dp.meta$GenGroup==6, i])) / 
    (length(dp.meta[dp.meta$GenGroup==3, i]) + length(dp.meta[dp.meta$GenGroup==6, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v36 <- mean(v$V2, na.rm = T)

#Vst Pop 3 and 7
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==3 | dp.meta$GenGroup==7, i])
  b <- (var(dp.meta[dp.meta$GenGroup==3, i])*length(dp.meta[dp.meta$GenGroup==3, i]) + 
          var(dp.meta[dp.meta$GenGroup==7, i])*length(dp.meta[dp.meta$GenGroup==7, i])) / 
    (length(dp.meta[dp.meta$GenGroup==3, i]) + length(dp.meta[dp.meta$GenGroup==7, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v37 <- mean(v$V2, na.rm = T)

#Vst Pop 3 and 8
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==3 | dp.meta$GenGroup==8, i])
  b <- (var(dp.meta[dp.meta$GenGroup==3, i])*length(dp.meta[dp.meta$GenGroup==3, i]) + 
          var(dp.meta[dp.meta$GenGroup==8, i])*length(dp.meta[dp.meta$GenGroup==8, i])) / 
    (length(dp.meta[dp.meta$GenGroup==3, i]) + length(dp.meta[dp.meta$GenGroup==8, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v38 <- mean(v$V2, na.rm = T)

#Vst Pop 4 and 5
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==4 | dp.meta$GenGroup==5, i])
  b <- (var(dp.meta[dp.meta$GenGroup==4, i])*length(dp.meta[dp.meta$GenGroup==4, i]) + 
          var(dp.meta[dp.meta$GenGroup==5, i])*length(dp.meta[dp.meta$GenGroup==5, i])) / 
    (length(dp.meta[dp.meta$GenGroup==4, i]) + length(dp.meta[dp.meta$GenGroup==5, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v45 <- mean(v$V2, na.rm = T)

#Vst Pop 4 and 5
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==4 | dp.meta$GenGroup==6, i])
  b <- (var(dp.meta[dp.meta$GenGroup==4, i])*length(dp.meta[dp.meta$GenGroup==4, i]) + 
          var(dp.meta[dp.meta$GenGroup==6, i])*length(dp.meta[dp.meta$GenGroup==6, i])) / 
    (length(dp.meta[dp.meta$GenGroup==4, i]) + length(dp.meta[dp.meta$GenGroup==6, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v46 <- mean(v$V2, na.rm = T)

#Vst Pop 4 and 7
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==4 | dp.meta$GenGroup==7, i])
  b <- (var(dp.meta[dp.meta$GenGroup==4, i])*length(dp.meta[dp.meta$GenGroup==4, i]) + 
          var(dp.meta[dp.meta$GenGroup==7, i])*length(dp.meta[dp.meta$GenGroup==7, i])) / 
    (length(dp.meta[dp.meta$GenGroup==4, i]) + length(dp.meta[dp.meta$GenGroup==7, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v47 <- mean(v$V2, na.rm = T)

#Vst Pop 4 and 8
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==4 | dp.meta$GenGroup==8, i])
  b <- (var(dp.meta[dp.meta$GenGroup==4, i])*length(dp.meta[dp.meta$GenGroup==4, i]) + 
          var(dp.meta[dp.meta$GenGroup==8, i])*length(dp.meta[dp.meta$GenGroup==8, i])) / 
    (length(dp.meta[dp.meta$GenGroup==4, i]) + length(dp.meta[dp.meta$GenGroup==8, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v48 <- mean(v$V2, na.rm = T)

#Vst Pop 5 and 6
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==5 | dp.meta$GenGroup==6, i])
  b <- (var(dp.meta[dp.meta$GenGroup==5, i])*length(dp.meta[dp.meta$GenGroup==5, i]) + 
          var(dp.meta[dp.meta$GenGroup==6, i])*length(dp.meta[dp.meta$GenGroup==6, i])) / 
    (length(dp.meta[dp.meta$GenGroup==5, i]) + length(dp.meta[dp.meta$GenGroup==6, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v56 <- mean(v$V2, na.rm = T)


#Vst Pop 5 and 7
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==5 | dp.meta$GenGroup==7, i])
  b <- (var(dp.meta[dp.meta$GenGroup==5, i])*length(dp.meta[dp.meta$GenGroup==5, i]) + 
          var(dp.meta[dp.meta$GenGroup==7, i])*length(dp.meta[dp.meta$GenGroup==7, i])) / 
    (length(dp.meta[dp.meta$GenGroup==5, i]) + length(dp.meta[dp.meta$GenGroup==7, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v57 <- mean(v$V2, na.rm = T)


#Vst Pop 5 and 8
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==5 | dp.meta$GenGroup==8, i])
  b <- (var(dp.meta[dp.meta$GenGroup==5, i])*length(dp.meta[dp.meta$GenGroup==5, i]) + 
          var(dp.meta[dp.meta$GenGroup==8, i])*length(dp.meta[dp.meta$GenGroup==8, i])) / 
    (length(dp.meta[dp.meta$GenGroup==5, i]) + length(dp.meta[dp.meta$GenGroup==8, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v58 <- mean(v$V2, na.rm = T)

#Vst Pop 6 and 7
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==6 | dp.meta$GenGroup==7, i])
  b <- (var(dp.meta[dp.meta$GenGroup==6, i])*length(dp.meta[dp.meta$GenGroup==6, i]) + 
          var(dp.meta[dp.meta$GenGroup==7, i])*length(dp.meta[dp.meta$GenGroup==7, i])) / 
    (length(dp.meta[dp.meta$GenGroup==6, i]) + length(dp.meta[dp.meta$GenGroup==7, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v67 <- mean(v$V2, na.rm = T)

#Vst Pop 6 and 8
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==6 | dp.meta$GenGroup==8, i])
  b <- (var(dp.meta[dp.meta$GenGroup==6, i])*length(dp.meta[dp.meta$GenGroup==6, i]) + 
          var(dp.meta[dp.meta$GenGroup==8, i])*length(dp.meta[dp.meta$GenGroup==8, i])) / 
    (length(dp.meta[dp.meta$GenGroup==6, i]) + length(dp.meta[dp.meta$GenGroup==8, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v68 <- mean(v$V2, na.rm = T)

#Vst Pop 7 and 8
v <- numeric(0)
for (i in 2:(ncol(dp.meta)-2)) {
  a <- var(dp.meta[dp.meta$GenGroup==7 | dp.meta$GenGroup==8, i])
  b <- (var(dp.meta[dp.meta$GenGroup==7, i])*length(dp.meta[dp.meta$GenGroup==7, i]) + 
          var(dp.meta[dp.meta$GenGroup==8, i])*length(dp.meta[dp.meta$GenGroup==8, i])) / 
    (length(dp.meta[dp.meta$GenGroup==7, i]) + length(dp.meta[dp.meta$GenGroup==8, i]))
  v <- rbind(v, cbind(colnames(dp.meta)[i], (a-b)/a))
}
v <- as.data.frame(v)
v$V2 <- as.numeric(as.character(v$V2))
v$V2 <- replace(v$V2, v$V2 < 0, 0)
v78 <- mean(v$V2, na.rm = T)

