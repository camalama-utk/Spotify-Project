HitPredictionTotal = read.csv("Case5-HitPrediction.csv")

HitPredictionTotal = Case5_HitPrediction

library(tidyverse)
library(dplyr)
library(multcompView)
library(stringr)
library(regclass)

replace_rare_levels <- function(x, threshold = 20, newname = "Other")  {
  x <- factor(x)
  rare.levels <- names(which(sort(table(x)) <= threshold))
  if (length(rare.levels) == 0) { return(x) }
  levels(x)[which(levels(x) %in% rare.levels)] <- newname
  ST <- sort(table(x))
  if (ST[newname] <= threshold) {
    levels.to.combine <- which(levels(x) %in% c(newname, 
                                                names(ST)[2]))
    levels(x)[levels.to.combine] <- newname
    rare.levels <- c(rare.levels, names(ST)[2])
  }
  return(x)
}

examine_driver_Ycat <- function(formula,data,sort=TRUE,inside=TRUE,equal=TRUE) { 
  require(regclass)
  require(multcompView)
  require(stringr)
  
  FORM <- as.character(formula)
  temp <- str_extract(FORM,"^.*==")
  temp <- temp[!is.na(temp)]
  y.label <- gsub(" ==","",temp)
  temp <- str_extract(FORM,'\\"[:alnum:]*')
  temp <- temp[!is.na(temp)]
  level <- substr(temp,2,nchar(temp))
  FORM <- as.formula(formula)
  
  variables <- as.character(attr(terms(FORM), "variables"))[-1]
  x.label <- variables[2]
  data <- data[,c(x.label,y.label)]
  if (head(class(data[,y.label]),1) != "ordered") {
    data[,y.label] <- factor(data[,y.label])
  }
  data[,y.label] <- factor(data[,y.label],ordered=TRUE,levels=c(level,setdiff(levels(data[,y.label]),level)))
  if(nlevels(data[,y.label])>2) { levels(data[,y.label]) <- c(levels(data[,y.label])[1],rep(paste("Not",level),nlevels(data[,y.label])-1)) }
  x <- data[,x.label]
  y <- data[,y.label]
  
  color <- FALSE
  labelat=c(); xlab=c(); ylab=c(); magnification=1
  complete.x <- which(!is.na(x))
  complete.y <- which(!is.na(y))
  complete.cases <- intersect(complete.x, complete.y)
  x <- x[complete.cases]
  y <- y[complete.cases]
  if (head(class(x),1) != "ordered") {
    x <- factor(x)
  }
  if (head(class(y),1) != "ordered") {
    y <- factor(y)
  }
  data[,1] <- x
  data[,2] <- y
  n <- length(x)
  nx.levels <- length(unique(x))
  ny.levels <- length(unique(y))
  if (nx.levels < 2 | ny.levels < 2) {
    stop(paste("Error:  need at least 2 levels to proceed.  x has", 
               nx.levels, "and y has", ny.levels))
  }
  if (nx.levels > 100) {
    stop(paste("Error:  function not designed for more than 100 levels of x"))
  }
  
  xlevel.names <- levels(x)
  if(length(labelat)>0) { xlevel.names[!(xlevel.names %in% labelat)] <- "" }
  
  ylevel.names <- levels(y)
  CONT.TAB <- table(x, y)
  CONT.TAB <- addmargins(CONT.TAB)
  rownames(CONT.TAB)[nx.levels + 1] <- "Total"
  colnames(CONT.TAB)[ny.levels + 1] <- "Total"
  O <- matrix(table(x, y), nrow = nx.levels, ncol = ny.levels)
  E <- (apply(O, 1, sum) %o% apply(O, 2, sum))/n
  plot(0, 0, col = "white", xlim = c(0, 1.3), ylim = c(-0.05, 
                                                       1), xlab=ifelse(equal==TRUE,"",x.label),cex.main=0.7,main = "", ylab = paste("Probability",y.label,"is",level), axes = FALSE)
  if(equal==FALSE) { axis(1, at = seq(0, 1, 0.05)) }
  axis(2, at = seq(0, 1, 0.05))
  COLORS <- grey(seq(0.1, 0.7, length = ny.levels))
  marginal.y <- apply(O, 2, sum)/n
  break.y <- c(0, cumsum(marginal.y))
  for (i in 1:ny.levels) {
    rect(1.01, break.y[i], 1.11, break.y[i + 1], col = COLORS[i])
    text(1.1, (break.y[i + 1] + break.y[i])/2, ylevel.names[i], 
         srt = 0, pos = 4, cex = 1)
  }
  marginal.x <- apply(O, 1, sum)/n
  
  
  if(equal==FALSE) { break.x <- c(0, cumsum(marginal.x)) } else { break.x <- seq(0,1,length=1+length(xlevel.names)) }
  
  for (i in 1:nx.levels) {
    marginal.y <- O[i, ]/sum(O[i, ])
    break.y <- c(0, cumsum(marginal.y))
    for (j in 1:ny.levels) {
      rect(break.x[i], break.y[j], break.x[i + 1], break.y[j + 
                                                             1], col = COLORS[j])
    }
    if(inside==TRUE) { 
      text((break.x[i + 1] + break.x[i])/2, 0.5, xlevel.names[i],cex=magnification,srt=90,col="white")
    } else {
      text((break.x[i + 1] + break.x[i])/2, -0.05, xlevel.names[i],cex=magnification) }
    
  }
  lines(c(0,1),rep(mean(data[,y.label]==level),2),lwd=2,col="red")
  if(equal==TRUE) { text(0.5,-0.05,x.label)   }
  FORM3 <- formula(paste(y.label,"=='",level,"'~", x.label,sep=""))
  SUMMARY <- aggregate(FORM3,data=data,FUN=mean)
  AOV <- aov(FORM3,data=data)
  TUKEY <- TukeyHSD(AOV)
  LETTERS <- multcompLetters4(AOV,TUKEY) 
  SUMMARY$letters <- LETTERS[[1]][1]$Letters[match(SUMMARY[,1],names( LETTERS[[1]][1]$Letters ) )]
  SUMMARY$n <- as.numeric(table(data[,x.label]))
  names(SUMMARY)[2] <- paste("Prob",level,sep="")
  if(sort==TRUE) { SUMMARY <- SUMMARY[order(SUMMARY[,2],decreasing=TRUE),] }
  rownames(SUMMARY) <- NULL
  print(SUMMARY)
  R2 <- summary(lm( formula, data ) )[8]$r.squared
  cat(paste("\nDriver Score:",round(R2,digits=3),"\n"))
  cat(paste("Scores range between 0 and 1.  Larger scores = stronger driver.\n"))
  cat(paste("Although context dependent, values above 0.02 or so are 'reasonably strong' drivers.\n"))
}

examine_driver_Ynumeric <- function(formula,data,sort=TRUE,lambda=NA) { 
  require(regclass)
  require(multcompView)
  require(stringr)
  
  FORM <- as.formula(formula)
  variables <- as.character(attr(terms(FORM), "variables"))[-1]
  x <- data[,variables[2]]
  y <- data[,variables[1]]
  complete.x <- which(!is.na(x))
  complete.y <- which(!is.na(y))
  complete.cases <- intersect(complete.x, complete.y)
  x <- x[complete.cases]
  y <- y[complete.cases]
  plot(y~x,xlab=variables[2],ylab=variables[1])
  M <- lm(FORM,data)
  if(class(x)[1] %in% c("integer","numeric")) { 
    #if(is.na(lambda)) { SS <- smooth.spline(x,y) } else { SS <- smooth.spline(x,y,lambda) }
    #lines(SS,col="red",lwd=3)
    abline(M,col="blue",lwd=3)
  }
  abline(h=mean(y),col="red")
  if(class(x)[1]  %in% c("integer","numeric") ) {
    print( summary(M) )
  } else { 
    
    SUMMARY <- aggregate(FORM,data=data,FUN=mean)
    AOV <- aov(FORM,data=data)
    TUKEY <- TukeyHSD(AOV)
    LETTERS <- multcompLetters4(AOV,TUKEY) 
    SUMMARY$letters <- LETTERS[[1]][1]$Letters[match(SUMMARY[,1],names( LETTERS[[1]][1]$Letters ) )]
    SUMMARY$n <- as.numeric(table(x))
    names(SUMMARY)[2] <- paste("Avg",variables[1],sep="")
    if(sort==TRUE) { SUMMARY <- SUMMARY[order(SUMMARY[,2],decreasing=TRUE),] }
    rownames(SUMMARY) <- NULL
    print(SUMMARY) }
  
  cat(paste("\nDriver Score:",round(summary(M)[8]$r.squared,digits=3),"\n"))
  cat(paste("Caution: misleading if driver is numerical and trend isn't linear.\n"))
  cat(paste("Scores range between 0 and 1.  Larger scores = stronger driver.\n"))
  cat(paste("Although context dependent, values above 0.02 or so are 'reasonably strong' drivers.\n"))
}

summary(HitPredictionTotal)


HitPredictionTotal[which(HitPredictionTotal$danceability == 0), c("track", "artist", "danceability")]
HitPredictionTotal = HitPredictionTotal %>% filter( danceability != 0)

HitPredictionTotal[which(HitPredictionTotal$chorus_hit == 0), c("track", "chorus_hit")]
HitPredictionTotal = HitPredictionTotal %>% filter( chorus_hit != 0  )

HitPredictionTotal[which(HitPredictionTotal$instrumentalness == 0), c("track", "instrumentalness")]

summary(HitPredictionTotal)


write.csv(HitPredictionTotal,file="Case5-cleanHitPrediction.csv",row.names=FALSE)

HitPredictionTotal = read.csv("Case5-cleanHitPrediction.csv")


HitPredictionClean = HitPredictionTotal

HitPredictionClean$track = NULL
HitPredictionClean$artist = NULL
HitPredictionClean$uri = NULL

mean(HitPredictionClean$target == 1)

HitPredictionFit1 = HitPredictionClean

HitPredictionFit1$Outcome = factor( ifelse(HitPredictionClean$target ==1,"Success","Fail"))

HitPredictionFit1$energy = NULL
HitPredictionFit1$target = NULL
HitPredictionFit1$instrumentalness = NULL

HitPredictionFit1$track = NULL
HitPredictionFit1$artist = NULL
HitPredictionFit1$uri = NULL

TREE = rpart(Outcome ~ . , data=HitPredictionFit1, minbucket=1000, cp=0.0065699157 )
TREE$cptable
visualize_model(TREE)

TREE


mean(HitPredictionClean$target == 1 & HitPredictionClean$duration_ms < 60000*3)

HitPredictionFit2 = HitPredictionClean

HitPredictionFit2$Outcome = factor( ifelse(HitPredictionClean$target == 1 & HitPredictionClean$duration_ms < 60000*3.25,"Success","Fail"))

HitPredictionFit1$energy = NULL
HitPredictionFit1$target = NULL
HitPredictionFit1$instrumentalness = NULL

HitPredictionFit2$duration_ms = NULL
HitPredictionFit2$target = NULL

TREE = rpart(Outcome ~ . , data=HitPredictionFit2, minbucket=800, cp=0.009353997)
TREE$cptable
visualize_model(TREE)

HitPredictionClean = HitPredictionTotal

mean(HitPredictionClean$target == 1)
HitPredictionClean$Outcome = factor( ifelse(HitPredictionClean$target == 1,"Success","Fail"))


summary(HitPredictionClean$target)
summary(HitPredictionClean$energy)
summary(cut(HitPredictionClean$energy, seq(from = 0, to = 1, length.out = 10),include.lowest = TRUE))
HitPredictionClean$energy = cut(HitPredictionClean$energy, seq(from = 0, to = 1, length.out = 10),include.lowest = TRUE)
examine_driver_Ycat(Outcome == "Success" ~ energy, data=HitPredictionClean ) 


summary(HitPredictionClean$danceability)
HitPredictionClean$danceability = cut(HitPredictionClean$danceability, seq(from = 0, to = 1, length.out = 10),include.lowest = TRUE)
examine_driver_Ycat(Outcome=="Success" ~ danceability, data=HitPredictionClean ) 


summary(HitPredictionClean$key)
examine_driver_Ycat(Outcome=="Success" ~ key, data=HitPredictionClean )

summary(HitPredictionClean$loudness)
HitPredictionClean$loudness = cut(HitPredictionClean$loudness, seq(from = -60, to = 0, length.out = 15),include.lowest = TRUE)
examine_driver_Ycat(Outcome=="Success" ~ loudness, data=HitPredictionClean )

summary(HitPredictionClean$mode)
examine_driver_Ycat(Outcome=="Success" ~ mode, data=HitPredictionClean )

summary(HitPredictionClean$speechiness)
HitPredictionClean$speechiness = cut(HitPredictionClean$speechiness, seq(from = 0, to = 1, length.out = 10),include.lowest = TRUE)
examine_driver_Ycat(Outcome=="Success" ~ speechiness, data=HitPredictionClean ) 

summary(HitPredictionClean$acousticness)
HitPredictionClean$acousticness = cut(HitPredictionClean$acousticness, seq(from = 0, to = 1, length.out = 10),include.lowest = TRUE)
examine_driver_Ycat(Outcome=="Success" ~ acousticness, data=HitPredictionClean ) 

summary(HitPredictionClean$instrumentalness)
HitPredictionClean$instrumentalness = cut(HitPredictionClean$instrumentalness, seq(from = 0, to = 1, length.out = 10),include.lowest = TRUE)
examine_driver_Ycat(Outcome=="Success" ~ instrumentalness, data=HitPredictionClean ) 

summary(HitPredictionClean$liveness)
HitPredictionClean$liveness = cut(HitPredictionClean$liveness, seq(from = 0, to = 1, length.out = 10),include.lowest = TRUE)
examine_driver_Ycat(Outcome=="Success" ~ liveness, data=HitPredictionClean ) 

summary(HitPredictionClean$valence)
HitPredictionClean$valence = cut(HitPredictionClean$valence, seq(from = 0, to = 1, length.out = 10),include.lowest = TRUE)
examine_driver_Ycat(Outcome=="Success" ~ valence, data=HitPredictionClean ) 

summary(HitPredictionClean$tempo)
HitPredictionClean$tempo = cut(HitPredictionClean$tempo, seq(from = 30, to = 242, length.out = 15),include.lowest = TRUE)
examine_driver_Ycat(Outcome=="Success" ~ tempo, data=HitPredictionClean ) 

summary(HitPredictionClean$duration_ms)
HitPredictionClean$duration_ms = cut(HitPredictionClean$duration_ms, seq(from = 0, to = 60000*10, length.out = 10),include.lowest = TRUE)
examine_driver_Ycat(Outcome=="Success" ~ duration_ms, data=HitPredictionClean ) 

summary(HitPredictionClean$time_signature)
examine_driver_Ycat(Outcome=="Success" ~ time_signature, data=HitPredictionClean ) 

summary(HitPredictionClean$chorus_hit)
HitPredictionClean$chorus_hit = cut(HitPredictionClean$chorus_hit, seq(from = 0, to = 450, length.out = 10),include.lowest = TRUE)
examine_driver_Ycat(Outcome=="Success" ~ chorus_hit, data=HitPredictionClean ) 

summary(HitPredictionClean$sections)
examine_driver_Ycat(Outcome=="Success" ~ sections, data=HitPredictionClean ) 

write.csv(HitPredictionClean,file="Case5-cleanHitPrediction2.csv",row.names=FALSE)
mean(HitPredictionClean$Outcome == "Success")


Case5_HitPrediction[which(Case5_HitPrediction$energy >= .7 & Case5_HitPrediction$target == 1),]
TREE
