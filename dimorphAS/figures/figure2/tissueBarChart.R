
## Data for number of tissues affected by an alternative splicing event


npgBlue<- rgb(60/256,84/256,136/256,1)
npgRed <- rgb(220/256,0,0,0.8)
npgGreen <- rgb(0,160/256,135/256,1)
npgBrown <- rgb(126/256,97/256,72/256,1)

colrs <- c(npgRed,npgBlue,npgGreen)

slices <-c(3273,459,9)
lbls <-c("1 tissue","2-5 tissues",">5 tissues")
pct <- round(slices/sum(slices)*100, digits=1)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=colrs,main="",cex=1.5,xaxs="r",yaxs="r")


