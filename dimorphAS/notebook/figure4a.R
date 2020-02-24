# This script creates figure 4a.
# please run the following command first
#$ perl parseMT.pl
# this creates the files needed for figure 4a and 4b


# Removal all variables from workspace
rm(list=ls())


dat <- read.table("lv.txt", header=FALSE, sep = "\t")
colnames(dat) <- c("RBP","Expression")

library(ggplot2)
library(ggsci)
library(grid)
mypal = pal_npg("nrc", alpha = 0.7)(9)
d2<-dat[dat$Expression!=0,]
lm_fit <- lm(d2$Expression ~ d2$RBP, data=d2)
LM<-summary(lm_fit)
rsquared<-round(LM$r.squared,digits=2)  

# save predictions of the model in the new data frame 
# together with variable you want to plot against
predicted_df <- data.frame(expr_pred = predict(lm_fit, d2), RBP=d2$RBP)


p<-ggplot(dat, aes(x=RBP, y=Expression)) + geom_point(shape=21,fill = mypal[3],size=3) +  theme_bw()
#+ scale_fill_npg() 
p <- p + theme(axis.text = element_text(size=32, hjust=0.5),
               axis.title.x=element_text(size=24),
               axis.title.y = element_text(size=24),
               axis.text.y = element_text(size=32),
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank()) 
p <- p +  geom_hline(yintercept=0, linetype="dashed", color = mypal[4])
p <- p +xlab('\U27F6 \n Sum of RBP effect magnitude')+ylab('Expression\ninclusion \U27F5 effect \U27F6 skipping')
p <- p+ geom_line(color='red',data = predicted_df, aes(y=expr_pred, x=RBP))
mylabel<-paste(italic(r)^2~"="~rsquared) 
p <- p+ geom_text(x = 3, y = 0.45, label = as.character(paste( "r^2==",rsquared)), size=7, parse = TRUE)
p


ggsave(file = "fig4a.pdf", plot = p)
