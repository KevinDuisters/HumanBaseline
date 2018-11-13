# ggplot


library(readxl)
library(ggplot2)

# read data
D <- as.data.frame(read_excel("dat.xlsx"))
mysort <- c(which(D$panel=="SomaLogic"),which(D$panel=="Metabolon"),which(D$panel=="BMFL"))
mylabs <- c("Proteins","Plasma Metabolon","Urine BMFL")
mycols <- c("lightblue","white","navy")
Platform <- numeric(length(D$panel))
Platform[which(D$panel=="SomaLogic")] <- 1
Platform[which(D$panel=="Metabolon")] <- 2
Platform[which(D$panel=="BMFL")] <- 3
Platform <- Platform[mysort]
Platform <- factor(Platform,levels=c(1,2,3),labels=mylabs)
str(Platform)

V <- data.frame(sigmab=D$sigmab[mysort],sigmaw=D$sigmaw[mysort],sigmatau=D$sigmatau[mysort])
rm(D)

pdf("testplot.pdf",width=9)
ggplot(V, aes(x=sigmaw, y=sigmab, size=(sigmatau))) +
  xlim(0,2)+ylim(0,2)+
  geom_point(shape=21,colour="black",alpha=0.8,aes(fill=Platform))+
  scale_fill_manual(name="Platform",values=mycols) +
  scale_size_continuous(range=c(0, 13)) +    # range can be used to scale the bubbles
  ylab(expression(paste("Between ",(sigma[b])))) +
  xlab(expression(paste("Noise ",(sigma[epsilon])))) +
  labs( size = expression(paste("Time ",(sigma[tau])))) +
  theme_bw() +
  theme(
    text = element_text(size=11),
    legend.position = c(.85, .3),
    legend.justification = c("center", "bottom"),
    panel.border = element_blank()
  )+ 
  guides(fill = guide_legend(override.aes = list(shape=21,size=12)))
dev.off()

#----------------------------------------------------------------------#
# Lines
low <- rnorm(10,1,1)
mid <- low + 0.5
high <- mid + 0.5
dynamicname <- "Prolylglycine"
dynamicage <- 22
dynamicBMI <- 20

D <- data.frame(low=low,mid=mid,high=high,time=1:10)

#pdf("testprofile.pdf",width=9)
ggplot(D,aes(x=time,y=mid))+
  ylab("Value")+
  xlab("Time (hrs) since breakfast")+
  ylim(-2,4)+
  ggtitle(paste0(dynamicname,"\n","Age = ",dynamicage,", BMI = ",dynamicBMI))+
  geom_ribbon(aes(x=time,ymax=high,ymin=low),fill="grey",alpha=0.5)+
  geom_segment(aes(x=time,xend=time,y=low,yend=high))+
  geom_point()+
  geom_line(colour="black")+
  theme_bw()+
  theme(
    plot.title = element_text(size=15,hjust = 0.5),
    panel.border=element_blank()
  )
#dev.off()