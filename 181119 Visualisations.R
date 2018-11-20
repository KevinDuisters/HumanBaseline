#---------------------------------------------------------------------------------------------------------#
# TITLE:
# Human Baseline Study
# Kevin Duisters, Mathematical Institute Universiteit Leiden
#---------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------#
# DISCLAIMER:
# Data is confidential and property of LACDR - Universiteit Leiden
# All rights reserved
#---------------------------------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------------------------------#
# load packages
library(readxl)
library(ggplot2)
library(gplots) # heatmap
library(RColorBrewer)


#---------------------------------------------------------------------------------------------------------#
# Load data (published set from Human Baseline Tool)
ST1<-read_excel("data/ST1.xlsx")

# Tables and mention in abstract
ST1[order(ST1$sigma_b+ST1$sigma_eps)[1:10],c(1,2,18,19,20)] # Table 2
ST1[order(ST1$sigma_time)[1:10],c(1,2,18,19,20)] # Table 3
ST1[order(ST1$sigma_time,decreasing=T)[1:10],c(1,2,18,19,20)] # for heatmap Fig 3


#-----------------------------------------------------------------------------------------------------------

# Figures Variance and subVariance (supplemental fig) 


#sub <- F;pdf("Figures/Variance.pdf",width=9)
#sub<-T;k = 1;pdf("Figures/subVariance1.pdf",width=9)
sub<-T;k = 2;pdf("Figures/subVariance2.pdf",width=9)
#sub<-T;k = 3;pdf("Figures/subVariance3.pdf",width=9)


# sort to make bubbles look nice on top of each other
panelnames <- c("SomaLogic","Metabolon","BMFL")
mysort <- c(which(ST1$panel==(panelnames[1])),which(ST1$panel==(panelnames[2])),which(ST1$panel==(panelnames[3])))
mylabs <- c("Plasma Proteins","Plasma Metabolon","Urine Amines")
if(sub==F){dynamictitle<-""}else{dynamictitle<-paste0(mylabs[k]," [",sum(ST1$panel==(panelnames[k])),"]")} # plot title
mycols <- c("lightblue","white","navy") # rgb Leiden Univ color #rgb(0,17,88,maxColorValue = 255)
Platform <- numeric(length(ST1$panel))
Platform[which(ST1$panel=="SomaLogic")] <- 1
Platform[which(ST1$panel=="Metabolon")] <- 2
Platform[which(ST1$panel=="BMFL")] <- 3
Platform <- Platform[mysort]
Platform <- factor(Platform,levels=c(1,2,3),labels=mylabs)


mysb <- ST1$sigma_b[mysort] + (sub==T)*(Platform!=mylabs[k])*100 # move VCs not in subplot out of plotted region
V <- data.frame(sigmab=mysb,sigmaeps=ST1$sigma_eps[mysort],sigmatime=ST1$sigma_time[mysort],name=ST1$compound[mysort],platformV=Platform)


gout <- ggplot(V, aes(x=sigmaeps, y=sigmab, size=(sigmatime))) +
          labs(title=dynamictitle,subtitle="")+
          xlim(0,2.2)+ylim(0,2.2)+
            geom_point(shape=21,colour="black",alpha=0.8,aes(fill=platformV))+
            scale_fill_manual(name="Platform",values=mycols) +
          scale_size_continuous(range=c(0, 13)) +    # range can be used to scale the bubbles
          ylab(expression(paste("Between ",(sigma[b])))) +
          xlab(expression(paste("Noise ",(sigma[epsilon])))) +
          labs( size = expression(paste("Time ",(sigma[tau])))) +
          theme_bw() +
          theme(
            plot.title = element_text(size=15,hjust = 0.5,face="bold"),
            text = element_text(size=11),
            legend.position = c(.85, .1),
            legend.justification = c("center", "bottom"),
            panel.border = element_blank()
          ) +
          guides(fill = guide_legend(override.aes = list(shape=21,size=12)))
print(gout)
dev.off()
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Fig 2: Profiles
#-----------------------------------------------------------------------------------------------------------
# prep
numtimes  <- c(-2,-0.5,0.5,1.5,3,4,4.5,6,8,9.5,11,12)
shortnames <- ST1$compound
for(j in 1:nrow(ST1)){
  if(nchar(ST1$compound[j])>20){shortnames[j] <- paste0(substr(ST1$compound[j],1,20),"...")}
}
plat<-ST1$panel
for(j in 1:length(plat)){
  if(plat[j]=="Metabolon"){plat[j] <- "Plasma Metabolites"}
  if(plat[j]=="BMFL"){plat[j] <- "Urine Amines"}
  if(plat[j]=="SomaLogic"){plat[j] <- "Plasma Proteins"}
}

# Some profile illustrations
N <- 10
agei <- 22
bmii <- 20
alpha <- 0.05
# one may do this for other ages/bmi as well (Human Baseline Tool)



sel <- c(206,1941,1817,798) # this is sensitive, look up manually
dynamic <- paste0("Figures/profiles",1:length(sel),".pdf")

for(k in 1:length(sel)){
pdf(dynamic[k],width=9)

  j <- sel[k]
coefs <- as.numeric(ST1[j,4:17])
means <- coefs[-c(1,2)] + (agei)*coefs[1] + (bmii)*coefs[2]
sds <- as.numeric(ST1[j,21:32])

localpha <- max(1e-6,alpha)
qa <- sqrt(1+1/N)*qt(1-localpha/2,df=N-1)
high <- means + qa*sds
low <- means-qa*sds

isna <- (is.na(means)==F)
aux <- data.frame(low=low[isna],mid=means[isna],high=high[isna],time=numtimes[isna])


pnascol <- c(207,214,235)/256


print(ggplot(aux,aes(x=time,y=mid))+
  #ylab("normalized relative concentration")+
    ylab("")+
  xlab("Time (hrs) since breakfast")+
  xlim(-2,12)+
  ylim(-5,5)+
  labs(title=paste0(shortnames[j]," (",plat[j],")"),
       subtitle=bquote(paste(alpha,"=",.(localpha),", ","Age=",.(agei),", BMI=",.(bmii))))+
  geom_line(aes(x=time,y=high),linetype="dashed")+
  geom_line(aes(x=time,y=low),linetype="dashed")+
  geom_ribbon(aes(x=time,ymax=high,ymin=low),fill=rgb(pnascol[1],pnascol[2],pnascol[3]),alpha=0)+
  geom_segment(aes(x=time,xend=time,y=low,yend=high),size=2)+
  geom_line(colour="black",linetype="dashed")+
  geom_point(colour="black",size=5)+
  theme_bw()+
  theme(
    axis.text=element_text(size=25),
    axis.title=element_text(size=25),
    plot.title = element_text(size=30,hjust = 0.5,face="bold"),
    plot.subtitle = element_text(size=25,hjust = 0.5,face="bold"),
    panel.border=element_blank()
  ))


dev.off()
}


#-----------------------------------------------------------------------------------------------------------
# Discussion example
# Lactate
M <- read_excel("data/data.xlsx",sheet=1)
lactate <- scale(M[,which(colnames(M)=="lactate")])
rm(M)
q90 <- quantile(lactate,0.9)



agei <- 22
bmii <- 20
j<-which(ST1$compound=="Lactate")
colnames(ST1)[4:17]
coefs <- as.numeric(ST1[j,4:17])
mj <- coefs[-c(1,2)] + (agei)*coefs[1] + (bmii)*coefs[2]
sj <- as.numeric(ST1[j,21:32])

1-pnorm((0.7-mj[4])/sj[4])
1-pnorm((2-mj[8])/sj[8])
1/(dnorm((2-mj[8])/sj[8])/dnorm((0.7-mj[4])/sj[4]))

dnorm((2-1.42)/1.01)/dnorm((0.7+0.584)/0.415)


qnorm(0.9,-0.584,0.415)
qnorm(0.9,mj[5],sj[5])
qnorm(0.9,mj[11],sj[11])




#-----------------------------------------------------------------------------------------------------------



#-----------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------
# Fig 3: heatmap
totaltimes <- c("-2 hr","-0.5 hr","0.5 hr","1.5 hr","3 hr","4 hr","4.5 hr","6 hr","8 hr","9.5 hr","11 hr","12 hr") # excl 24 hr
top <- order(ST1$sigma_time,decreasing=T)[1:10]
set <- matrix(NA,length(top),ncol=12,dimnames = list(row=top,col=totaltimes))
for(i in 1:length(top)){
set[i,] <-  as.numeric(ST1[top[i],6:17])+23.8*as.numeric(ST1[top[i],4])+23.3*as.numeric(ST1[top[i],5])
}

labs <- as.data.frame(ST1)[order(ST1$sigma_time,decreasing=T)[1:10],1]


pdf("Figures/heat.pdf")
hm <- heatmap.2(set,trace="none",scale="row",density.info="none",adjCol=c(0.5,1),
                                sepcolor="white",rowsep=0:10,colsep=0:12,
                                col=brewer.pal(10,"RdBu"),srtCol=0,dendrogram="none",Rowv=F,Colv=F,key=F,cexRow = 1,
                                lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(2,5,2),lwid=c(2,7,2),
                cexCol=1,labRow=labs
                )


#locator(1)
ep<-0.18
for(j in 1:length(hm$col)){
  rect(-0.2+ep,0.4+(j-1)*0.025,-0.175+ep,0.425+(j-1)*0.025,col=hm$col[j],xpd=T)
}

text(x=-0.22+ep,y=0.52,"row zscore",xpd=T,srt=90,cex=0.9)
text(x=-0.15+ep,y=0.4,"-2",xpd=T,srt=90,cex=0.9)
text(x=-0.15+ep,y=0.525,"0",xpd=T,srt=90,cex=0.9)
text(x=-0.15+ep,y=0.65,"2",xpd=T,srt=90,cex=0.9)

dev.off()
#-----------------------------------------------------------------------------------------------------------
# Figure 4 (p values) in separate R file "... pvalue visuals.R"
#-----------------------------------------------------------------------------------------------------------
#
#-----------------------------------------------------------------------------------------------------------
# Figure 5 boxplot VC
# Repro Kim et al (2014)
pdf("Figures/relVCbox.pdf",width=20)
panelnames <- c("SomaLogic","Metabolon","BMFL") # C, A, B


pi.tab <- matrix(NA,3,5)
#par(mfrow=c(1,3))
par(mfrow=c(1,3),mar=c(6,8,6,2)+0.1)
for(k in c(2,3,1)) {
  loc <- (ST1$panel==(panelnames[k]))
  VCk <- data.frame(between=(ST1$sigma_b[loc])^2,time=(ST1$sigma_time[loc])^2,noise=(ST1$sigma_eps[loc])^2)
  boxplot(VCk/rowSums(VCk),ylim=c(0,1),main="",xaxt="n",cex.axis=3) # relative proportions
  axis(side=1,at=c(1,2,3),colnames(VCk),cex.axis=3,line=1,tick=F)
  title(paste0((LETTERS[c(3,1,2)])[k],": ",mylabs[k]," [",sum(ST1$panel==(panelnames[k])),"]"),cex.main=4,line=3)
  pi.tab[k,] <- c(colMeans(VCk/rowSums(VCk)),median((VCk$between)/((VCk$between)+VCk$time+VCk$noise)),median((VCk$between)/(VCk$between+VCk$noise)))
}
dev.off()


# Urine / plasma comparison table
round(100*pi.tab[,c(2,1,3,4,5)],0)


#-----------------------------------------------------------------------------------------------------------
# Figure 6
# Variance visualisation per pathway
paths <- as.data.frame(read_excel("data/pathways.xlsx")) # this is Metabolon only!
sel <- (ST1$panel=="Metabolon" ) 
superpath <- paths[match(ST1$compound[sel],paths[,1]),2]
subpath <- paths[match(ST1$compound[sel],paths[,1]),2]



pdf("Figures/famVariance.pdf",width=9)

sort <- c(1,3,4,5,7,8,2,6,5)
us <- unique(superpath)[sort]
mycols <- brewer.pal(length(us),"RdBu")


V <- data.frame(familyV=superpath,sigmab=ST1$sigma_b[sel],sigmaeps=ST1$sigma_eps[sel],sigmatime=ST1$sigma_time[sel])
Vm <- matrix(NA,length(us),4)
for(i in 1:length(us)){
  Vm[i,2:4] <- apply(V[V$familyV==us[i],2:4],2,median)
}  
V <- data.frame(familyV=us,sigmab=Vm[,2],sigmaeps=Vm[,3],sigmatime=Vm[,4])

gout <- ggplot(V, aes(x=sigmaeps, y=sigmab, size=(sigmatime))) +
  labs(title="",subtitle="")+
  xlim(0,1.0)+ylim(0,1.0)+
  geom_point(shape=21,colour="black",alpha=0.8,aes(fill=familyV),show.legend = FALSE)+
  geom_point(shape=21,colour="black",alpha=0.8)+
  scale_fill_manual(name="Family",values=mycols) +
  scale_size_continuous(range=c(0, 13)) +    # range can be used to scale the bubbles
  geom_text(size=4,aes(label=familyV) ,nudge_x = 0.0,nudge_y = 0.0)+
  ylab(expression(paste("Between ",(sigma[b])))) +
  xlab(expression(paste("Noise ",(sigma[epsilon])))) +
  labs( size = expression(paste("Time ",(sigma[tau])))) +
  theme_bw() +
  theme(
    plot.title = element_text(size=15,hjust = 0.5,face="bold"),
    text = element_text(size=11),
    legend.position = c(.15, .1),
    legend.justification = c("center", "bottom"),
    panel.border = element_blank()
  ) +
  guides(fill = guide_legend(override.aes = list(shape=21,size=12)))
print(gout)


dev.off()
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Supplemental Figures
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Repro Maitre et al (2017)
# ICC on urine
# TO BE UPDATED
#-----------------------------------------------------------------------------------------------------------
#pdf("Figures/ICC.pdf")
icccols <- c("black",cols[3],cols[5],cols[6])
swmat<-sqrt(slist[[3]]^2-s2blist[[3]])
ICC <- sqrt(s2blist[[3]])/(sqrt(s2blist[[3]]) + swmat)
srtd <- order(rowMeans(ICC))

par(xpd=F)
par(mar=c(5,13.5,4,2)+0.1)
plot(0,xlim=c(0,1),ylim=c(1,63),type="n",ylab="",xlab="",yaxt="n",xaxt="n",bty="n")
axis(side=1,seq(0,1,0.1),tick=F,line=-1)
mtext(side=1,at=0.5,"Intra-class correlation coefficient (ICC)",line=1)
segments(y0=rep(1,10),y1=rep(63,10),x0=seq(0,1,0.1),x1=seq(0,1,0.1),col="grey")
segments(x0=rep(0,63),x1=rep(1,63),y0=seq(1,62,1)-1/2,y1=seq(1,62,1)-1/2,col="grey")
polygon(x=c(0,1,1,0,0),y=c(0,0,63,63,0),lwd=2)
points(seq(1,62,1),x=ICC[srtd,1],col=icccols[1])
points(seq(1,62,1),x=ICC[srtd,2],pch=22,bg=icccols[2])
points(seq(1,62,1),x=ICC[srtd,3],pch=21,bg=icccols[3])
points(seq(1,62,1),x=ICC[srtd,4],pch=24,bg=icccols[4])
par(xpd=T)
text(x=rep(-0.05,62),y=seq(1,62,1),namelist[[3]][srtd],adj=1,cex=0.6)

legend(x=0.58,y=8.5,c("-2 hr","4 hr","8 hr","12 hr"),pt.bg=icccols,pch=c(1,22,21,24),
       ncol=2,bg="white")


#dev.off()

#-----------------------------------------------------------------------------------------------------------
# TO BE UPDATED
#Relative urine proportions
#pdf("Figures/relbars.pdf")
barcols <- c(cols[6],cols[2])
VC <- data.frame(between=s2blist[[3]],within=s2wlist[[3]])
rownames(VC) <- namelist[[3]]
par(mar=c(5,13.5,4,2)+0.1)
barplot(t(as.matrix(VC/rowSums(VC))[order((VC/rowSums(VC))[,1]),]),beside=F,horiz=T,xaxt="n",las=1,col=barcols) # relative proportions
legend(x=0.45,y=8,c("between","within"),horiz=T,pch=c(15,15),col=barcols,bg="white")
axis(side=1,at=seq(0,1,0.2),labels=paste0(seq(0,100,20),"%"))
dev.off()


#-----------------------------------------------------------------------------------------------------------