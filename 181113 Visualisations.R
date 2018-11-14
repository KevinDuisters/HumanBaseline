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


#---------------------------------------------------------------------------------------------------------#
# Load data (published set from Human Baseline Tool)
ST1<-read_excel("data/ST1.xlsx")

# General info
#totaltimes <- c("-2 hr","-0.5 hr","0.5 hr","1.5 hr","3 hr","4 hr","4.5 hr","6 hr","8 hr","9.5 hr","11 hr","12 hr") # excl 24 hr




#-----------------------------------------------------------------------------------------------------------

# Figures Variance and subVariance (supplemental fig) 


sub <- F;pdf("Figures/Variance.pdf",width=9)
sub<-T;k = 1;pdf("Figures/subVariance1.pdf",width=9)
sub<-T;k = 2;pdf("Figures/subVariance2.pdf",width=9)
sub<-T;k = 3;pdf("Figures/subVariance3.pdf",width=9)


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
  if(nchar(ST1$compound[j])>30){shortnames[j] <- paste0(substr(ST1$compound[j],1,30),"...")}
}

# Some profile illustrations
N <- 10
agei <- 22
bmii <- 20
alpha <- 0.05

# 4x4 plot; new version ggplot, save them separately and make 4x4 in LaTeX using minipage




# Prolylglycine 194 (in Metab panel)
# Lactate 206 (in Metab panel)


# Discussion example
#1-pnorm((0.7-mi[4])/s[4])
#1-pnorm((2-mi[8])/s[8])
#dnorm((2-mi[8])/s[8])/dnorm((0.7-mi[4])/s[4])
#dnorm((2-1.41)/1.01)/dnorm((0.7+0.584)/0.415)
#qnorm(0.9,mi[4],s[4])
#qnorm(0.9,mi[8],s[8])

# Valine (Urine) (50 in urine panel)

# CD27 (1130 in protein panel)

# one may do this for other ages/bmi as well (Human Baseline Tool)

sel <- c(14,1,32,28) # this is sensitive, look up manually
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
  #ylab("standardized intensity")+
    ylab("")+
  xlab("Time (hrs) since breakfast")+
  xlim(-2,12)+
  ylim(-5,5)+
  labs(title=paste0(shortnames[j]," (",ST1$panel,", ",ST1$sample[j],")"),
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



#-----------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------
# Fig 3: heatmap
# To BE UPDATED

top <- (stau > quantile(stau,0.99))
subtop <- ((1:sum(pvec))[top])[order(stau[top],decreasing=T)][1:10]
cumsum(pvec)


# This is also in Supp Table 1
coeftotmat <- matrix(NA,nrow=sum(pvec),ncol=2+length(totaltimes),dimnames=list(fullnames,c("age","bmi",totaltimes)))
coeftotmat[1:pvec[1],c(1,2,4,5,6,7,9,10,12,13,14)] <- unlist(coeflist[[1]])[,c(10,11,1:9)]
coeftotmat[(cumsum(pvec)[1]+1):(cumsum(pvec)[2]),c(1,2,4,5,6,7,9,10,12,13,14)] <- unlist(coeflist[[2]])[,c(10,11,1:9)]
coeftotmat[(cumsum(pvec)[2]+1):(cumsum(pvec)[3]),c(1,2,3,8,11,14)] <- unlist(coeflist[[3]])[,c(5,6,1:4)]
coeftotmat[(cumsum(pvec)[3]+1):(cumsum(pvec)[4]),c(1,2,4,5,6,7,10,12,14)] <- unlist(coeflist[[4]])[,c(8,9,1:7)]

mmat <- (coeftotmat+23.8*coeftotmat[,1]+23.3*coeftotmat[,2])[,-c(1,2)] 




set <- mmat[subtop,]
labs <- c("2-piperidinone","tartronate (hydroxymalonate)","prolylglycine","N6,N6,N6-Trimethyl-lysine","N-(2-furoyl)glycine",
          "acetylcarnitine","3-hydroxybutyrylcarnitine", "5,6-dihydrothymine","5-hydroxylysine","N-palmitoylglycine",
          "palmitoleoylcarnitine (C16:1)","sucrose","adrenate (22:4n6)","fructose","N2,N2-dimethylguanosine",
          "N4-acetylcytidine","1-arachidonoyl-GPC (20:4)",expression(paste(5,alpha,"-pregnan-",3,beta,20,alpha,"-diol disulfate")), 
          "palmitoleate (16:1n7)","1-lignoceroly-GPC (24:0)"
          )[1:10]
          # manual input! clean names


#pdf("Figures/heat.pdf")
#hm <- heatmap.2(set,trace="none",scale="row",density.info="none",adjCol=c(0.5,1),
#                sepcolor="white",rowsep=1:20,colsep=1:9,
#                col=brewer.pal(10,"RdBu"),srtCol=0,dendrogram="none",Rowv="",Colv="",key=F,cexRow = 1,
#                lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(1.5,25,0.5),lwid=c(2,7,2),
#                cexCol=1,labRow=labs)

hm <- heatmap.2(set,trace="none",scale="row",density.info="none",adjCol=c(0.5,1),
                                sepcolor="white",rowsep=0:10,colsep=0:12,
                                col=brewer.pal(10,"RdBu"),srtCol=0,dendrogram="none",Rowv="",Colv="",key=F,cexRow = 1,
                                lmat=rbind(c(0,3,0),c(2,1,0),c(0,4,0)),lhei=c(2,5,2),lwid=c(2,7,2),
                                cexCol=1,labRow=labs)
?heatmap.2

#locator(1)
ep<-0.18
for(j in 1:length(hm$col)){
  rect(-0.2+ep,0.4+(j-1)*0.025,-0.175+ep,0.425+(j-1)*0.025,col=hm$col[j],xpd=T)
}

text(x=-0.22+ep,y=0.52,"row zscore",xpd=T,srt=90,cex=0.9)
text(x=-0.15+ep,y=0.4,"-2",xpd=T,srt=90,cex=0.9)
text(x=-0.15+ep,y=0.525,"0",xpd=T,srt=90,cex=0.9)
text(x=-0.15+ep,y=0.65,"2",xpd=T,srt=90,cex=0.9)

#dev.off()

#-----------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------
# Advanced covariance visualisation in profile
# Instead of polygon connecting pointwise prediction intervals
library(mvtnorm)
cj <- coeflist[[1]][10,]
s <- slist[[1]][10,]
nm <- namelist[[1]][10]
tvec <- numtimes[which(totaltimes%in%names(cj))]
mi <- cj[-which(names(cj)%in%c("AGE","BMI"))] + agei*cj["AGE"] + bmii*cj["BMI"]

plot(tvec,tvec,ylim=c(-3,4.5),ylab="",xlab="Time (hrs) since breakfeast",type="n",main=nm)
segments(x0=tvec,x1=tvec,y0=mi-ca*s,y1=mi+ca*s,lwd=2)
lines(tvec,mi,lwd=2,type="b",pch=16)


mat <- rmvnorm(1e3,mi,Vj) # generate this Vj manually!
lines(tvec,apply(mat,2,function(x)quantile(x,1-alpha/2)))

tgrid <- seq(6,9.5,by=0.1)

imat <- matrix(NA,nrow(mat),length(tgrid))
for(r in 1:nrow(mat)){
  imat[r,] <- mat[r,6]+(mat[r,7]-mat[r,6])/3.5*(tgrid-6)
}

lines(tgrid,apply(imat,2,function(x)quantile(x,1-alpha/2)),col="red")
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Figure 4 (p values) in separate R file "... pvalue visuals.R"
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
# Figure 5 boxplot VC
# Repro Kim et al (2014)
pdf("Figures/relVCbox.pdf",width=20)
panelnames <- c("SomaLogic","Metabolon","BMFL") # C, A, B


pi.tab <- matrix(NA,K,5)
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


#-----------------------------------------------------------------------------------------------------------
# this is outdated, remove this


BMFL.id.all <- which(panels=="BMFL"&samples=="plasma") #(pvec[1]+1):(pvec[1]+pvec[2])
#cbind(1:pvec[2],fullnames[BMFL.id.all])

# manual! Amy has made match based on urine names, KD converted to plasma (minor loss). 38 out of 55 plasma amines BMFL measured in both
Amy.match.plasma <- c("glycine", "serine",	"threonine",		"alanine",		"aspartate",	"asparagine",	"glutamate",		"glutamine",		"histidine",	"3-methylhistidine",
                      "lysine"	,"N6,N6,N6-trimethyllysine"	,"5-hydroxylysine",	"2-aminoadipate",	"phenylalanine",	"tyrosine",	"3-methoxytyrosine",
                      "tryptophan",	"kynurenine",		"serotonin",	"leucine",	"isoleucine",		"valine",	"methionine","methionine sulfone",	"methionine sulfoxide",	"cysteine",	"taurine",
                      "arginine",	"ornithine",	"citrulline",		"homocitrulline",	"proline",	"dimethylarginine (ADMA + SDMA)",		"hydroxyproline",	"gamma-glutamylalanine",	
                      "phosphoethanolamine (PE)",		"beta-alanine")  
dup <- c(15,39,40,22,26,25,27,28,29,3,34,45,19,20,37,42,2,41,32,54,33,31,43,35,44,36,8,55,24,49,6,17,38,52,21,13,47,5)
BMFL.id <- BMFL.id.all[dup]
Metabolon.id <- sapply(Amy.match.plasma,function(i)which(fullnames[panels=="Metabolon"]==i))
checknames <- cbind(fullnames[BMFL.id],fullnames[Metabolon.id]) # check


# Plots
#pdf("SI Plasma Amines.pdf")
par(xpd=F,mfrow=c(3,3))
cx <- 1
plot(x=sb[BMFL.id],y=sb[Metabolon.id],ylab="Metabolon",xlab="BMFL",main=expression(paste(sigma[b]," (Between)")),xlim=c(0,1),ylim=c(0,1),pch=16,col=ifelse(correl>0.6,"navy","red"),asp=1)
lines(x=c(-100,100),c(-100,100),lty=3)
plot(x=sw[BMFL.id],y=sw[Metabolon.id],ylab="Metabolon",xlab="BMFL",main=expression(paste(sigma[epsilon]," (Noise)")),xlim=c(0,1),ylim=c(0,1),pch=16,col=ifelse(correl>0.6,"navy","red"),asp=1)
lines(x=c(-100,100),c(-100,100),lty=3)
plot(x=stau[BMFL.id],y=stau[Metabolon.id],ylab="Metabolon",xlab="BMFL",main=expression(paste(sigma[tau]," (Time)")),xlim=c(0,1),ylim=c(0,1),pch=16,col=ifelse(correl>0.6,"navy","red"),asp=1)
lines(x=c(-100,100),c(-100,100),lty=3)

# simple correlation analysis on estimated means
#meandif.panels <-  coeftotmat[Metabolon.id,c(4,5,6,7,9,10,12,13,14)] - coeftotmat[BMFL.id,c(4,5,6,7,9,10,12,13,14)]
metabmeans <- coeftotmat[Metabolon.id,c(4,5,6,7,9,10,12,13,14)]+matrix(rep(sapply(Metabolon.id,function(k) c(mean(AGE),mean(BMI))%*%(coeftotmat[k,c(1,2)])),each=9),ncol=9,byrow=T) #metabmeans[1,];mlist[[1]][Metabolon.id[1],] # check
bmflmeans <- coeftotmat[BMFL.id,c(4,5,6,7,9,10,12,13,14)]+matrix(rep(sapply(BMFL.id,function(k) c(mean(AGE),mean(BMI))%*%(coeftotmat[k,c(1,2)])),each=9),ncol=9,byrow=T)
meandif.panels <- metabmeans - bmflmeans

# correl
correl <- matrix(sapply(1:length(BMFL.id),function(k) cor(coeftotmat[Metabolon.id[k],c(4,5,6,7,9,10,12,13,14)],coeftotmat[BMFL.id[k],c(4,5,6,7,9,10,12,13,14)])),ncol=1,dimnames=list(fullnames[BMFL.id],"cor")) 



plot(correl,col=ifelse(correl>0.6,"navy","red"),pch=16,xlab="Amine",ylab="Pearson",main="Correlation between mean profiles",cex.main=cx,xaxt="n")
abline(h=0.6,lty=3)
par(xpd=T)
text(x=which(correl<0.6)-2,y=correl[correl<0.6]+0.05,fullnames[Metabolon.id][correl<0.6],cex=0.5)
par(xpd=F)

# good
#plot(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),type="n",ylab="Mean profile",ylim=c(-2,2),xlab="Time (hrs) since breakfast",main="Asparagine (correlation 0.94)",cex.main=cx)
#lines(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),metabmeans[6,],col="navy")
#lines(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),bmflmeans[6,],lty=2,lwd=2,col="navy")
#text(x=2.25,y=-1,"Metabolon")
#text(x=8.5,y=0.3,"BMFL")

# bad
#plot(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),type="n",ylab="Mean profile",ylim=c(-2,2),xlab="Time (hrs) since breakfast",main="Beta-Alanine (correlation -0.68)",cex.main=cx)
#lines(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),metabmeans[38,],col="red")
#lines(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),bmflmeans[38,],lty=2,lwd=2,col="red")
#text(x=2.25,y=1.2,"Metabolon")
#text(x=1.5,y=-1.25,"BMFL")

# plot raw data
D <- read_xlsx("data/Final data normalized.xlsx",sheet=1) 
D <- D[D$TIME!="24 hr",]
ID <- D$ID
Yj <- scale(D[,4+628])
plot(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),type="n",ylab="Raw rescaled data",ylim=c(-4,4),xlab="Time (hrs) since breakfast",main="Beta-Alanine Metabolon",cex.main=cx)
for(i in 1:10){
  lines(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),Yj[ID==i],col=1+i)
}
lines(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),metabmeans[38,],lwd=2,col="black")

D <- read_xlsx("data/Final data normalized.xlsx",sheet=1) 
D <- D[D$TIME!="24 hr",]
ID <- D$ID
colnames(D)


#---------------------------------------------------------------------------------------------------------#
# Load data (published set from Human Baseline Tool)
BM <- read_xlsx("data/Final data normalized.xlsx",sheet=2) 
BM <- BM[BM$TIME!="24 hr",]
ID <- BM$ID
colnames(BM)

sum(BM[,47])/90


sqrt(var(BM[,47]))
sqrt(var(D[,118]))


Yj <- scale(BM[,47])
Yj<-scale(BM[,47],center=F,scale=F)
Yj<-scale(D[,118],center=F,scale=F)


plot(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),type="n",ylab="Raw rescaled data",ylim=c(0,2),xlab="Time (hrs) since breakfast",main="Valine",cex.main=cx)
for(i in 1:10){
  lines(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),Yj[ID==i],col=1+i)
}
lines(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),bmflmeans[38,],lwd=2,col="black")

plot(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),ylim=range(meandif.panels),xlab="Time (hrs) since breakfast", ylab="Difference in mean profiles",main="Metabolon mean - BMFL mean",type="n",cex.main=cx)
for(k in 1:length(BMFL.id)){
  lines(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),meandif.panels[k,],pch=16,col=ifelse((correl[k])>0.6,"navy","red"))
}

plot(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),type="n",ylim=c(-2,2),xlab="Time (hrs) since breakfast",main="Metabolon",ylab="Mean profile",cex.main=cx)
for(k in 1:length(Metabolon.id)){
  lines(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),metabmeans[k,],col=k)
}

plot(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),type="n",ylim=c(-2,2),xlab="Time (hrs) since breakfast",main="BMFL",ylab="Mean profile",cex.main=cx)
for(k in 1:length(BMFL.id)){
  lines(c(-0.5,0.5,1.5,3,4.5,6,9.5,11,12),bmflmeans[k,],col=k)
}
#dev.off()

# Maybe I should check influence of age, bmi on final results