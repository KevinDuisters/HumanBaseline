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
# p value analysis
#---------------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------------#
# load pvallist, utlist and pvec from analysis (saved as Rdata)
load("data/pvals.Rdata")

#---------------------------------------------------------------------------------------------------------#
# Unpack and reformat p-values to prepare for visualization
qthresh <- 0.1
signiflist <- signifnamelist <- signifcountlist <- lapply(1:3,function(x)NULL)
for(k in 1:3){
  
  ut <- as.character(utlist[[k]])
  ni <- length(ut)
  signifmat <- matrix(0,ni,ni,dimnames=list(ut,ut))
  signifnames <- logical(pvec[k])
  signifnamecount <- numeric(pvec[k])
  
  
  for(r in 1:(ni-1)){
    for(c in (r+1):ni){
      # Bonf within, BH across
      pvals <- apply(t(1:pvec[k]),2,function(j)pvallist[[k]][[j]][r,c])
      qvals <- p.adjust((ni-1)*ni/2*pvals,"BH")  # Bonf within, BH across
      for(j in 1:pvec[k]){
        if(qvals[j] < qthresh){ 
          signifmat[r,c] <- signifmat[r,c] + 1
          signifnames[j] <- T
          signifnamecount[j] <- signifnamecount[j] + 1
        }
      }
      signifmat[c,r] <- signifmat[r,c]
    }
  }
  signiflist[[k]] <- signifmat
  signifnamelist[[k]] <- namelist[[k]][signifnames]
  signifcountlist[[k]] <- signifnamecount
}

# Manual adjust to ordering of ut
signiflist[[1]] # ok
signiflist[[2]] # ok

signiflist[[3]] # not ok
out <- aux <- signiflist[[3]]
out[,1:7] <- aux[,c(4,3,2,7,5,1,6)]
aux <- out
out[1:7,] <- aux[c(4,3,2,7,5,1,6),]
colnames(out) <- rownames(out) <- c("-0.5 hr","0.5 hr","1.5 hr","3 hr","6 hr","9.5 hr","12 hr")
signiflist[[3]] <- out
rm(out,aux)
signiflist[[3]] #ok now

#-----------------------------------------------------------------------------------------------------------
# Make visual
library(RColorBrewer)
histnames <- c("A: Plasma Metabolites","B: Urine Amines","C: Plasma Proteins")
cols <- c("white",brewer.pal(9,"Blues")[-1],"navy") 
cx <- 3

pdf("Figures/Contrasts.pdf",width=20)
par(mfrow=c(1,3),mar=c(6,8,6,2)+0.1)

for(k in 1:3){
  sm <- signiflist[[k]]
  ni <- ncol(sm)
  par(xpd=T)
  plot(0:(nrow(sm)),0:(ncol(sm)),type="n",axes=F,xlab="",ylab="",main="")
  for(r in 1:(ni-1)){
    for(c in (r+1):(ni)){
      polygon(x=c(-1.5+c,-0.5+c,-0.5+c,-1.5+c,-1.5+c),y=c(-0.5+r,-0.5+r,0.5+r,0.5+r,-0.5+r),col=
                cols[sum( (0:10)/10 <= 3.33*sm[r,c]/pvec[k])]) # multiplied by 3.33: range up to 30%
    }}
  short.hrs <- sapply(colnames(sm),function(s)substr(s,1,nchar(s)-3))
  title(paste0(histnames[k]," [",pvec[k],"]"),cex.main=4,line=3)
  axis(side=1,short.hrs,at=0:(ni-1),tick=F,cex.axis=cx,line=-1,las=2)
  axis(side=2,short.hrs,at=(1:ni),tick=F,las=1,line=1,cex.axis=cx)
  
  
}

dev.off()

#-----------------------------------------------------------------------------------------------------------
# custom color legend
pdf("Figures/Contrastslegend.pdf",width=20)
par(mfrow=c(1,1))
plot.new()
yc <- 0.4
xc <- 0.2
for(j in 1:length(cols)){
  rect(xc+0.05*(j-1),yc,xc+0.05*j,yc+0.1,col=cols[j])
}
text(x=0.45,y=0.6,"% significant biochemicals in panel",cex=cx)
text(x=xc-0.05,y=yc+0.05,"0 %",cex=cx)
text(x=xc+length(cols)*0.05+0.07,y=yc+0.05,"30 %",cex=cx)
dev.off()

