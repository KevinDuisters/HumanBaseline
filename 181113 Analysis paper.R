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
library(writexl)
library(lsmeans)
library(nlme)

# General info
totaltimes <- c("-2 hr","-0.5 hr","0.5 hr","1.5 hr","3 hr","4 hr","4.5 hr","6 hr","8 hr","9.5 hr","11 hr","12 hr") # excl 24 hr
numtimes  <- c(-2,-0.5,0.5,1.5,3,4,4.5,6,8,9.5,11,12)
lc <- lmeControl(maxIter=1e3,msMaxIter=1e3,msMaxEval=1e3,niterEM=100,returnObject=T,rel.tol=1e-10,xf.tol=1e-18)

# Create result matrices
pvec <- numeric(3)
namelist <- utlist <- coeflist <- mlist <- slist <- s2taulist <- s2blist <- s2wlist <- choicelist <- pvallist <- errorlist <- lapply(1:3,function(x)NULL)

# loop over entities
for(k in 1:3){
  
  # localize excel sheet of paneldata
  D <- as.data.frame(read_xlsx("data/data.xlsx",sheet=k) )
  
  # short version for testing
  #D <- as.data.frame(read_xlsx("data/data.xlsx",sheet=k) )[,1:20]
  
  
  
  p <- ncol(D)-4
  pvec[k] <- p
  
  # exclude 24 hr measurement for simplicity (no day reproduction)
  D <- D[D$TIME!="24 hr",]
  
  ID <- D$ID
  AGE <- D$Age
  BMI <- D$BMI
  TIME <- factor(D$TIME,levels=totaltimes)
  ut <- unique(TIME)  # levels that actually are present. ordering influences p value columns!!! Manual adjust
  utlist[[k]] <- ut
  namelist[[k]] <- names(D)[-(1:4)]
  ni <- length(unique(D$TIME))
  X<-cbind(kronecker(rep(1,10),diag(ni)),AGE,BMI)
  
  # create panel results
  errorlist[[k]] <- matrix(0,p,6,dimnames=list(names(D)[-c(1:4)],paste0("M",1:6)))
  coeflist[[k]] <- matrix(NA,p,2+length(unique(TIME)),dimnames=list(row=names(D)[-c(1:4)],col=c(levels(TIME)[sort(unique(TIME))],"AGE","BMI")))
  mlist[[k]] <- slist[[k]] <- matrix(NA,p,length(unique(TIME)),dimnames=list(row=names(D)[-c(1:4)],col=c(levels(TIME)[sort(unique(TIME))])))
  s2taulist[[k]] <- s2blist[[k]] <- s2wlist[[k]] <- choicelist[[k]] <- numeric(p)
  pvallist[[k]] <- lapply(1:p,function(x) NULL)
  names(pvallist[[k]]) <- namelist[[k]]
  
  options(warn=2) # upgrade convergence warning to error inside lme
  
  for(j in 1:p){
    
    # normalize
    Yj <- scale(D[,4+j])
    
    # Step 1: covariance matrix selection
    M1 <- tryCatch(lme(Yj~-1+TIME+AGE + BMI,random=~1|ID,control=lc,method="REML"),error=function(e){3.14})  # simple structure (diag)
    M2 <- tryCatch(lme(Yj~-1+TIME+AGE + BMI,random=~1|ID,correlation=corCAR1(0.1,~TIME|ID),control=lc,method="REML"),error=function(e){3.14}) # CAR1
    M3 <- tryCatch(lme(Yj~-1+TIME+AGE + BMI,random=~1|ID,correlation=corGaus(form=~TIME|ID),control=lc,method="REML"),error=function(e){3.14})
    M4 <- tryCatch(lme(Yj~-1+TIME+AGE + BMI,random=~1|ID,weights=varIdent(form=~1|TIME),control=lc,method="REML"),error=function(e){3.14}) # diag, random int, heterosked error variance per time 
    M5 <- tryCatch(lme(Yj~-1+TIME+AGE + BMI,random=~1|ID,correlation=corCAR1(0.1,~TIME|ID),weights=varIdent(form=~1|TIME),control=lc,method="REML"),error=function(e){3.14})
    M6 <- tryCatch(lme(Yj~-1+TIME+AGE + BMI,random=~1|ID,correlation=corGaus(form=~TIME|ID),weights=varIdent(form=~1|TIME),control=lc,method="REML"),error=function(e){3.14})
   
    
    if(mode(M2)!="list"){errorlist[[k]][j,2] <- 1; M2 <- M1}
    if(mode(M3)!="list"){errorlist[[k]][j,3] <- 1; M3 <- M1}
    if(mode(M4)!="list"){errorlist[[k]][j,4] <- 1; M4 <- M1}
    if(mode(M5)!="list"){errorlist[[k]][j,5] <- 1; M5 <- M4}
    if(mode(M6)!="list"){errorlist[[k]][j,6] <- 1; M6 <- M4}
    Mlist<-list(M1,M2,M3,M4,M5,M6)  
    
    
    # while loop to exlude choices that deliver ML issues (I want model to work both for REML and ML)
    MLproblem <- T
    while(MLproblem==T){
      
      AICs <- sapply(Mlist,function(x)summary(x)$AIC) # model selection based on AIC
      choice <- which.min(AICs)
      full.model <-Mlist[[choice]]
      choicelist[[k]][j] <- choice
      
      
      # Step 2: get coefs and variance components in optimal (full) model
      coefsj <- summary(full.model)$tTable[,1]
      coeflist[[k]][j,] <- coefsj
      s2taulist[[k]][j] <- var(coefsj[-which(names(coefsj)%in%c("AGE","BMI"))])
      s2blist[[k]][j] <- getVarCov(full.model)[[1]]
      Vj <- getVarCov(full.model)[[1]] + getVarCov(full.model,type="conditional")[[1]]
      trSigmaw <- diag(getVarCov(full.model,type="conditional")[[1]])
      s2wlist[[k]][j] <- mean(trSigmaw)
      
      
      mlist[[k]][j,] <- coefsj[-which(names(coefsj)%in%c("AGE","BMI"))]+mean(AGE)*coefsj["AGE"]+mean(BMI)*coefsj["BMI"]
      slist[[k]][j,] <- sqrt(diag(Vj))
      
      # Step 3: LRT tests on 'meal effects'
      # FULL LOOP OVER LRT TESTS ON ALL PAIRWISE COMPARISONS (SLOW BUT CORRECT)
      pvallist[[k]][[j]] <- matrix(NA,ni,ni,dimnames=list(ut,ut))
      counter <- counter1 <- 0
      ll.full <- tryCatch(as.numeric(summary(update(full.model,method="ML"))$logLik),error=function(e)3.14) # REML = FALSE for fixed effect comparison
      for(r in 1:(ni-1)){
        for(c in (r+1):ni){
          if(counter==counter1){
            counter1 <- counter1 + 1  
            BTIME <- TIME
            BTIME[TIME==ut[r]] <- ut[c]
            ll.reduced <- tryCatch(as.numeric(summary(update(full.model,fixed=Yj~-1+BTIME+AGE+BMI,method="ML"))$logLik),error=function(e)3.14) # REML = FALSE for fixed effect comparison
            if(ll.full == 3.14 | ll.reduced==3.14){
              Mlist[[choice]] <- M1 # reset erroneous models to default (M1)
              errorlist[[k]][j,choice] <- 1
            }else{
              pvallist[[k]][[j]][r,c] <- (1-pchisq(2*(ll.full-ll.reduced),df=1)) # reject for signif contrast
              counter <- counter+1
            }
          }
        }
      }
      if(counter==counter1){MLproblem<-FALSE} # all went well
    } # end while
    
    
    rm(M1,M2,M3,M4,M5,M6,Mlist) # clear memory
    cat("\r",paste("k=",k," ,j=",j))
  }
  
  options(warn=0) # reset warnings
} # end panel loop



#-----------------------------------------------------------------------------------------------------------
# Store data and generate Supplemental Information result tables 
save.now<-F
if(save.now==T){

  # create file for pvalue analysis (p-value figure)
  save(pvallist,pvec,utlist,namelist,file="data/pvals.Rdata")

  # Create ST1 excel file (Supplemental Information) for Figures and Tables
  sigmab <- sqrt(unlist(s2blist))
  sigmaeps <- sqrt(unlist(s2wlist))
  sigmatime <- sqrt(unlist(s2taulist))
  fullnames <- unlist(namelist)
  panels <- c(rep("Metabolon",pvec[1]),rep("BMFL",pvec[2]),rep("SomaLogic",pvec[3]))
  samples <- c(rep("plasma",pvec[1]),rep("urine",pvec[2]),rep("plasma",pvec[3]))

  # Supplemental table of coefficients
  coeftotmat <- matrix(NA,nrow=sum(pvec),ncol=2+length(totaltimes),dimnames=list(fullnames,c("age","bmi",totaltimes)))
  coeftotmat[1:pvec[1],c(1,2,4,5,6,7,9,10,12,13,14)] <- unlist(coeflist[[1]])[,c(10,11,1:9)]
  coeftotmat[(cumsum(pvec)[1]+1):(cumsum(pvec)[2]),c(1,2,3,8,11,14)] <- unlist(coeflist[[2]])[,c(5,6,1:4)]
  coeftotmat[(cumsum(pvec)[2]+1):(cumsum(pvec)[3]),c(1,2,4,5,6,7,10,12,14)] <- unlist(coeflist[[3]])[,c(8,9,1:7)]

  stotmat <- matrix(NA,nrow=sum(pvec),ncol=length(totaltimes),dimnames=list(fullnames,c(totaltimes)))
  stotmat[1:pvec[1],c(2,3,4,5,7,8,10,11,12)] <- unlist(slist[[1]])
  stotmat[(cumsum(pvec)[1]+1):(cumsum(pvec)[2]),c(1,6,9,12)] <- unlist(slist[[2]])
  stotmat[(cumsum(pvec)[2]+1):(cumsum(pvec)[3]),c(2,3,4,5,8,10,12)] <- unlist(slist[[3]])


  ST1 <- data.frame(fullnames,panels,samples,coeftotmat,sigmab,sigmatime,sigmaeps,stotmat)
  colnames(ST1) <- c("compound","panel","sample","Age","BMI",totaltimes,"sigma_b","sigma_time","sigma_eps",paste("sigmat",totaltimes))
  #write.table(ST1,file="data/ST1.txt",na=" ",row.names=F,sep=";",dec=",",qmethod="double")
  write_xlsx(ST1,"data/ST1.xlsx")

}


#-----------------------------------------------------------------------------------------------------------
