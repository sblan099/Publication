library(asreml)

################################### FIELD DATA  #########################

dat<-read.table(file="Mixed Model Correlation Data.600.csv",sep=",",header=T)

#Prepare Data

summary(dat)
dat$SITE<-factor(dat$Site)
dat$YEAR<-factor(dat$Year)
dat$SEX<-factor(dat$Sex)
dat$RECAP<-factor(dat$Recapture)
dat$ID<-factor(dat$Code) 
dat$trait.1<-scale(dat$Aggression)
dat$trait.2<-scale(dat$Bin.Shell)
dat$trait.3<-scale(log(1+dat$Start))
dat$trait.4<-scale(dat$Move) 
#Make Model

MODEL<-asreml (cbind(trait.1,trait.2,trait.3,trait.4 )~trait+at(trait):SEX+at(trait):Time ,
               random=~corgh(trait):ID+idh(trait):SITE, 
               residual=~units:corgh(trait), 
               data=dat,na.action=na.method(x="include")) 
#Plot Model

summary(MODEL) 
#$varcomp
#                                                component  std.error     z.ratio bound %ch
#trait:SITE!trait_trait.1                      0.044183719 0.02747224  1.60830394     P 0.5
#trait:SITE!trait_trait.2                      0.103504849 0.04453042  2.32436251     P 0.5
#trait:SITE!trait_trait.3                      0.160313575 0.06176607  2.59549593     P 0.3
#trait:SITE!trait_trait.4                      0.162541222 0.06797533  2.39117967     P 0.5
#trait:ID!trait!trait.2:!trait!trait.1.cor    -0.513767527 0.07798190 -6.58829207     U 0.0
#trait:ID!trait!trait.3:!trait!trait.1.cor    -0.851124398 0.83870110 -1.01481254     U 0.4
#trait:ID!trait!trait.3:!trait!trait.2.cor     0.631097875 0.53779173  1.17349867     U 0.2
#trait:ID!trait!trait.4:!trait!trait.1.cor    -0.009805403 0.13024041 -0.07528694     U 0.5
#trait:ID!trait!trait.4:!trait!trait.2.cor    -0.038956429 0.15574327 -0.25013235     U 0.4
#trait:ID!trait!trait.4:!trait!trait.3.cor    -0.352064332 0.57074441 -0.61685113     U 0.0
#trait:ID!trait_trait.1                        0.765500325 0.08181050  9.35699332     P 0.0
#trait:ID!trait_trait.2                        0.593409229 0.09217222  6.43804868     P 0.0
#trait:ID!trait_trait.3                        0.081069508 0.15634049  0.51854455     P 0.3
#trait:ID!trait_trait.4                        0.427449356 0.13755218  3.10754321     P 0.0
#units:trait!R                                 1.000000000         NA          NA     F 0.0
#units:trait!trait!trait.2:!trait!trait.1.cor -0.260701576 0.18087207 -1.44135895     U 0.1
#units:trait!trait!trait.3:!trait!trait.1.cor -0.050295147 0.19406501 -0.25916648     U 1.4
#units:trait!trait!trait.3:!trait!trait.2.cor  0.369474031 0.16211891  2.27903109     U 0.0
#units:trait!trait!trait.4:!trait!trait.1.cor  0.182694566 0.19374762  0.94295129     U 0.3
#units:trait!trait!trait.4:!trait!trait.2.cor -0.038102992 0.19267921 -0.19775352     U 0.3
#units:trait!trait!trait.4:!trait!trait.3.cor -0.260125453 0.17492437 -1.48707382     U 0.1
#units:trait!trait_trait.1                     0.203528333 0.05389888  3.77611406     P 0.1
#units:trait!trait_trait.2                     0.277230293 0.07461029  3.71571137     P 0.0
#units:trait!trait_trait.3                     0.728462177 0.16162707  4.50705560     P 0.1
#units:trait!trait_trait.4                     0.453922482 0.12816623  3.54166987     P 0.0


wald.asreml(MODEL,denDF = "default",ssType = "conditional")# 
#                        Df denDF  F.inc  F.con Margin      Pr 
#trait                    4  38.5 0.5763 0.5763        0.68154
#at(trait, trait.1):SEX   2 408.8 0.2325 0.2842      B 0.75273
#at(trait, trait.2):SEX   2 375.4 1.9230 3.0620      B 0.04796
#at(trait, trait.3):SEX   2 456.0 3.4430 4.2710      B 0.01453
#at(trait, trait.4):SEX   2 387.7 2.1870 2.2110      B 0.11097
#at(trait, trait.1):Time  1 370.9 1.2130 0.0695      B 0.79217
#at(trait, trait.2):Time  1 475.6 5.8890 4.5250      B 0.03391
#at(trait, trait.3):Time  1 471.0 0.2741 0.2022      B 0.65318
#at(trait, trait.4):Time  1 452.2 0.0535 0.0535      B 0.81721#

################################### PRODUCE HEAT MAP  #########################
#Set Up Variables

N.traits<-4 
N.cor<-(N.traits*(N.traits-1))/2 
#Create Among-Individual Correlation Matrix and Standard Error Matrix

R.IND<-summary(MODEL)$varcomp$component[4+(1:N.cor)] 
M.IND<-matrix(nrow=N.traits,ncol=N.traits,1)
M.IND[upper.tri(M.IND)]<-R.IND
M.IND2<-t(M.IND)
M.IND[lower.tri(M.IND, diag=F)==T]<-M.IND2[lower.tri(M.IND2, diag=F)==T] 
#
R.IND.se<-summary(MODEL)$varcomp$std.error[4+(1:N.cor)] 
M.IND.se<-matrix(nrow=N.traits,ncol=N.traits,1)
M.IND.se[upper.tri(M.IND.se)]<-R.IND.se
M.IND2.se<-t(M.IND.se)
M.IND.se[lower.tri(M.IND.se, diag=F)==T]<-M.IND2.se[lower.tri(M.IND2.se, diag=F)==T] 
#
#Create Residual Correlation Matrix and Standard Error Matrix

R.RES<-summary(MODEL)$varcomp$component[16:21] 
M.RES<-matrix(nrow=N.traits,ncol=N.traits,1)
M.RES[upper.tri(M.RES)]<-R.RES
M.RES2<-t(M.RES)
M.RES[lower.tri(M.RES, diag=F)==T]<-M.RES2[lower.tri(M.RES2, diag=F)==T] 
#
R.RES.se<-summary(MODEL)$varcomp$std.error[16:21] 
M.RES.se<-matrix(nrow=N.traits,ncol=N.traits,1)
M.RES.se[upper.tri(M.RES.se)]<-R.RES.se
M.RES2.se<-t(M.RES.se)
M.RES.se[lower.tri(M.RES.se, diag=F)==T]<-M.RES2.se[lower.tri(M.RES2.se, diag=F)==T] 
#Set Up Color Palette and Labels

library(raster) 
rowcol <- expand.grid(row = seq(1, N.traits, 1), col = seq(N.traits, 1, -1)) 
colin <- colorRampPalette(colors = c("red", "white", "blue"), bias = 1)(64) 
#breaksin <- seq(min(M.to.plot),  max(M.to.plot), 0.01)
#zlimin <- c(-max(M.to.plot),  max(M.to.plot))
breaksin <- seq(-1,  1, 0.01)
zlimin <- c(-1,  1)

LABELS<-c("Handling reaction",
          "Latency shell emergence",
          "Latency initial movement",
          "Time spet moving") 
#Generate Among-Individual Correlation Heatmap

M.IND.no.diag<-M.IND
diag(M.IND.no.diag)<-NA
diag(M.IND.se)<-NA
#
rowcol2<-rowcol
rowcol2$sum<-rowcol$row+rowcol$col
rowcol2<-rowcol2[-which(rowcol2$sum==6),]
par(mfrow=c(1,1),mar=c(8,9,1,1))
plot(rasterFromXYZ(cbind(rowcol, z = c(M.IND))),
     col=colin,lab.breaks=breaksin,zlim=zlimin,
     axes=FALSE,legend=FALSE,box=F,main="among-individual correlations") 
text(x=rowcol$row, y=rowcol$col+0.25, labels=c(round(M.IND.no.diag,digits=3)),cex=0.85)
text(x=rowcol2$row, y=rowcol2$col     , labels="?",cex=0.85)
text(x=rowcol$row, y=rowcol$col-0.25, labels=c(round(M.IND.se,digits=3)),cex=0.85)
axis(1,at=seq(1,N.traits),    LABELS, tick=F,line=-3,las=2)
axis(2,at=seq(1,N.traits),rev(LABELS),tick=F,line=-1,las=1)
#text(1.33,2.3,"**", col=1,cex=1.2)
text(1.33,3.3,"**", col=1,cex=1.2)
#text(2.33,2.3,"**", col=1,cex=1.2)

#Generate Within-Individual Correlation Heatmap

#the residual correlations
M.RES.no.diag<-M.RES
diag(M.RES.no.diag)<-NA
diag(M.RES.se)<-NA
#
plot(rasterFromXYZ(cbind(rowcol, z = c(M.RES))),
     col=colin,lab.breaks=breaksin,zlim=zlimin,
     axes=FALSE,legend=FALSE,box=F,main="within-individual correlations")
#Annotations and Legends

text(x=rowcol$row, y=rowcol$col+0.25, labels=c(round(M.RES.no.diag,digits=3)),cex=0.85)
text(x=rowcol2$row, y=rowcol2$col     , labels="?",cex=0.85)
text(x=rowcol$row, y=rowcol$col-0.25, labels=c(round(M.RES.se,digits=3)),cex=0.85)
axis(1,at=seq(1,N.traits),    LABELS, tick=F,line=-3,las=2)
axis(2,at=seq(1,N.traits),rev(LABELS),tick=F,line=-1,las=1)
text(2.33,2.3,"**", col=1,cex=1.2)

######################################## LAB DATA ############################# 

dat<-read.table(file="Combined Repeatability Dataset.600.csv",sep=",",header=T)
#Prepare Data

summary(dat)
dat$SITE<-factor(dat$Site)
dat$YEAR<-factor(dat$Year)
dat$SEX<-factor(dat$Sex)
dat$DAY<-factor(dat$Day)
dat$ID<-factor(dat$Code)
dat$trait.1<-scale(dat$Aggression)
dat$trait.2<-scale(dat$BIN.Shell)
dat$trait.3<-scale(log(1+dat$Start))
dat$trait.4<-scale(dat$Move)
hist(dat$trait.3)
#Make Model

MODEL<-asreml(cbind(trait.1,trait.2,trait.3,trait.4)~trait+at(trait):SEX+at(trait):DAY,
              random=~corgh(trait):ID+idh(trait):SITE,
              residual=~units:corgh(trait),
              data=dat,na.action=na.method(x="include"))
#Plot Model

summary(MODEL)
#$varcomp
#                                                 component  std.error    z.ratio bound %ch
#trait:SITE!trait_trait.1                      6.739390e-02 0.06535546  1.0311900     P 0.0
#trait:SITE!trait_trait.2                      5.049882e-08         NA         NA     B  NA
#trait:SITE!trait_trait.3                      2.989333e-01 0.17185670  1.7394337     P 0.1
#trait:SITE!trait_trait.4                      1.374939e-07         NA         NA     B 0.0
#trait:ID!trait!trait.2:!trait!trait.1.cor    -5.843874e-01 0.10747390 -5.4374821     U 0.0
#trait:ID!trait!trait.3:!trait!trait.1.cor    -6.246128e-01 0.11940660 -5.2309737     U 0.1
#trait:ID!trait!trait.3:!trait!trait.2.cor     6.295772e-01 0.12154789  5.1796637     U 0.1
#trait:ID!trait!trait.4:!trait!trait.1.cor     6.098126e-01 0.12678928  4.8096545     U 0.0
#trait:ID!trait!trait.4:!trait!trait.2.cor    -4.779230e-01 0.15072665 -3.1707927     U 0.0
#trait:ID!trait!trait.4:!trait!trait.3.cor    -6.057586e-01 0.14006985 -4.3246897     U 0.0
#trait:ID!trait_trait.1                        5.520284e-01 0.12693375  4.3489485     P 0.0
#trait:ID!trait_trait.2                        6.975085e-01 0.14860619  4.6936705     P 0.0
#trait:ID!trait_trait.3                        6.362560e-01 0.17321324  3.6732525     P 0.1
#trait:ID!trait_trait.4                        6.710113e-01 0.17694088  3.7922907     P 0.0
#units:trait!R                                 1.000000e+00         NA         NA     F 0.0
#units:trait!trait!trait.2:!trait!trait.1.cor -9.604931e-02 0.08034875 -1.1954051     U 0.0
#units:trait!trait!trait.3:!trait!trait.1.cor -1.132253e-01 0.08739106 -1.2956161     U 0.0
#units:trait!trait!trait.3:!trait!trait.2.cor  3.401568e-01 0.08462386  4.0196326     U 0.0
#units:trait!trait!trait.4:!trait!trait.1.cor -1.196354e-01 0.08800091 -1.3594787     U 0.0
#units:trait!trait!trait.4:!trait!trait.2.cor  6.098250e-02 0.09579126  0.6366187     U 0.0
#units:trait!trait!trait.4:!trait!trait.3.cor  2.668853e-01 0.08651429  3.0848692     U 0.0
#units:trait!trait_trait.1                     2.279373e-01 0.02607241  8.7424711     P 0.0
#units:trait!trait_trait.2                     3.136403e-01 0.03598788  8.7151645     P 0.0
#units:trait!trait_trait.3                     3.628706e-01 0.04769230  7.6085782     P 0.0
#units:trait!trait_trait.4                     4.908437e-01 0.06443708  7.6174111     P 0.0


wald.asreml(MODEL,denDF = "default",ssType = "conditional")
#                       Df denDF   F.inc   F.con Margin      Pr
#trait                   4  27.8  0.9375  0.9375        0.45676
#at(trait, trait.1):SEX  1  53.1 13.6100 10.3000      B 0.00226
#at(trait, trait.2):SEX  1  55.9  1.9530  1.0380      B 0.31269
#at(trait, trait.3):SEX  1  43.8  0.9932  0.3159      B 0.57695
#at(trait, trait.4):SEX  1  46.6  0.5791  0.6321      B 0.43060
#at(trait, trait.1):DAY  3 153.7  8.7550  8.8600      B 0.00002
#at(trait, trait.2):DAY  3 153.3  0.0525  0.6275      B 0.59832
#at(trait, trait.3):DAY  3 124.2  6.8610  4.5460      B 0.00465
#at(trait, trait.4):DAY  3 121.8  3.4070  3.4070      B 0.01985
#

################################### PRODUCE HEAT MAP  #########################
#Set Up Variables


N.traits<-4
N.cor<-(N.traits*(N.traits-1))/2
#Create Among-Individual Correlation Matrix and Standard Error Matrix

R.IND<-summary(MODEL)$varcomp$component[4+(1:N.cor)]
M.IND<-matrix(nrow=N.traits,ncol=N.traits,1)
M.IND[upper.tri(M.IND)]<-R.IND
M.IND2<-t(M.IND)
M.IND[lower.tri(M.IND, diag=F)==T]<-M.IND2[lower.tri(M.IND2, diag=F)==T]
#
R.IND.se<-summary(MODEL)$varcomp$std.error[4+(1:N.cor)]
M.IND.se<-matrix(nrow=N.traits,ncol=N.traits,1)
M.IND.se[upper.tri(M.IND.se)]<-R.IND.se
M.IND2.se<-t(M.IND.se)
M.IND.se[lower.tri(M.IND.se, diag=F)==T]<-M.IND2.se[lower.tri(M.IND2.se, diag=F)==T]
#
#Create Residual Correlation Matrix and Standard Error Matrix

R.RES<-summary(MODEL)$varcomp$component[16:21]
M.RES<-matrix(nrow=N.traits,ncol=N.traits,1)
M.RES[upper.tri(M.RES)]<-R.RES
M.RES2<-t(M.RES)
M.RES[lower.tri(M.RES, diag=F)==T]<-M.RES2[lower.tri(M.RES2, diag=F)==T]
#
R.RES.se<-summary(MODEL)$varcomp$std.error[16:21]
M.RES.se<-matrix(nrow=N.traits,ncol=N.traits,1)
M.RES.se[upper.tri(M.RES.se)]<-R.RES.se
M.RES2.se<-t(M.RES.se)
M.RES.se[lower.tri(M.RES.se, diag=F)==T]<-M.RES2.se[lower.tri(M.RES2.se, diag=F)==T]
#Set Up Color Palette and Labels

library(raster)
rowcol <- expand.grid(row = seq(1, N.traits, 1), col = seq(N.traits, 1, -1))
colin <- colorRampPalette(colors = c("red", "white", "blue"), bias = 1)(64)
#breaksin <- seq(min(M.to.plot),  max(M.to.plot), 0.01)
#zlimin <- c(-max(M.to.plot),  max(M.to.plot))
breaksin <- seq(-1,  1, 0.01)
zlimin <- c(-1,  1)

LABELS<-c("Handling reaction",
          "Latency shell emergence",
          "Latency initial movement",
          "Time spet moving")
#Generate Among-Individual Correlation Heatmap

M.IND.no.diag<-M.IND
diag(M.IND.no.diag)<-NA
diag(M.IND.se)<-NA
#
rowcol2<-rowcol
rowcol2$sum<-rowcol$row+rowcol$col
rowcol2<-rowcol2[-which(rowcol2$sum==6),]

par(mfrow=c(1,1),mar=c(8,9,1,1))
plot(rasterFromXYZ(cbind(rowcol, z = c(M.IND))),
     col=colin,lab.breaks=breaksin,zlim=zlimin,
     axes=FALSE,legend=FALSE,box=F,main="among-individual correlations")
text(x=rowcol$row, y=rowcol$col+0.25, labels=c(round(M.IND.no.diag,digits=3)),cex=0.85)
text(x=rowcol2$row, y=rowcol2$col     , labels="?",cex=0.85)
text(x=rowcol$row, y=rowcol$col-0.25, labels=c(round(M.IND.se,digits=3)),cex=0.85)
axis(1,at=seq(1,N.traits),    LABELS, tick=F,line=-3,las=2)
axis(2,at=seq(1,N.traits),rev(LABELS),tick=F,line=-1,las=1)
text(1.33,2.3,"**", col=1,cex=1.2)
text(1.33,3.3,"**", col=1,cex=1.2)
text(2.33,2.3,"**", col=1,cex=1.2)
text(1.33,1.3,"**", col=1,cex=1.2)
text(2.33,1.3,"**", col=1,cex=1.2)
text(3.33,1.3,"**", col=1,cex=1.2)
#Generate Within-Individual Correlation Heatmap

#the residual correlations
M.RES.no.diag<-M.RES
diag(M.RES.no.diag)<-NA
diag(M.RES.se)<-NA
#
plot(rasterFromXYZ(cbind(rowcol, z = c(M.RES))),
     col=colin,lab.breaks=breaksin,zlim=zlimin,
     axes=FALSE,legend=FALSE,box=F,main="within-individual correlations")
#Annotations and Legends

text(x=rowcol$row, y=rowcol$col+0.25, labels=c(round(M.RES.no.diag,digits=3)),cex=0.85)
text(x=rowcol2$row, y=rowcol2$col     , labels="?",cex=0.85)
text(x=rowcol$row, y=rowcol$col-0.25, labels=c(round(M.RES.se,digits=3)),cex=0.85)
axis(1,at=seq(1,N.traits),    LABELS, tick=F,line=-3,las=2)
axis(2,at=seq(1,N.traits),rev(LABELS),tick=F,line=-1,las=1)
text(2.33,2.3,"**", col=1,cex=1.2)
text(3.33,1.3,"**", col=1,cex=1.2)
