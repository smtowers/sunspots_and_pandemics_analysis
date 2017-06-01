##################################################################################
##################################################################################
##################################################################################
##################################################################################
set.seed(183672)
mybinom=function(x,n,p){
  a = pbinom(x,n,p,lower.tail=F)+dbinom(x,n,p)
  return(a)
}

require("sfsmisc")
require("kSamples")

##################################################################################
##################################################################################
# first run make_table.R to collate sources
##################################################################################
lpandemic = 1
if (lpandemic){
 pdat = read.table("summary_pandemic_data.txt",header=T,as.is=T)
}else{
 pdat = read.table("summary_serious_outbreak_data.txt",header=T,as.is=T)
}

nagree = seq(1,max(pdat$nreview_agree))
cat(nagree,"\n")
vp_ks_phi = numeric(0)
vp_ad_phi = numeric(0)
vp_ks_Q = numeric(0)
vp_ad_Q = numeric(0)
vp_ks_sunspot = numeric(0)
vp_ad_sunspot = numeric(0)
vagree = numeric(0)
vtype = numeric(0)
vpbinom_max = numeric(0)
vpbinom_extrema = numeric(0)
for (itype in 1:2){
for (iagree in nagree){

vtype = append(vtype,itype)
vagree = append(vagree,iagree)
vyear = pdat$year[which(pdat$nreview_agree>=iagree)]
#vyear[vyear==1800] = 1799

##################################################################################
##################################################################################
# read in the sunspot data
# note that Tapping phase isn't very well defined for dates from 2008
# onwards because we haven't yet gotten to the next minimum in solar
# activity (2008 was the year of a minimum)
##################################################################################
sdat = read.table("sunspot_wolf_and_group_1700_to_2014.txt",header=T,as.is=T)
if (itype==1){
  sdat$sunspot = sdat$wolf
  sdat$extrema = sdat$wolf_extrema
  #sdat$phi = sdat$tapping_phase_group
  #sdat$phi = sdat$ertel_Q_wolf
  sdat$phi = sdat$tapping_phase_wolf
  sdat$Q = sdat$ertel_Q_wolf
}else{
  sdat$sunspot = sdat$group
  sdat$extrema = sdat$group_extrema
  #sdat$phi = sdat$tapping_phase_group
  #sdat$phi = sdat$ertel_Q_group
  sdat$phi = sdat$tapping_phase_group
  sdat$Q = sdat$ertel_Q_group
}
sdat = subset(sdat,year<=2014)
sdat$near_extrema = rep(0,nrow(sdat))
sdat$near_max = rep(0,nrow(sdat))

l = which(sdat$extrema==(-1)|sdat$extrema==(+1))
sdat$near_extrema[l] = 1
sdat$near_extrema[l-1] = 1
sdat$near_extrema[pmin(l+1,nrow(sdat))] = 1

l = which(sdat$extrema==(+1))
sdat$near_max[l] = 1
sdat$near_max[l-1] = 1
sdat$near_max[pmin(l+1,nrow(sdat))] = 1

cat("min and max phi",min(sdat$phi),max(sdat$phi),"\n")

mult.fig(4)
m = which(sdat$year%in%vyear)
amain = "Pandemic years"
if (!lpandemic) amain = "Serious outbreak and pandemic years"
vbreaks = seq(-1.05,1.05,0.1)

##################################################################################
##################################################################################
# work out the p-value of the phi statistic
##################################################################################
a=ks.test(sdat$phi,sdat$phi[m])
vD = numeric(0)
for (iter in 1:1000){
  vD = append(vD,ks.test(sdat$phi,sdat$phi[sample(nrow(sdat),length(m),replace=T)])$statistic)
}
a$p.value = sum(vD>=a$statistic)/length(vD)

b=ad.test(sdat$phi,sdat$phi[m])
vD = numeric(0)
for (iter in 1:100){
  if (iter%%10==0) cat(iter,1000,"\n")
  vD = append(vD,ad.test(sdat$phi,sdat$phi[sample(nrow(sdat),length(m),replace=T)])$ad[1,1])
}
b$ad[1,3] = sum(vD>=b$ad[1,1])/length(vD)
vp_ks_phi = append(vp_ks_phi,a$p.value)
vp_ad_phi = append(vp_ad_phi,b$ad[1,3])
#print(a)
#print(b)

##################################################################################
##################################################################################
# work out the p-value of the Q statistic
##################################################################################
a=ks.test(sdat$Q,sdat$Q[m])
vD = numeric(0)
for (iter in 1:1000){
  vD = append(vD,ks.test(sdat$Q,sdat$Q[sample(nrow(sdat),length(m),replace=T)])$statistic)
}
a$p.value = sum(vD>=a$statistic)/length(vD)

b=ad.test(sdat$Q,sdat$Q[m])
vD = numeric(0)
for (iter in 1:100){
  if (iter%%10==0) cat(iter,1000,"\n")
  vD = append(vD,ad.test(sdat$Q,sdat$Q[sample(nrow(sdat),length(m),replace=T)])$ad[1,1])
}
b$ad[1,3] = sum(vD>=b$ad[1,1])/length(vD)
vp_ks_Q = append(vp_ks_Q,a$p.value)
vp_ad_Q = append(vp_ad_Q,b$ad[1,3])

##################################################################################
##################################################################################
# work out the p-value of the sunspot distribution comparison
##################################################################################
a=(ks.test(sdat$sunspot,sdat$sunspot[m]))
b=(ad.test(sdat$sunspot,sdat$sunspot[m]))
vp_ks_sunspot = append(vp_ks_sunspot,a$p.value)
vp_ad_sunspot = append(vp_ad_sunspot,b$ad[1,3])
print(a)
print(b)


##################################################################################
##################################################################################
# work out the Binomial p-values
##################################################################################
imax = which(sdat$near_max[m]==1)
imin = which(sdat$near_extrema[m]==1)
cat("Number of reference years = ",nrow(sdat),"\n")
cat("\n")
cat("******************************\n")
cat("******************************\n")
cat("Fraction of reference years at max = ",sum(sdat$near_max==1)/nrow(sdat),"\n")
h=(binom.test(length(imax),length(m),sum(sdat$near_max==1)/nrow(sdat)))
abinom = h$p.value
vpbinom_max = append(vpbinom_max,abinom)
print(h)
cat("Binomial probability:",abinom,"\n")

cat("\n")
cat("******************************\n")
cat("******************************\n")
cat("Fraction of reference years at min = ",sum(sdat$near_extrema==1)/nrow(sdat),"\n")
h=(binom.test(length(imin),length(m),sum(sdat$near_extrema==1)/nrow(sdat)))
abinom = h$p.value
vpbinom_extrema = append(vpbinom_extrema,abinom)
print(h)
cat("Binomial probability:",abinom,"\n")
} # loop over # reviewers agreeing
} # loop over Wolf and Group sunspots

##################################################################################
##################################################################################
# plot everything
##################################################################################

amain = "Pandemic years 1700 to 2014"
if (!lpandemic) amain = "Serious outbreak and pandemic years"
mult.fig(6,main=amain,mfrow=c(2,4),mar = c(3,3,3,1))
i = which(vtype==1)
amain = "Wolf sunspot \043's:\n Binomial probability"
plot(vagree[i],vpbinom_max[i],type="l",lwd=4,ylim=c(0,1),main=amain,ylab="p-value",xlab="Min \043 reviewers agreeing",yaxs="i",col=4,xaxs="i")
lines(vagree[i],vpbinom_extrema[i],lwd=4,col=3)
lines(c(-1,1000),c(0.05,0.05),lty=3)
legend("topleft",legend=c("\043 years near max","\043 years near extrema","","","","",""),col=c(4,3,0,0,0,0,0),lwd=4,bty="n")

amain = "Wolf sunspot \043's:\n Test of Q distribution"
plot(vagree[i],vp_ks_Q[i],type="l",lwd=4,ylim=c(0,1),main=amain,ylab="p-value",xlab="Min \043 reviewers agreeing",yaxs="i",xaxs="i")
lines(vagree[i],vp_ad_Q[i],lwd=4,col=2)
lines(c(-1,1000),c(0.05,0.05),lty=3)
legend("topleft",legend=c("K-S test","A-D test"),col=c(1,2),lwd=4,bty="n")

amain = "Wolf sunspot \043's:\n Test of phi distribution"
plot(vagree[i],vp_ks_phi[i],type="l",lwd=4,ylim=c(0,1),main=amain,ylab="p-value",xlab="Min \043 reviewers agreeing",yaxs="i",xaxs="i")
lines(vagree[i],vp_ad_phi[i],lwd=4,col=2)
lines(c(-1,1000),c(0.05,0.05),lty=3)
legend("topleft",legend=c("K-S test","A-D test"),col=c(1,2),lwd=4,bty="n")

amain = "Wolf sunspot \043's:\n Test of sunspot \043 distribution"
plot(vagree[i],vp_ks_sunspot[i],type="l",lwd=4,ylim=c(0,1),main=amain,ylab="p-value",xlab="Min \043 reviewers agreeing",yaxs="i",xaxs="i")
lines(vagree[i],vp_ad_sunspot[i],lwd=4,col=2)
lines(c(-1,1000),c(0.05,0.05),lty=3)
legend("topleft",legend=c("K-S test","A-D test"),col=c(1,2),lwd=4,bty="n")


i = which(vtype==2)
amain = "Group sunspot \043's:\n Binomial probability"
plot(vagree[i],vpbinom_max[i],type="l",lwd=4,ylim=c(0,1),main=amain,ylab="p-value",xlab="Min \043 reviewers agreeing",yaxs="i",col=4,xaxs="i")
lines(vagree[i],vpbinom_extrema[i],lwd=4,col=3)
lines(c(-1,1000),c(0.05,0.05),lty=3)
legend("topleft",legend=c("\043 years near max","\043 years near extrema","","","","",""),col=c(4,3,0,0,0,0,0),lwd=4,bty="n")

amain = "Group sunspot \043's:\n Test of Q distribution"
plot(vagree[i],vp_ks_Q[i],type="l",lwd=4,ylim=c(0,1),main=amain,ylab="p-value",xlab="Min \043 reviewers agreeing",yaxs="i",xaxs="i")
lines(vagree[i],vp_ad_Q[i],lwd=4,col=2)
lines(c(-1,1000),c(0.05,0.05),lty=3)
legend("topleft",legend=c("K-S test","A-D test"),col=c(1,2),lwd=4,bty="n")

amain = "Group sunspot \043's:\n Test of phi distribution"
plot(vagree[i],vp_ks_phi[i],type="l",lwd=4,ylim=c(0,1),main=amain,ylab="p-value",xlab="Min \043 reviewers agreeing",yaxs="i",xaxs="i")
lines(vagree[i],vp_ad_phi[i],lwd=4,col=2)
lines(c(-1,1000),c(0.05,0.05),lty=3)
legend("topleft",legend=c("K-S test","A-D test"),col=c(1,2),lwd=4,bty="n")

amain = "Group sunspot \043's:\n Test of sunspot \043 distribution"
plot(vagree[i],vp_ks_sunspot[i],type="l",lwd=4,ylim=c(0,1),main=amain,ylab="p-value",xlab="Min \043 reviewers agreeing",yaxs="i",xaxs="i")
lines(vagree[i],vp_ad_sunspot[i],lwd=4,col=2)
lines(c(-1,1000),c(0.05,0.05),lty=3)
legend("topleft",legend=c("K-S test","A-D test"),col=c(1,2),lwd=4,bty="n")


