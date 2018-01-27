library(VGAM)
rm(list=ls())
graphics.off()
theta<-seq(-pi,pi,0.0001)
col_leg<-NULL
par(mfrow=c(1,2))
R<-15
sigma<-15
for (I in seq(10,100,20)){
     c<-2*I*cos(theta)+2*R*sin(theta)
     f1<-1/(pi*sigma^2)*(sigma^2/2+c/2*exp((c^2)/(4*sigma^2))*(sqrt(pi)*sigma/2-sqrt(pi)*sigma/2*erf(-c/(2*sqrt(sigma)))))*exp(-(I^2+R^2)/(sigma^2))
     if (I==10){p1<-plot(theta, f1, type="n", main=expression(paste(sigma,'=15 ',',R=15')),cex.main=1,ylim=c(0,3.5),ylab="Prob Density",xlab='theta(rad)')}
     lines(theta, f1, type="l",col=rgb(0.5,1/90*I,1/90*I))
     col_leg<-c(col_leg,rgb(0.5,1/90*I,1/90*I))
}
legend(-3,3,legend=c('I=10','30','50','70','90'),
       lty=1, cex=0.7,col=col_leg)
R<-70
for (I in seq(10,100,20)){
  c<-2*I*cos(theta)+2*R*sin(theta)
  f1<-1/(pi*sigma^2)*(sigma^2/2+c/2*exp((c^2)/(4*sigma^2))*(sqrt(pi)*sigma/2-sqrt(pi)*sigma/2*erf(-c/(2*sqrt(sigma)))))*exp(-(I^2+R^2)/(sigma^2))
  if (I==10){p2<-plot(theta, f1, type="n", main=expression(paste(sigma,'=15 ',',R=70')),cex.main=1,ylim=c(0,4.5),ylab="Prob Density",xlab='theta(rad)')}
  lines(theta, f1, type="l",col=rgb(0.5,1/90*I,1/90*I))
  col_leg<-c(col_leg,rgb(0.5,1/90*I,1/90*I))
}
legend(-3,4,legend=c('I=10','30','50','70','90'),
       lty=1, cex=0.7,col=col_leg)
title('Effect of MRI Noise on CDI Phase Distribution',outer=TRUE,line=-1)