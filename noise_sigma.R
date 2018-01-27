library(VGAM)
rm(list=ls())
graphics.off()
theta<-seq(-pi,pi,0.0001)
col_leg<-NULL
par(mfrow=c(1,2))
R<-70
I<-50
for (sigma in seq(4,24,4)){
     c<-2*I*cos(theta)+2*R*sin(theta)
     f1<-1/(pi*sigma^2)*(sigma^2/2+c/2*exp((c^2)/(4*sigma^2))*(sqrt(pi)*sigma/2-sqrt(pi)*sigma/2*erf(-c/(2*sqrt(sigma)))))*exp(-(I^2+R^2)/(sigma^2))
     if (sigma==4){p1<-plot(theta, f1, type="n", main="Pixel Phase Prob Distribution",cex.main=0.8,ylim=c(0,14),ylab="Prob Density",xlab='theta(rad)',yaxt="n")
     axis(2,cex.axis=0.7)
     }
     lines(theta, f1, type="l",col=rgb(1/40*sigma,1/24*sigma,1/24*sigma))
     col_leg<-c(col_leg,rgb(1/40*sigma,1/24*sigma,1/24*sigma))
}
legend(-3, 12, legend=c(expression(paste(sigma,'=4 ')),'8','12','16','20','24'),
       lty=1, cex=0.7,col=col_leg)
for (sigma in seq(4,24,4)){
  c<-2*I*cos(theta)+2*R*sin(theta)
  f1<-1/(pi*sigma^2)*(sigma^2/2+c/2*exp((c^2)/(4*sigma^2))*(sqrt(pi)*sigma/2-sqrt(pi)*sigma/2*erf(-c/(2*sqrt(sigma)))))*exp(-(I^2+R^2)/(sigma^2))
  if (sigma==4){p2<-plot(theta, f1, type="n", main="Zoomed Pixel Phase Prob Distribution",cex.main=0.8,ylim=c(0,14),xlim=c(0,2),ylab="Prob Density",xlab='theta(rad)',yaxt="n")
  axis(2,cex.axis=0.7)
  }
    lines(theta, f1, type="l",col=rgb(1/40*sigma,1/24*sigma,1/24*sigma))
  col_leg<-c(col_leg,rgb(1/40*sigma,1/24*sigma,1/24*sigma))
}
legend(0.1, 12, legend=c(expression(paste(sigma,'=4 ')),'8','12','16','20','24'),
       lty=1, cex=0.7,col=col_leg)
title('Effect of MRI Noise on CDI Phase Distribution',outer=TRUE,line=-1)