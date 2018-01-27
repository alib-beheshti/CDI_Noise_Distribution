library(VGAM)
rm(list=ls())
graphics.off()
theta<-seq(-pi,pi,0.0001)
i<-0
par(mfrow=c(1,2))
for (mag in c(1600,2500)){
  for (I in seq(10,30,10)){
     sigma<-8
     R<-sqrt(mag-I^2)
     c<-2*I*cos(theta)+2*R*sin(theta)
     f1<-1/(pi*sigma^2)*(sigma^2/2+c/2*exp((c^2)/(4*sigma^2))*(sqrt(pi)*sigma/2-sqrt(pi)*sigma/2*erf(-c/(2*sqrt(sigma)))))*exp(-(I^2+R^2)/(sigma^2))
     if (mag==1600 && I==10){p1<-plot(theta, f1, type="n", main="Pixel Prob Distribution",cex.main=0.8,ylim=c(0,4),ylab="Prob Density")}
     if (mag==1600){lines(theta, f1, type="l",col="red")}
     if (mag==2500){lines(theta, f1, type="l",col="blue")}
     i<-i+1
  }   
}
legend(-3, 4, legend=c("Mag=40", "Mag=50"),col=c("red", "blue"),
       lty=1, cex=0.7)
for (mag in c(1600,2500)){
  for (I in seq(10,30,10)){
    sigma<-8
    R<-sqrt(mag-I^2)
    c<-2*I*cos(theta)+2*R*sin(theta)
    f1<-1/(pi*sigma^2)*(sigma^2/2+c/2*exp((c^2)/(4*sigma^2))*(sqrt(pi)*sigma/2-sqrt(pi)*sigma/2*erf(-c/(2*sqrt(sigma)))))*exp(-(I^2+R^2)/(sigma^2))
    if (mag==1600 && I==10){p1<-plot(theta, f1, type="n", main="Zoomed Pixel Prob Distribution",cex.main=0.8,ylim=c(0,4),xlim=c(0,2),ylab="Prob Density")}
    if (mag==1600){lines(theta, f1, type="l",col="red")}
    if (mag==2500){lines(theta, f1, type="l",col="blue")}
    i<-i+1
  }   
}
legend(-3, 4, legend=c("Mag=40", "Mag=50"),col=c("red", "blue"),
       lty=1, cex=0.7)
