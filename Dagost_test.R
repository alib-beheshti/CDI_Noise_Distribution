library(VGAM)#erf function
library(moments)#for d agostino test
library(spatstat)#for im
rm(list=ls())
graphics.off()
theta<-seq(-pi,pi,0.0001);
kurt_cell<-NULL;
quant_sample_vect_index<-array(0,dim=c(length(seq(1,100,5)),length(seq(1,100,5)),length(seq(-pi,pi,0.002))))
vect_sample_freq<-array(0,dim=c(length(seq(1,100,5)),length(seq(1,100,5)),length(seq(-pi,pi,0.002))))
#quantize the distribution values :
sigma<-10
j<-0
for (R in seq(1,100,5)){
    j<-j+1
    k<-0
    for (I in seq(1,100,5)){
        k<-k+1
        l<-0
        for (theta_sample in seq(-pi,pi,0.002)){
            l<-l+1
            c<-2*I*cos(theta_sample)+2*R*sin(theta_sample)
            f<-1/(pi*sigma^2)*(sigma^2/2+c/2*exp((c^2)/(4*sigma^2))*(sqrt(pi)*sigma/2-sqrt(pi)*sigma/2*erf(-c/(2*sqrt(sigma)))))*exp(-(I^2+R^2)/(sigma^2))
            vect_sample_freq[j,k,l]<-f
        }
        vect_sample_freq[j,k,]<-round(vect_sample_freq[j,k,],1)
        quant_sample_vect_index[j,k,]=vect_sample_freq[j,k,]*10
    }
}
##Make the samples:
kurt_cube=array(0,dim=c(length(seq(1,100,5)),length(seq(1,100,5))))
skew_cube=array(0,dim=c(length(seq(1,100,5)),length(seq(1,100,5))))
test_cube=array(0,dim=c(length(seq(1,100,5)),length(seq(1,100,5))))
sigma<-10
j<-0
for (R in seq(1,100,5)){
    j<-j+1
    k<-0
    for (I in seq(1,100,5)){
        k<-k+1
        l<-0
        sample_vec_atpoint=NULL
        for (theta_sample in seq(-pi,pi,0.002)){
            l<-l+1
            val<-quant_sample_vect_index[j,k,l]
            if (val>=1){
               for (kk in seq(1,val)){
                   sample_vec_atpoint=cbind(sample_vec_atpoint,theta_sample)
               }
            }  
        }    
        kurt_cube[j,k]<-kurtosis(sample_vec_atpoint[1,])
        skew_cube[j,k]<-skewness(sample_vec_atpoint[1,])
        test_cube[j,k]<- agostino.test(sample_vec_atpoint[1,],alternative="two.sided")$p.value
    }    
}
a=array(0,dim=c(length(seq(1,100,5)),length(seq(1,100,5))))
a=kurt_cube
x=seq(1,100,5)
y=seq(1,100,5)
par(mfrow=c(1,3))
image(im(a[nrow(a):1,]),main=NULL)
title(main='Kurtosis',line=-2,xlab='Imaginary')
mtext("Real", side=2, line=3,cex=0.7)
ticks = c(1,5,10,15,20)
axis(1, at=ticks,labels=c('1','25','50','75','100'),line=-5.5)
axis(2, at=ticks,labels=c('100','75','50','25','1'))
a=skew_cube
x=seq(1,100,5)
y=seq(1,100,5)
image(im(a[nrow(a):1,]),main=NULL)
title(main='Skewness',line=-2,xlab='Imaginary')
mtext("Real", side=2, line=3,cex=0.7)
axis(1, at=ticks,labels=c('1','25','50','75','100'),line=-5.5)
axis(2, at=ticks,labels=c('100','75','50','25','1'))
a=test_cube
x=seq(1,100,5)
y=seq(1,100,5)
a[a>0.05]=1
a[a<=0.05]=0
image(im(a[nrow(a):1,]),main=NULL)
title(main='Agostino Test\nGaussian:1,Non Gaussian:0',line=-2,xlab='Imaginary')
mtext("Real", side=2, line=3,cex=0.7)
axis(1, at=ticks,labels=c('1','25','50','75','100'),line=-5.5)
axis(2, at=ticks,labels=c('100','75','50','25','1'))
mtext(expression(bold(paste('Agostino Test Based Sensitivity Analysis for ',sigma,'=10'))),side=3,line=-2,outer=TRUE)
