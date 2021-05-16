######  MRes graphs for thesis real data ########

source('~/R/rfunctions.R');
source('~/R/Change_GGS.R');

## real data ##
df_data<-read.delim(header = TRUE, "C:/MRes/Data/fast5s/test.txt");  
my_data<-df_data[,1];
data_real0<-my_data[200:600];

#real data graphs
par(mfrow=c(1,2)); 
plot(1:5000, my_data[1:5000], t='l', ylab="Raw signal", xlab="Data points", main="(a) Sample nanopore signal");
plot(200:600, my_data[200:600], t='l', ylab="Raw signal", xlab="Data points", main="(b) A section of the same data");
par(mfrow=c(1,1)); 
plot(1:401, data_real0, t='l');

ite=200;
ite01=600;
mu_ini_real= c(700, 600, 500, 400);
tau_ini_real= c(15, 15, 15, 15);
pi_ini= c(0.25, 0.25, 0.25, 0.25);
MaxSeg_real=100; #30
minLen_real=10;  #10
result0=Change_GGS(data_real0, ite, mu_ini_real, tau_ini_real, pi_ini, MaxSeg_real, minLen_real);
#result01=Change_GGS(data_real0, ite01, mu_ini_real, tau_ini_real, pi_ini, MaxSeg_real, minLen_real);
avg_range0=c(100,199) #range of iterations to average
average0=average_profile(avg_range0, result0$CK, result0$CH_POS, result0$Segs);
par(mfrow=c(1,1)); 
plot(1:length(data_real0), data_real0, xlab="Data points", ylab="Mean", main="A section of raw nanopre data");
lines(1:length(average0), average0, type="l", col="red")


### plot 3 ranges of iterations ###
avg_range01=c(100,250) #range of iterations to average
average_real1=average_profile(avg_range01, result01$CK, result01$CH_POS, result01$Segs);

avg_range02=c(250,450) #range of iterations to average
average_real2=average_profile(avg_range02, result01$CK, result01$CH_POS, result01$Segs);

avg_range03=c(400,550) #range of iterations to average
average_real3=average_profile(avg_range03, result01$CK, result01$CH_POS, result01$Segs);

par(mfrow=c(1,1)); 
plot(1:length(data_real0), data_real0, xlab="Data points", ylab="Mean", main="Comparison of methods - Bayesian sampler vs changepoint.np");
lines(1:length(average_real1), average_real1, type="l", col="red");
lines(1:length(average_real2), average_real2, type="l", col="blue");
lines(1:length(average_real3), average_real3, type="l", col="green");

########## plot result$Segs ###############
par(mfrow=c(1,2)); 
plot(1:200, result0$Segs, xlab="No of iteration", ylab="No Segments", main="Real nanopore data");
lines(1:200, result0$Segs, type="l", col="red");
plot(1:200, result1$Segs, xlab="No of iteration", ylab="No Segments", main="Simulated data - one change");
lines(1:200, result1$Segs, type="l", col="blue")
#legend("bottomright", legend=c("Real nanopore data", "Simulated data - one change point"),
#       col=c("red", "blue"), lty=1:1, cex=0.8)


################# compare with other packages #################

# changepoint
library(changepoint)

fit_changepoint = cpt.meanvar(data_real0, penalty="Manual",pen.value="2*log(n)",method="PELT",Q=10,class=FALSE);
plot(fit_changepoint)

#changepoin.np
library(changepoint.np)
cps_changepoint=changepoint.np::cpt.np(data_real0)
abline(v = cps_changepoint@cpts, col="blue")
legend("bottomleft", legend=c("Bayesian sampler", "changepoint.np"),
       col=c("red", "blue"), lty=1:1, cex=0.7)

#bcp
library(bcp)
fit_bcp = bcp(data_real0, d = 1000)
par(mfrow=c(2,1)); 
plot(1:length(data_real0), data_real0, xlab="Data points", ylab="Mean", main="Comparison of methods - Bayesian sampler vs package");
lines(1:length(average_real1), average_real1, type="l", col="red");
plot(fit_bcp)

library(cpm)
#fit_cpm = detectChangePoint(df$y, cpmType = "Student")  # a single change point
#fit_cpm = processStream(df$y, cpmType = "Mann-Whitney")  # Detects three
fit_cpm = processStream(data_real0, cpmType = "Student")  # Multiple change points
fit_cpm$changePoints
par(mfrow=c(1,1)); 
plot(1:length(data_real0), data_real0, xlab="Data points", ylab="Mean", main="Comparison of methods - Bayesian sampler vs cpm");
lines(1:length(average_real1), average_real1, type="l", col="red");

abline(v = fit_cpm$changePoints, col="blue")
legend("bottomleft", legend=c("Bayesian sampler", "cpm package"),
       col=c("red", "blue"), lty=1:1, cex=0.7)

#strucchange
library(strucchange)
fit_bp = breakpoints(y ~ 1, data = data_real0, breaks = 20)
