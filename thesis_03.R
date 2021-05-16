############ small segment #################
##### change point locaion 18 ########

#data_loc18

########### prior ini values ########
ite=200;
mu_ini= c(32, 23, 15, 7);
tau_ini= c(2, 2, 1, 1);
pi_ini= c(0.25, 0.25, 0.25, 0.25);
MaxSeg=20; #30
minLen=50;  #10
avg_range=c(100,199) #range of iterations to average

#return values result$CK, result$GK, result$CH_POS, result$Segs

result1=Change_GGS(data_loc18, ite, mu_ini, tau_ini, pi_ini, MaxSeg, minLen);
average1=average_profile(avg_range, result1$CK, result1$CH_POS, result1$Segs);

par(mfrow=c(1,1)); 
plot(1:length(data1), data1, xlab="Data points", ylab="Mean", main="One change point at location 18");
#lines(1:length(data1), data1);
lines(1:length(average1), average1, type="l", col="red")
