######### Call Change_GGS.R Function ############
source('~/R/rfunctions.R');
source('~/R/Change_GGS.R');

########### data ################
#data1=generate_data_1cp_smallsig();
data1=generate_data_1cp_smallsig();

########### prior ini values ########
ite=200;
mu_ini= c(32, 23, 15, 7);
tau_ini= c(2, 2, 1, 1);
pi_ini= c(0.25, 0.25, 0.25, 0.25);
MaxSeg=10; #30
minLen=100;  #10
MaxSeg2=50;
minLen2=20
avg_range=c(100,199) #range of iterations to average

#return values result$CK, result$GK, result$CH_POS, result$Segs

result_a1=Change_GGS(data_nocp, ite, mu_ini, tau_ini, pi_ini, MaxSeg, minLen);
average_a1=average_profile(avg_range, result_a1$CK, result_a1$CH_POS, result_a1$Segs);
result_a2=Change_GGS(data_nocp, ite, mu_ini, tau_ini, pi_ini, MaxSeg2, minLen2);
average_a2=average_profile(avg_range, result_a2$CK, result_a2$CH_POS, result_a2$Segs);

par(mfrow=c(1,2)); 
plot(1:length(data_nocp), data_nocp, xlab="Data points", ylab="Mean", main="No change-point");
lines(1:length(average_a1), average_a1, type="l", col="red");
lines(1:length(average_a2), average_a2, type="l", col="blue");
legend("bottomleft", legend=c("Max no=10, Min len=100", "Max no=50, Min len=20"),
       col=c("red", "blue"), lty=1:1, cex=1.2)

result_b1=Change_GGS(data_1cp, ite, mu_ini, tau_ini, pi_ini, MaxSeg, minLen);
average_b1=average_profile(avg_range, result_b1$CK, result_b1$CH_POS, result_b1$Segs);
result_b2=Change_GGS(data_1cp, ite, mu_ini, tau_ini, pi_ini, MaxSeg2, minLen2);
average_b2=average_profile(avg_range, result_b2$CK, result_b2$CH_POS, result_b2$Segs);

#par(mfrow=c(1,1)); 
plot(1:length(data_1cp), data_1cp, xlab="Data points", ylab="Mean", main="One change-point");
lines(1:length(average_b1), average_b1, type="l", col="red");
lines(1:length(average_b2), average_b2, type="l", col="blue");
legend("bottomleft", legend=c("Max no=10, Min len=100", "Max no=50, Min len=20"),
       col=c("red", "blue"), lty=1:1, cex=1.2)

data_4cp=c(rnorm(400,mean=10,sd=1),rnorm(50,mean=20,sd=1),rnorm(350,mean=5,sd=1),rnorm(50,mean=30,sd=1),rnorm(150,mean=20,sd=1));

########################################
mu_ini1= c(28, 18, 8, 4);
result2=Change_GGS(data1, ite, mu_ini1, tau_ini, pi_ini, MaxSeg, minLen);
average2=average_profile(avg_range, result2$CK, result2$CH_POS, result2$Segs);

############ plot graphs ############
par(mfrow=c(1,1)); 
plot(1:length(data1), data1, xlab="Data points", ylab="Mean", main="One change point at location 18");
#lines(1:length(data1), data1);
lines(1:length(average1), average1, type="l", col="red");
lines(1:length(average2), average2, type="l", col="blue");
legend("bottomleft", legend=c("Initial levels (15,7)", "Initial levels (8,4)"),
       col=c("red", "blue"), lty=1:1, cex=1.2)


########### 3 change points ##############
data3=generate_data_3cp();
mu_ini1= c(33, 21, 12, 4);
ite=500
result3=Change_GGS(data3, ite, mu_ini1, tau_ini, pi_ini, MaxSeg, minLen);

data3_GK=c(rep(2,250), rep(3, 250), rep(4, 250), rep(1, 250));

per_correct3=numeric();

for (i in 1:500){
  it_temp=numeric();
  it_G=numeric();
  for (m in 1:result3$Segs[i]){
    if (result3$Segs[i]==1){
      it_temp=rep(result3$GK[i, 1], (result3$CH_POS[i,2]-result3$CH_POS[i, 1]));
    }
    else{
      it_temp= c(it_temp, rep(result3$GK[i, m], (result3$CH_POS[i,m+1]-result3$CH_POS[i, m])));
    }
  }
  it_G=c(it_temp, result3$GK[i,result3$Segs[i]]);
  #per_correct[i]=length(which(it_G==data_GK))/1000*100
  per_correct3[i]=sum(it_G == data3_GK)/1000*100
}  

avg_range=c(100,200) #range of iterations to average
average3=average_profile(avg_range, result3$CK, result3$CH_POS, result3$Segs);

par(mfrow=c(1,2)); 
plot(1:length(data3), data3, xlab="Data points", ylab="Mean", main="3 change point at even intervals");
#lines(1:length(data1), data1);
lines(1:length(average3), average3, type="l", col="red");
plot(1:500, per_correct3, xlab="Iterations", ylab="% of correct group assignments", main="Group assignments vs iteration")

par(mfrow=c(2,2)); 
plot(1:500, result3$MUg[1:500,1], xlab="Iteration", ylab="Mean level", main="Group 1 mean vs interation (init=33)")
plot(1:500, result3$MUg[1:500,2], xlab="Iteration", ylab="Mean level", main="Group 2 mean vs interation (init=21)")
plot(1:500, result3$MUg[1:500,3], xlab="Iteration", ylab="Mean level", main="Group 3 mean vs interation (init=12)")
plot(1:500, result3$MUg[1:500,4], xlab="Iteration", ylab="Mean level", main="Group 4 mean vs interation (init=4)")


#plot(1:500, result3$MUg[1:500,1])
#plot(1:500, result3$CK[1:500,2])
#plot(1:500, result3$CK[1:500,3])
#plot(1:500, result3$CK[1:500,4])


#########################################################
#################### NOT USES ###########################
result2=Change_GGS(data1, ite, mu_ini, tau_ini, pi_ini, MaxSeg, minLen);

tau_ini= c(1, 1, 1, 1);
pi_ini= c(0.25, 0.25, 0.25, 0.25);
ite=50;

result1=Change_GGS(data1, ite, mu_ini, tau_ini, pi_ini, MaxSeg, minLen);


legend("bottomright", legend=c("Min Length=8", "Min Length=5"),
       col=c("blue", "green"), lty=1:1, cex=0.8)

lines(1:length(profile_data), profile_data, type="l", col="blue");
lines(1:length(profile_data), profile_data, type="l", col="green");
lines(1:length(profile_data), profile_data, type="l", col="yellow");
lines(1:length(profile_data), profile_data, type="l", col="orange")

########################################################


# number range to average
lower=200;
matrix_index=199; # ones less than lower
#length = (lower-matrix_index) : (higher-matrix_index)
higher=400;
best50=50;

p_data_m2 <- matrix(0, nrow=higher-lower+1, ncol=length(data1));
best_100_2 <- matrix(0, nrow=best50, ncol=length(data1));

mserror2<- matrix(0, nrow=higher-lower+1, ncol=2);
mse_best2<-matrix(0, nrow=higher-lower+1, ncol=2);

for (i in lower:higher) {
  
  p_tmp2=numeric();
  for (m in 1:result2$Segs[i]){
    p_tmp2= c(p_tmp2, rep(result2$CK[i, m], (result2$CH_POS[i,m+1]-result2$CH_POS[i, m])));
  }
  p_tmp2=c(p_tmp2, result2$CK[i,result2$Segs[m]]);
  p_data_m2[i-matrix_index, 1:length(p_tmp2)]=p_tmp2;
  
  mserror2[i-matrix_index,1]=sum((p_tmp2-data1)^2);
  mserror2[i-matrix_index,2]=i;
}

mse_best2=mserror2[order(mserror2[,1]),];
profile_data2= colMeans(p_data_m2);
#choose best 50
for (k in 1:50){
  best_100_2[k,1:length(data1)]= p_data_m2[mse_best2[k,2]-lower+1,1:length(data1)];
}

best_profile2=colMeans(best_100_2);

par(mfrow=c(1,1)); 
plot(1:length(data1), data1, xlab="Data points", ylab="Mean", main="Raw nanopore data between location 200-600");
#lines(1:length(data1), data1);
lines(1:length(profile_data2), profile_data2, type="l", col="yellow");
lines(1:length(best_profile2), best_profile2, type="l", col="green")


