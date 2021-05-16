#########################################################
###################  Generate Data ######################
#########################################################
generate_data<-function() {
  Levels4=c(10, 15, 25, 20);
  datavar=c(2, 2, 2, 2);
  probs4=c(0.4, 0.1, 0.2, 0.3);

  no_segment=5;
  seg_length=100;
  data=numeric();
  
  for (i in 1:no_segment) {
    no=sample(1:4, 1, FALSE, probs4);
    dmean=Levels4[no];
    dvar=datavar[no];
    oneseg=rnorm(seg_length, mean=dmean, sd=sqrt(dvar));
    data=c(data,oneseg);
  }
  plot(1:(no_segment*seg_length), data);
  #hist(data)
  return(data)  
}

## generate data with 1 change point
generate_data_1cp_smallsig<-function() {
  
  data=numeric();
  
  CK=c(30, 20, 10, 5);
  sig2=c(2, 2, 1, 1 ); #fixed
  #SIG=c(8, 6, 5, 3);
  dataLen=1000;   #fixed
  
  C1=sample(1:4, 1);
  C2=sample(1:4, 1);
  s=sample(1:dataLen, 1)

  data=c(rnorm(s-1, mean=CK[C1], sd=sig2[C1]), rnorm(dataLen-s+1, mean=CK[C2], sd=sig2[C2]));

  print(paste("Change point s is at ", s));
  plot(1:dataLen, data);
  
  return(data);
  #return(list("data"=data, "CPLoc"=s));
}

## generate data with 1 change point
generate_data_1cp_largesig<-function() {
  
  data=numeric();
  
  CK=c(30, 20, 10, 5);
  #sig2=c(5, 3, 2, 1 ); #fixed
  SIG=c(8, 6, 5, 3);
  dataLen=1000;   #fixed
  
  C1=sample(1:4, 1);
  C2=sample(1:4, 1);
  s=sample(1:dataLen, 1)
  
  data=c(rnorm(s-1, mean=CK[C1], sd=SIG[C1]), rnorm(dataLen-s+1, mean=CK[C2], sd=SIG[C2]));
  
  print(paste("Change point s is at ", s));
  plot(1:dataLen, data);
  
  return(data);
  #return(list("data"=data, "CPLoc"=s));
}

## generate data with 2 change points
generate_data_2cp<-function() {
  
  data=numeric();
  
  CK=c(30, 20, 10, 5);
  sig=c(1, 1, 1, 1);
  seg_hei=sample(CK, 3);
  seg_sig=sample(sig, 3); #fixed
  dataLen=1000;   #fixed
  
  data=c(rnorm(300, mean=seg_hei[1], sd=sig), rnorm(550, mean=seg_hei[2], sd=seg_sig[2]), rnorm(150, mean=seg_hei[3], sd=seg_sig[3]));
  plot(1:dataLen, data);
  
  return(data);
  #return(list("data"=data, "CPLoc"=s));
}

## generate data with 3 change points
generate_data_3cp<-function() {
  
  data=numeric();
  
  CK=c(30, 20, 10, 5);
  seg_hei=sample(CK, 4);
  sig2=1; #fixed
  dataLen=1000;   #fixed
  
  data=c(rnorm(250, mean=seg_hei[1], sd=sig2), rnorm(250, mean=seg_hei[2], sd=sig2), rnorm(250, mean=seg_hei[3], sd=sig2), rnorm(250, mean=seg_hei[4], sd=sig2));
  plot(1:dataLen, data);
  
  return(data);
  #return(list("data"=data, "CPLoc"=s));
}

## generate data iwth 3 to 6 change points
generate_data_3to6cp<-function() {
  
  data=numeric();
  
  hei=c(30, 20, 10, 5);
  sig2=1; #fixed
  dataLen=1000;   #fixed
  no_cp=sample(3:6, 1);
  hei_index=sample(1:4, no_cp+1, TRUE);
  seg_hei=hei[hei_index];
  cp_pos=c(1, sort(sample(2:(dataLen-1), no_cp)), dataLen);
  
  for (i in 1:(no_cp+1)){
    if (i== (no_cp)){
      seg_len=cp_pos[i+1]-cp_pos[i]+1;
    }
    else{
     seg_len=cp_pos[i+1]-cp_pos[i];
    }
    seg=rnorm(seg_len, mean=seg_hei[i], sd=sig2);
    data=c(data, seg);
  }
  #print(paste("Data length:", length(data)));
  print(paste("Change points at: ", cp_pos));
  plot(1:dataLen, data);
  
  return(data);
  #return(list("data"=data, "CPLoc"=s));
}

#########################################################
###################  Other Functions ####################
#########################################################
choose_a_group <-function(probs, num_probs){
  ran<-runif(1);
  sum=0.0;
  
  #print(paste("probs=",probs));
  #print(paste("num_probs=",num_probs));
  #print(paste("ran=",ran));
  

  for(i in 1: num_probs){

    sum = sum+probs[i];
    #print(paste("sum=",sum));

    if (sum >= ran) { return (i);}
  }
  if(sum(probs) != 1){
    print("ERROR - probability sum not = ONE");
  }
  return (i);
  
}


box_muller_rnorm <-function(mean, stddev){
  
  while(TRUE){
    x<-2*runif(1)-1;
    #print(paste("x=",x));
    y<-2*runif(1)-1;
    #print(paste("y=",y));
    r=x*x+y*y;
    #print(paste("r=",r));
    if ((r >=1.0) | (r==0.0)){
      #print("ERROR in box_muller_rnorm function");      
    }
    else {
      s=sqrt(-2*log(r)/r);
      break;
    }
  }
  return (mean+x*s*stddev);
}

marsaglia_rgamma<-function(alpha0, beta0){
  # Marsaglia's (1961) method: gam(a) = gam(a+1)*U^(1/a)
  
  if (alpha0 >=1) {
    d = alpha0-1/3
    c = 1/sqrt(9*d);
    while (TRUE){
      x = rnorm(1);
      v = 1 + c*x;
      while(TRUE){
        if (v <= 0){
          x = rnorm(1);
          v = 1.0 + c*x;
        }
        else{break;}
      }
      v = v*v*v;
      u = runif(1);
      #print(paste("u=",u));
      #print(paste("v=",v));
      #print(paste("x=",x));
      #print(paste("d=",d));
      if (u < 1.0-0.0331*(x*x)*(x*x)) return (d*v/beta0);
      if (log(u) < 0.5*x*x+d*(1.0-v+log(v))) return (d*v/beta0);
    }
  }
  else {
    print("else");
    x=marsaglia_rgamma(alpha0+1,beta0);
    x=x*runif(1)^(1/alpha0);
    return (x)
  }

}
  
ran_beta<-function(alpha0, beta0){
  prob1=marsaglia_rgamma(alpha0, 1);
  prob2=marsaglia_rgamma(beta0, 1);
  return (prob1/(prob1+prob2))
  
}

ran_dir<-function(num, num_group){
  total=0;
  prob <- vector();
  group_probs <- vector();
  
  for (i in 1:num_group) {
    prob[i]=marsaglia_rgamma(num[i]+1,1);
    #print(paste("prob=",prob[i]));
    total=total+ prob[i];
    #print(paste("total=",total));
  }
  for (j in 1:num_group) {
    group_probs[j]=prob[j]/total;
  }
  return (group_probs)
}


##################################################################################
# plot average change point profile

average_profile<-function(range, mean_levels, ch_positions, Seg_nos){
  # number range to average
  lower=range[1];
  matrix_index=lower-1; # ones less than lower
  higher=range[2];
  #length = (lower-matrix_index) : (higher-matrix_index)
  
  p_data_m <- matrix(0, nrow=higher-lower+1, ncol=length(data1));

  for (i in lower:higher) {
    
    p_tmp=numeric();
    for (m in 1:Seg_nos[i]){
      p_tmp= c(p_tmp, rep(mean_levels[i, m], (ch_positions[i,m+1]-ch_positions[i, m])));
    }
    p_tmp=c(p_tmp, mean_levels[i, Seg_nos[m]]);
    p_data_m[i-matrix_index, 1:length(p_tmp)]=p_tmp;
  }

  profile_data= colMeans(p_data_m);
  return(profile_data)
}

