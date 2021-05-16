##################################################################
############ Change Point Generalised Gibbs sampler ##############
##################################################################

Change_GGS<-function(data2, max_ite, mu_g_ite, tau_g_ite, gpi_ite, MaxSeg, minLen) {

  #library(MCMCpack);
  #library(Rmpfr);
  source('~/R/rfunctions.R');
  
  #############  initialisation ##############
  #max_ite=200;
  N_types=4;  # 4 types of segment
  #MaxSeg=10;
  #minLen=15;
  
  #MaxSeg=20; #30
  #minLen=50;  #10
  offset=2;
  ex=2;  # extra space to store variable values 
  minLenCounter=0;
  
  CK <- matrix(0, nrow=max_ite+ex, ncol=MaxSeg+ex);
  GK <- matrix(0, nrow=max_ite+ex, ncol=MaxSeg+ex);
  CH_POS <- matrix(0, nrow=max_ite+ex, ncol=MaxSeg+ex);
  MUg <- matrix(0, nrow=max_ite+ex, ncol=4);
  TAUg2 <- matrix(0, nrow=max_ite+ex, ncol=4);
  gPI <- matrix(0, nrow=max_ite+ex, ncol=4);
  del_accept_prob <- matrix(0, nrow=max_ite+ex, ncol=MaxSeg+ex);
  del_Accept <- matrix(0, nrow=max_ite+ex, ncol=MaxSeg+ex);
  insert_accept_prob <- matrix(0, nrow=max_ite+ex, ncol=MaxSeg+ex);
  ins_Accept <- matrix(0, nrow=max_ite+ex, ncol=MaxSeg+ex);
  
  #local variables - values not kept
  cpos_ite=numeric();
  ck_ite=numeric();
  gk_ite=numeric();
  #mu_g_ite=numeric();
  #tau_g_ite=numeric();
  #gpi_ite=numeric();
  phi=numeric();
  sig=numeric();
  Segs=numeric();
  
  ### prior values ###
  a0=1;
  b0=1;
  alpha0=3;
  beta0=3;
  u0=3;
  v0=3;
  phi[1]=ran_beta(a0, b0);
  sig[1]=1/marsaglia_rgamma(u0, v0);
  #tau_g_ite=1/rgamma(4, alpha0, scale=beta0);  #c(1, 1, 1, 1)
  #tau_g_ite=c(10, 10, 10, 10); 
  #gpi_ite=c(0.1, 0.5, 0.2, 0.1); #gPI=rdirichlet(1, c(1, 1, 1, 1));
  #mu_g_ite=c(900, 700, 500, 400);  # real values = (30, 20, 10, 5);
  #values from input values
  MUg[1, 1:4]<-matrix(mu_g_ite, nrow=1, ncol=4, byrow = TRUE); 
  TAUg2[1, 1:4]<-matrix(tau_g_ite, nrow=1, ncol=4, byrow = TRUE);
  gPI[1, 1:4]<-matrix(gpi_ite, nrow=1, ncol=4, byrow = TRUE);
  
  ###### starting value - number of segment #######
  NoSeg=1;  # assume no chang points to start with
  #NoSeg=3; 
  
  ## randomly select change point(s) initally 
  if (NoSeg==1) {
    cpos_ite=c(1, length(data2));
    gk_ite[1]=choose_a_group(gpi_ite,N_types); 
    ck_ite[1]=box_muller_rnorm(mu_g_ite[gk_ite[1]], sqrt(tau_g_ite[gk_ite[1]]));
  }
  if (NoSeg > 1) {
    cp_temp=sort(sample(2:length(data2), (NoSeg-1)), FALSE); #pick change point at random
    cpos_ite=c(1, cp_temp, length(data2));
    
    for (i in 1:NoSeg){
      #randomly select a group(Gk) and mean (Ck) for each segment #sample(1:4, 1, TRUE, gpi_ite);
      gk_ite[i]=choose_a_group(gpi_ite,N_types);   
      ck_ite[i]=box_muller_rnorm(mu_g_ite[gk_ite[i]], sqrt(tau_g_ite[gk_ite[i]])); 
    }
  }

  
  CH_POS[1, 1:length(cpos_ite)]<-matrix(cpos_ite, nrow=1, ncol=length(cpos_ite), byrow = TRUE);
  CK[1, 1:length(ck_ite)]<-matrix(ck_ite, nrow=1, ncol=length(ck_ite), byrow = TRUE); 
  GK[1, 1:length(gk_ite)]<-matrix(gk_ite, nrow=1, ncol=length(gk_ite), byrow = TRUE); 
  Segs[1]=NoSeg;
  
  #########################
  #### start segmention ###
  #########################
  #iteration counter l
  l=1;
  
  
  ## do max_ite number of iteration
  while (l < max_ite) {
    
    k=1;
    #print(paste("-------------iteration-------",l,  "- NoSeg - ", NoSeg));
    Segs[l]=NoSeg;
    
    #### start segmention ###
    
    while (k <= NoSeg){ 
      
      ########## Deletion Step ############## 
      
      # merged group
      if (k != 1) {
        G0=choose_a_group(gpi_ite,N_types);  
        C0=box_muller_rnorm(mu_g_ite[G0], sqrt(tau_g_ite[G0]));
        
        dlerr=numeric();
        se_d1=0;
        for (j in (cpos_ite[k-1]:(cpos_ite[k+1]-1))){
          #dlerr[j]=(-0.5)*log(2*pi*sig[l])-0.5/sig[l]*((data2[j]- C0)^2);
          dlerr[j]=-1/(2*sig[l])*((data2[j]- C0)^2);
          se_d1=se_d1+ dlerr[j];
        }
        tk=4*NoSeg-1+2*N_types+3-3; #### check if need to -3 as in Farhana's code?????
        tkm1=4*(NoSeg-1)-1+2*N_types+3;
        #delete_p1=log(1-phi[l])-log(tkm1)+log_d1-log(cpos_ite[k+1]-cpos_ite[k-1]);
        delete_p1=log(1-phi[l])/phi[l]-log(cpos_ite[k+1]-cpos_ite[k-1])-log(tkm1/tk)+se_d1;
        
        # keep 2 segments
        derr1=numeric();
        se_dleft=0;
        for (i in (cpos_ite[k-1]:(cpos_ite[k]-1))){
          derr1[i]=-1/(2*sig[l])*((data2[i]-ck_ite[k-1])^2);
          se_dleft=se_dleft+ derr1[i];
        }
        derr2=numeric();
        se_dright=0;
        for (i in cpos_ite[k]:(cpos_ite[k+1]-1)){
          derr2[i]=-1/(2*sig[l])*((data2[i]-ck_ite[k])^2);
          se_dright=se_dright+ derr2[i];
        }
        
        delete_p0=se_dleft+se_dright;
        
        #calculate delete probability
        
        max_prob=0.0;
        sum=0;
        choice=0;
        
        if (delete_p0>delete_p1){
          max_prob_del=delete_p0;
        }
        else {
          max_prob_del=delete_p1;
        }
        
        # fix precision of exp(x) when x is very small 
        #p0_del=mpfr((delete_p0-max_prob_del), precBits = 106);
        #p1_del=mpfr((delete_p1-max_prob_del), precBits = 106);
        
        p0_del_exp=exp(mpfr((delete_p0-max_prob_del), precBits = 106));
        p1_del_exp=exp(mpfr((delete_p1-max_prob_del), precBits = 106));
        sum_probs_del=p0_del_exp+p1_del_exp;
        probs0_del=p0_del_exp/sum_probs_del;
        probs1_del=p1_del_exp/sum_probs_del;
        
        #debug
        #print(paste("delete_p0 no change=", delete_p0));
        #print(paste("delete_p1 merge=", delete_p1));
        #print(paste("p0_del_exp =", p0_del_exp));
        #print(paste("p1_del_exp =", p1_del_exp));
        #print(paste("probs0_del =", probs0_del));
        #print(paste("probs1_del =", probs1_del));
        
        choice_del=choose_a_group(c(probs0_del, probs1_del), 2);
        
        #print(paste("choice_del =", choice_del));
        
        ############ DELETION acceptance #############
        cpos_ite_left=numeric();
        cpos_ite_right=numeric();
        gk_ite_left=numeric();
        gk_ite_right=numeric();
        ck_ite_left=numeric();
        ck_ite_right=numeric();
        
        # delete accepted
        if (choice_del==2){  
          
          #remove change point
          #print(paste("CH to be deleted=", cpos_ite[k]));
          remove_cp=rep(FALSE, NoSeg);
          remove_cp[k]=TRUE;
          cpos_ite=cpos_ite[!(remove_cp)];
          
          #update ck_ite & gk_ite
          
          if(NoSeg==2){
            gk_ite=G0;
            ck_ite=C0;
          }
          else {
            if(k==2){
              gk_ite_right=gk_ite[(k+1):NoSeg];
              gk_ite=c(G0, gk_ite_right);
              ck_ite_right=ck_ite[(k+1):NoSeg];
              ck_ite=c(C0, ck_ite_right);
            }
            else if(k==NoSeg){
              gk_ite_left=gk_ite[1:(NoSeg-2)];
              gk_ite=c(gk_ite_left, G0);
              ck_ite_left=ck_ite[1:(NoSeg-2)];
              ck_ite=c(ck_ite_left, C0);
            }
            else{
              gk_ite_left=gk_ite[1:(k-2)];
              gk_ite_right=gk_ite[(k+1):NoSeg];
              gk_ite=c(gk_ite_left, G0, gk_ite_right);
              ck_ite_left=ck_ite[1:(k-2)];
              ck_ite_right=ck_ite[(k+1):NoSeg];
              ck_ite=c(ck_ite_left, C0, ck_ite_right);
            }          
          } # else (NoSeg==2)
          
          #k=k-1;  
          NoSeg=NoSeg-1; 
          # k stays the same (ie go to next segment)
        }
        else if (choice_del==1){  
          #print(paste("Delete rejected - k=", k));
        }
        else {
          print(paste("Delete ERROR - k=", k));
        }
        
      } # deletion step - if (k!=1)
      
      ########### INSERATION ################
      
      if ((NoSeg < MaxSeg) && (k <= NoSeg)) { 
        
        if ((cpos_ite[k+1]-cpos_ite[k]) < minLen) {
          minLenCounter=minLenCounter+1;
          #k=k+1;
          #print(paste("Segment is less than minLen - length", cpos_ite[k+1]-cpos_ite[k]));
        }
        else {
          insert_flag=TRUE;
          max_loop=0;
          ## added while loop 22/02/2021
          while (insert_flag) {
            
            #break out of insert loop if no suitable point is found
            if (max_loop > 10)
            {
              break;
            }
            # no change point inserted
            lerr=numeric();
            se_p1=0;
            for (i in cpos_ite[k]:(cpos_ite[k+1]-1)){
              #lerr[i]=(-0.5)*log(2*pi*sig[l])-0.5/sig[l]*(data2[i]- ck_ite[k])^2;
              lerr[i]=-1/(2*sig[l])*(data2[i]- ck_ite[k])^2;
              se_p1=se_p1+lerr[i];
            }
            
            tk=4*NoSeg-1+2*N_types+3; #no of moves with no insertion
            tkp1=4*(NoSeg+1)-1+2*N_types+3; #no of moves with insertion(t(k+1))
            
            ### Equation 10 ### insertion P1 probability -inset_p1 is log scale
            #insert_p1=log((1-phi[l]))-log(tk)-log(cpos_ite[k+1]-cpos_ite[k]+1)+log_p1; 
            #######check length is +1 or -1 (Farhana's code) ???????????????
            insert_p1=log((1-phi[l])/phi[l])+log(tkp1/tk)-log(cpos_ite[k+1]-cpos_ite[k]-1)+se_p1;  
            
            # randomly pick a change point z between sk and dk-1
            z=sample((cpos_ite[k]+offset):(cpos_ite[k+1]-offset), 1, FALSE); #pick a change point
            G1=choose_a_group(gpi_ite,N_types);
            G2=choose_a_group(gpi_ite,N_types);
            #C1=rnorm(1, mean=gu_ite[G1], sd=sqrt(tau_g_ite[G1]));
            #C2=rnorm(1, mean=gu_ite[G2], sd=sqrt(tau_g_ite[G2]));
            C1=box_muller_rnorm(mu_g_ite[G1], sqrt(tau_g_ite[G1]));
            C2=box_muller_rnorm(mu_g_ite[G2], sqrt(tau_g_ite[G2]));
            
            lerr1=numeric();
            se_ileft=0;
            for (i in cpos_ite[k]:(z-1)){
              #lerr1[i]=(-0.5)*log(2*pi*sig[l])-0.5/sig[l]*(data2[i]-C1)^2;
              lerr1[i]=-1/(2*sig[l])*(data2[i]-C1)^2;
              se_ileft=se_ileft+ lerr1[i];
            }
            lerr2=numeric();
            se_iright=0;
            for (i in z:(cpos_ite[k+1]-1)){
              lerr2[i]=-1/(2*sig[l])*(data2[i]- C2)^2;
              se_iright=se_iright+ lerr2[i];
            }     
            
            insert_p0=+se_ileft+se_ileft;
            
            #calculate acceptane probability
            
            max_prob_ins=0.0;
            sum=0;
            choice_ins=0;
            
            if (insert_p0>insert_p1){
              max_prob_ins=insert_p0;
            }
            else {
              max_prob_ins=insert_p1;
            }
            
            p0_ins_exp=exp(mpfr((insert_p0-max_prob_ins), precBits = 106));
            p1_ins_exp=exp(mpfr((insert_p1-max_prob_ins), precBits = 106));
            sum_probs_ins=p0_ins_exp+p1_ins_exp;
            probs0_ins=p0_ins_exp/sum_probs_ins;
            probs1_ins=p1_ins_exp/sum_probs_ins;
            
            #print(paste("insert_p1 =", insert_p1));
            #print(paste("insert_p0 =", insert_p0));
            #print(paste("probs0_ins insert=", probs0_ins));
            #print(paste("probs1_ins NOT insert=", probs1_ins));
            
            choice_ins=choose_a_group(c(probs0_ins, probs1_ins), 2);
            
            cpos_ite_left=numeric();
            cpos_ite_right=numeric();
            gk_ite_left=numeric();
            gk_ite_right=numeric();
            ck_ite_left=numeric();
            ck_ite_right=numeric();
            
            # insert accepted
            if (choice_ins==1){
              insert_flag=FALSE;
              #print(paste("insert accepted - k=", k));
              #print(paste("CH to be inserted = ", z));
              #add change point
              cpos_ite_left=cpos_ite[1:k];
              cpos_ite_right=cpos_ite[(k+1):(NoSeg+1)];
              cpos_ite=c(cpos_ite_left, z, cpos_ite_right);
              
              #print(paste("cpos_ite_left =", cpos_ite_left));
              #print(paste("cpos_ite_right =", cpos_ite_right));
              #print(paste("cpos_ite =", cpos_ite));
              
              #update ck_ite & gk_ite
              if(NoSeg==1){
                gk_ite=c(G1, G2);
                ck_ite=c(C1, C2);
              }
              else{
                if (k==1){
                  gk_ite_right=gk_ite[(k+1):NoSeg];
                  gk_ite=c(G1, G2, gk_ite_right);
                  ck_ite_right=ck_ite[(k+1):NoSeg];
                  ck_ite=c(C1, C2, ck_ite_right);
                }
                else if (k==NoSeg){
                  gk_ite_left=gk_ite[1:(NoSeg-1)];
                  gk_ite=c(gk_ite_left, G1, G2);
                  ck_ite_left=ck_ite[1:(NoSeg-1)];
                  ck_ite=c(ck_ite_left, C1, C2);
                }
                else{
                  gk_ite_left=gk_ite[1:(k-1)];
                  gk_ite_right=gk_ite[(k+1):NoSeg];
                  gk_ite=c(gk_ite_left, G1, G2, gk_ite_right);
                  ck_ite_left=ck_ite[1:(k-1)];
                  ck_ite_right=ck_ite[(k+1):NoSeg];
                  ck_ite=c(ck_ite_left, C1, C2, ck_ite_right);
                }
              } #if(NoSeg==1) 
              
              NoSeg=NoSeg+1; 
              #k=k+1;
            } # if (insert accepted)
            else if (choice_ins==2){
              #k=k+1;
              #print(paste("Insert rejected - k=", k));
              max_loop=max_loop+1;
            }
            else
            { print(paste("Insert ERROR - k=", k)); }
            
          }  ### added end while loop 22/02/2021
        } # else (over min length) 
      }  # if (insert)
      
      k=k+1; 
    } # while (k <= NoSeg)
    
    l=l+1;  # increment counter before next iteration
    
    #############################################################
    ######## update parameters to use in next iteration #########
    
    #work out et2 (epsilon t square)
    et2=0;
    data_ck=numeric();
    for (m in 1:NoSeg){
      data_ck= c(data_ck, rep(ck_ite[m], (cpos_ite[m+1]-cpos_ite[m])));
    }
    data_ck=c(data_ck, ck_ite[NoSeg]);
    
    # sum of error (xt-ck) squared
    et2=sum((data2-data_ck)^2)
    #print(paste("ep =", ep));
    
    #update sig^2
    sig[l]=marsaglia_rgamma(u0+length(data2)/2, v0+0.5*et2);
    
    #update phi
    phi[l]=ran_beta(NoSeg,(length(data2)-NoSeg));   #rbeta(K, T-K) 
    prior_phi=phi[l]^(a0-1)*(1-phi[l])^(b0-1)/beta(a0, b0)
    
    # update pi - rdirichlet(b1+1, b2+11...)
    gpi_ite=ran_dir(c(length(which(gk_ite==1))+1, length(which(gk_ite==2))+1,length(which(gk_ite==3))+1, length(which(gk_ite==4))+1),N_types);
    
    
    ########## Update ck and gk ############## 
    
    # pr_ck=numeric();
    # max_prck=0;
    # #print(paste("ck_ite =", ck_ite));
    # 
    # for (j in 1:NoSeg){
    #   for (g in 1:N_types){
    #     se_temp=log(1/(sqrt(2*pi*tau_g_ite[g])))- (ck_ite[j]- mu_g_ite[g])^2/(2*tau_g_ite[g]);
    #     for(i in cpos_ite[j]:(cpos_ite[j+1]-1)){
    #       if (data2[i]<ck_ite[j]){
    #         tmp=pnorm(data2[i],mean=ck_ite[j], sd=sqrt(sig[l]));
    #         if (tmp==0){p_temp=-1250;}
    #         else{
    #           p_temp=log(tmp);           
    #         }
    #         pr_ck[g]=p_temp+se_temp; 
    #         #print(paste("p_temp (data2[i]<ck_ite[j]) =", p_temp, " i=", i));
    #       }
    #       else if (data2[i]>ck_ite[j]){
    #         tmp1=pnorm((2*ck_ite[j]-data2[i]),mean=ck_ite[j], sd=sqrt(sig[l])); 
    #         if (tmp1==0){p_temp=-1250;}
    #         else{
    #           p_temp=log(tmp1);           
    #         }         
    #         #print(paste("p_temp (data2[i]>ck_ite[j]) =", p_temp, " i=", i));
    #         pr_ck[g]=p_temp+se_temp; 
    #       }
    #       else {print(paste("error - p_temp=", p_temp));}
    #       
    #     } # for i loop
    #   }  # for g loop
    #   #print(paste("pr_ck =", pr_ck));
    #   max_prck=max(pr_ck);
    #   #print(paste("max_prck =", max_prck));
    #   prck_total=0.0;
    #   
    #   probs_ck1=exp(pr_ck[1]-max_prck);
    #   probs_ck2=exp(pr_ck[2]-max_prck);
    #   probs_ck3=exp(pr_ck[3]-max_prck);
    #   probs_ck4=exp(pr_ck[4]-max_prck);
    #   prck_total = probs_ck1+probs_ck2+probs_ck3+probs_ck4;
    #   #print(paste("prck_total =", prck_total));
    #   probs1=probs_ck1/prck_total;
    #   probs2=probs_ck2/prck_total;
    #   probs3=probs_ck3/prck_total;
    #   probs4=probs_ck4/prck_total;
    #   
    #   new_ck_probs=c(probs1, probs2, probs3, probs4);
    #   #print(paste("new_ck_probs =", new_ck_probs));
    #   
    #   new_ck=sum(new_ck_probs*mu_g_ite);
    #   #print(paste("new_ck =", new_ck, " for segment ", j));
    #   ck_ite[j]=new_ck;
    # }
    
    
    ###### Update Ck using slice sample ####
    
    density=numeric();
    new_density=numeric();
    
    for (j in 1:NoSeg){
      
      density_tmp=0;
      L_bound=0;
      R_bound=2*max(data2);
      
      for(i in cpos_ite[j]:(cpos_ite[j+1]-1)){
        se_temp=log(1/(sqrt(2*pi*sig[l])))-(data2[i]- ck_ite[j])^2/(2*sig[l]);
        density_tmp=density_tmp+se_temp;
      }
      density[j]=density_tmp+log(1/(sqrt(2*pi*tau_g_ite[gk_ite[j]])))-(ck_ite[j]- mu_g_ite[gk_ite[j]])^2/(2*tau_g_ite[gk_ite[j]]);
      rand1=log(runif(1));
      threshold=density[j]+rand1;
      #print(paste("segment j=", j)); 
      #print(paste("oldC=", ck_ite[j])); 
      #print(paste("density_tmp=", density_tmp));
      #print(paste("rand1=", rand1));
      #print(paste("density[j]=", density[j]));
      #print(paste("threshold=", threshold));   
      
      for(ii in 1:250) 
      {
        new_density_tmp=0;
        #print(paste("iteration=", ii));   
        newC = L_bound + runif(1) * (R_bound-L_bound);  #x1=L+uniform()*(R-L);
        for(i in cpos_ite[j]:(cpos_ite[j+1]-1)){
          se_temp=log(1/(sqrt(2*pi*sig[l])))-(data2[i]- newC)^2/(2*sig[l]);
          new_density_tmp=new_density_tmp+se_temp;
        }
        new_density[j]=new_density_tmp+log(1/(sqrt(2*pi*tau_g_ite[gk_ite[j]])))-(newC- mu_g_ite[gk_ite[j]])^2/(2*tau_g_ite[gk_ite[j]]);

        #print(paste("ii=", ii));    
        #print(paste("new_density[j]=", new_density[j]));
        #print(paste("newC=", newC));  
        #print(paste("new_density[j]=", new_density[j])); 
        
        if (threshold <= new_density[j]) 
        { 
          ck_ite[j] = newC; 
          #density_x = density_newx; 
          #print(paste("new C=", newC));
          break; 
        }  
        if (newC < ck_ite[j]) 
        {
          L_bound = newC;
        }
        else 
        {
        R_bound = newC;
        }
      } #for(i in 1:250) 
    } #for (j in 1:NoSeg)
    
    
    ########## Update gk ##############
    pr=numeric();
    max_pr=0;
    for (j in 1:NoSeg){
      for (g in 1:N_types){
        se_temp=log(1/(sqrt(2*pi*tau_g_ite[g])))- (ck_ite[j]- mu_g_ite[g])^2/(2*tau_g_ite[g]);
        pr[g]=log(gpi_ite[g])+se_temp; 
      }
      max_pr=max(pr);
      
      probs=numeric();
      pr_total=0;
      for (i in 1:N_types){
        probs[i]=exp(pr[i]-max_pr);
        pr_total=pr_total+probs[i];
      }
      # for debugging 
      pr1=probs[1]/pr_total;
      pr2=probs[2]/pr_total;
      pr3=probs[3]/pr_total;
      pr4=probs[4]/pr_total;
      new_pi=c(pr1, pr2, pr3, pr4);
      gk_ite[j]=choose_a_group(new_pi,N_types); 
      
      # ths is not correct - do it temperatorily for debugging
      #ck_ite[j]=box_muller_rnorm(mu_g_ite[gk_ite[j]], sqrt(tau_g_ite[gk_ite[j]]));
      
      #print(paste("new_pi =", new_pi));
      #print(paste("new group =", gk_ite[j]));
      
    }  # for (j in 1:NoSeg){
    
    
    #Update group tau^2 and group mu
    m1=mu_g_ite[1];
    m2=mu_g_ite[2];
    m3=mu_g_ite[3];
    m4=mu_g_ite[4];
    t1=tau_g_ite[1];
    t2=tau_g_ite[2];
    t3=tau_g_ite[3];
    t4=tau_g_ite[4];
    if (length(which(gk_ite==1)) > 1){
      t1=1/marsaglia_rgamma(alpha0+length(which(gk_ite==1))/2, beta0+sum((ck_ite[which(gk_ite==1)]-mu_g_ite[1])^2)/2);
      temp_mean=sum(ck_ite[which(gk_ite==1)])/length(which(gk_ite==1));
      temp_var=t1/length(which(gk_ite==1));
      #m1=rnorm(1, mean=temp_mean, sd=sqrt(temp_var));
      m1=box_muller_rnorm(temp_mean, sqrt(temp_var));
    }
    if (length(which(gk_ite==2)) > 1){
      t2=1/marsaglia_rgamma(alpha0+length(which(gk_ite==2))/2, beta0+sum((ck_ite[which(gk_ite==2)]-mu_g_ite[2])^2)/2);
      temp_mean=sum(ck_ite[which(gk_ite==2)])/length(which(gk_ite==2));
      temp_var=t2/length(which(gk_ite==2));
      #m2=rnorm(1, mean=temp_mean, sd=sqrt(temp_var));
      m2=box_muller_rnorm(temp_mean, sqrt(temp_var));
    }
    if (length(which(gk_ite==3)) > 1){
      t3=1/marsaglia_rgamma(alpha0+length(which(gk_ite==3))/2, beta0+sum((ck_ite[which(gk_ite==3)]-mu_g_ite[3])^2)/2);
      temp_mean=sum(ck_ite[which(gk_ite==3)])/length(which(gk_ite==3));
      temp_var=t3/length(which(gk_ite==3));
      #m3=rnorm(1, mean=temp_mean, sd=sqrt(temp_var));
      m3=box_muller_rnorm(temp_mean, sqrt(temp_var));
    }
    if (length(which(gk_ite==4)) > 1){
      t4=1/marsaglia_rgamma(alpha0+length(which(gk_ite==4))/2, beta0+sum((ck_ite[which(gk_ite==4)]-mu_g_ite[4])^2)/2);
      temp_mean=sum(ck_ite[which(gk_ite==4)])/length(which(gk_ite==4));
      temp_var=t4/length(which(gk_ite==4));
      #m4=rnorm(1, mean=temp_mean, sd=sqrt(temp_var));
      m4=box_muller_rnorm(temp_mean, sqrt(temp_var));
    }
    mu_g_ite=c(m1, m2, m3, m4);
    tau_g_ite=c(t1, t2, t3, t4);
    
    MUg[l, 1:N_types]<-matrix(mu_g_ite, nrow=1, ncol=4, byrow = TRUE); 
    TAUg2[l, 1:N_types]<-matrix(tau_g_ite, nrow=1, ncol=4, byrow = TRUE);
    gPI[l, 1:N_types]<-matrix(gpi_ite, nrow=1, ncol=4, byrow = TRUE);
    
    #store into global matrix
    CK[l, 1:length(ck_ite)]<-matrix(ck_ite, nrow=1, ncol=length(ck_ite), byrow = TRUE); 
    GK[l, 1:length(gk_ite)]<-matrix(gk_ite, nrow=1, ncol=length(gk_ite), byrow = TRUE); 
    CH_POS[l, 1:length(cpos_ite)]<-matrix(cpos_ite, nrow=1, ncol=length(cpos_ite), byrow = TRUE);
    Segs[l]=NoSeg;
    
    
    ####################################################
    # plot change point location and height for every iteration #
    # plot(1:1000, data2);
    # abline(v=cpos_ite, col="red"); # or CH_POS[l,]
    # for (i in 1:NoSeg){
    #     segments(x0=cpos_ite[i],y0=ck_ite[i],x1=cpos_ite[i+1],y1=ck_ite[i],col="blue");
    # }
    # 
    # Sys.sleep(0.5)
    ####################################################
    
    ######################################################
    # calculate posterior likelihood function:
    
    # post2=0.0;
    # for (j in 1:NoSeg){
    #   post1=0.0;
    #   se_temp=log(1/(sqrt(2*pi*tau_g_ite[gk_ite[j]])))-(ck_ite[j]- mu_g_ite[gk_ite[j]])^2/(2*tau_g_ite[gk_ite[j]]);
    #   for(i in cpos_ite[j]:(cpos_ite[j+1]-1)){
    #     density_x=log(1/(sqrt(2*pi*sig[l])))-(data2[i]-ck_ite[j])^2/(2*sig[l]);
    #     post1=post1+density_x;
    #   }
    #   post2=post1+post2+se_temp+log(gpi_ite[gk_ite[j]]);
    # }
    # post_dist[l]=post2+(Segs[l]-1)*log(phi[l]) + (length(data2)-Segs[l])*log(1-phi[l])+
    #   log(prior_phi);
    # 
    #post_prob[l]=log(prior_phi)+log(prior_pi)+log(prior_sig)+log(prior_tau)+log(prior_mu);
    
    # plot(1:length(post_prob), post_prob);
    

  } # while iteration 
  return(list("CK"=CK, "GK"=GK, "CH_POS"=CH_POS, "Segs"=Segs, "MUg"=MUg, "TAUg2"=TAUg2, "gPI"=gPI));
  
}  # end of function

