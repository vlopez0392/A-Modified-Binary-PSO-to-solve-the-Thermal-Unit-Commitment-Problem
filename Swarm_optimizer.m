%%%%%%%%%%%%%%%Optimization of the TUC problem by VS-BPSO%%%%%%%%%%%%%%%%%%
%%%%Coded by Vick @NSYSU%%%

clear;
tic
run('Swarm_generator.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%OPTIMIZATION BY VS-BPSO%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Swarm size
POPN=SWARM_SZ; %%%POP_N provided by swarm_generator.m

%%%Particles in the swarm 
swarm=population;

%%%%Parameters for the BPSO 
VMAX=2;
VMIN=-2;
WMAX=0.9;
WMIN=0.4;
C1=2;
C2=2;

v=zeros(POPN,N*T); %%%Velocity vector 
m=zeros(POPN,N*T); %%%Mapping function vector
x=zeros(POPN,N*T); %%%Position vector

MAX_K=100;
[F,idx_gb]=sort(total_COST);

conv_beh=zeros(1,MAX_K);

for i=1:MAX_K
    
    %%%Inertia weight update
    w=WMAX-i*((WMAX-WMIN)/MAX_K);
    
    if i==1
       x=swarm;      %%% Initial position is the population by swarm generator
       pbest=swarm;  %%% Previous best position is initial position 
       gbest=swarm(idx_gb(1),:); %%%Global best position in current swarm
       
       fit=1./F;
       pbest_v=fit;
       gbest_v=fit(idx_gb(1));
    else
       
        for m=1:POPN 
           if fit(m)>pbest_v(m)
               pbest(m,:)=x(m,:);
               
               pbest_v(m)=fit(m);
           end 
        end
        
        [gbest_tmp_v,idx_tmp]=max(fit);
        
        if gbest_tmp_v>gbest_v
            
           gbest_v= gbest_tmp_v;
           
           gbest=x(idx_tmp(1),:);
        end
    end
    
    for j=1:POPN 
       for k=1:N*T

           v(j,k)=w*v(j,k)+C1*rand*(pbest(j,k)-x(j,k))+C2*rand*(gbest(k)-x(j,k));
           
           if v(j,k)>VMAX
              v(j,k)=VMAX;
              
           elseif v(j,k)<VMIN
              v(j,k)=VMIN;
           end     

           m(j,k)=abs(tanh(v(j,k))); %%%Map continuous velocity values into probability values
        
           if m(j,k)>rand
               x(j,k)=~x(j,k);
           else
               x(j,k)=x(j,k);
           end
       end

       Init_SOL=reshape(x(j,:),N,T);

       [curr_SOL,status,trials] = constraint_repair(Init_SOL,I,Tdiff_UP,Tdiff_DW,P_D,SRREQ);

       x(j,:)=reshape(curr_SOL,1,N*T);
    
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       %%%%%%%%%%%%%%%%%%Economic-Dispatch%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

       %%%%Base case

       [P_SOL,P_srv,P_COST,tot_gen_COST,itt] = F_LIM_ED(curr_SOL,T,I,P_D,OPTS);

       %%%Ramping Rate constraints included
 
%        [P_SOL,P_srv,P_COST,tot_gen_COST,itt,curr_SOL] = F_LIM_ED_RR(curr_SOL,T,I,P_D,IDX_INS,OPTS,POZ);

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
       %%%%%%%%%%%%%%%%%%Start-up-cost computation%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
        [su_COST,tot_su_COST] = SU_COST(N,curr_SOL,init_status,SU_H,SU_C,CSH,MUT,MDT);
       
        fit(j)=1/(tot_gen_COST+tot_su_COST);          
    end
    
    conv_beh(i)=(1/gbest_v);
    disp(i);
    disp(gbest_v);
    
end    
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RESULTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%Optimal Solution Computation%%%%%%%%%%%%%%%%%%%%%%

%%%Optimal Schedule
OPT_SOL=reshape(gbest,N,T)';

%%%Economic Dispatch
% %Ramping rates and POZ neglected
[P_SOL_OPT,P_srv_opt,P_COST_opt,tot_gen_COST_opt,itt_opt] = F_LIM_ED(OPT_SOL',T,I,P_D,OPTS);

%Ramping Rate and POZ constraints included
% [P_SOL_OPT,P_srv_opt,P_COST_opt,tot_gen_COST_opt,itt_opt,~] = F_LIM_ED_RR(OPT_SOL',T,I,P_D,IDX_INS,OPTS,POZ);

%%%Startup-cost

[su_COST_opt,tot_su_COST_opt,~] = SU_COST(N,OPT_SOL',init_status,SU_H,SU_C,CSH,MUT,MDT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting the results%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(conv_beh(:),'k-o');
title('Gbest vs iteration');
xlabel('Iteration');
ylabel('Gbest');


%%%%%%%%%%%%%%Velocity mapping functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Copy & paste as required%%%%

%%%V-SHAPED

%m(j,k)=abs(erf(v(j,k)*((pi)^(1/2)/2)));  %%V1

%m(j,k)=abs(v(j,k)/(1+((v(j,k))^2)^(1/2))); %%%V2

% m(j,k)=abs(tanh(v(j,k))); %%%V2

% m(j,k)=abs((2/pi)*atan(v(j,k)*((pi)^(1/2)/2))); %%%V4


