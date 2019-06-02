%%%%Classical Economic Dispatch problem solved by Bee swarm optimization (BSO) 
%%%%Based on BSO by Akbari%%%%
%%%%Coded by Vick @NSYSU%%%%%
clear;
clc;
tic

%%%Based on the quadratic I/O cost function description of a generator with
%%%Valve point effect coefficients 
%%%% F(Pgi)=ai*Pgi^2+bi*Pgi+ci+|di*sin(ei*(Pgi_min-Pgi))|

%%%SAMPLE TEST I unit characteristics                   %V.P
%            %ai    bi      ci      %PGi min  %Pgi max  %di    %ei
P_gen= [0.001562    7.92    561     150        600      0      0     ;
        0.001940    7.85    310     100        400      0      0     ;
        0.004820    7.97    78      50         200      0      0     ];
%     
%%%SAMPLE TEST II unit characteristics                   %V.P
%            %ai    bi      ci      %PGi min  %Pgi max  %di    %ei
% P_gen= [0.001280    6.48    459     150        600      0      0     ;
%         0.001940    7.85    310     100        400      0      0     ;
%         0.004820    7.97    78      50         200      0      0     ];    

%%%SAMPLE TEST III unit characteristics    
%           %ai   bi      ci  %PGi min  %Pgi max  %di    %ei
% P_gen=[0.0063   5.32    500   100       600       0      0   ; 
%        0.0413   5.13    400   150       700       0      0   ;
%        0.0093   4.83    450   50        400       0      0   ;  
%        0.0073   4.61    550   50        200       0      0  ]; 

%%%SAMPLE TEST IV unit characteristics (Considering valve point effects)

%            %ai    bi      ci   %PGi min  %Pgi max      %di    %ei
% P_gen=  [0.001562  7.92     561	 100	    600      300	0.0315;
%         0.00194   7.85	    310	 100	    400 	 200	0.042;
%         0.00482   7.97	    78   50	        200      150	0.063];

%            %ai    bi      ci   %PGi min  %Pgi max      
% P_gen=  [0.00048  16.19  1000	 150	    455    0    0;                
%          0.00031  17.26  970	 150	    455    0    0]; 	 

a_P=P_gen(:,1);
b_P=P_gen(:,2);
c_P=sum(P_gen(:,3));
d_P=P_gen(:,6);
e_P=P_gen(:,7);
P_L=700; %%%Demand at time t

 %%%%BSO parameter tuning 
 N=size(P_gen,1); %%%Generator number  
 D=12*N; %%%Population number 
 E_O_pr=0.47; %%%Percent of experienced forager & onlooker bees 
 S_pr=0.06; %%%Percent of scout bees
 w_bmin=0.75;
 w_bmax=1.5;
 w_gmin=0.75;
 w_gmax=1.5;
 
 
 if E_O_pr*2+S_pr~=1
     disp('Wrong assignment of bee percentages, please retry.');
 else
 
 %%%%Arrays and vectors for Solution and Type of bees 
 P_SOL=zeros(D,N);
 OBJ_arr=zeros(D,1);
 
 EFB=zeros(floor(E_O_pr*D),1); %%%Experienced forager bee set
 OLK=zeros(floor(E_O_pr*D),1); %%%Onlooker bee set 
 SCT=zeros(D-sum(numel(EFB)*2),1); %%%Scout bee set
 
 s_EFB=size(EFB,1);
 s_OLK=size(OLK,1);
 s_SCT=size(SCT,1);
 
 best_food_EFB=zeros(size(EFB,1),N); %%%Cognitive component of EFB
 Gbest=zeros(1,N);                   %%%Global best solution found so far 
 
 %%%%Initialize the bees 
 
 %%%%Solution is initialized within (Pgmin,Pgmax) 
 %%%In the UC problem after predispatch this must be modified to include
 %%%UR/DR 

 for k=1:N 
     P_SOL(:,k)= (P_gen(k,4)-P_gen(k,5))*rand(D,1)+ P_gen(k,5);
 end
 
%%%Demand Violation correction (FUNCTION1) %%%Include losses with beta coef
[PV,P_SOL]= P_violation(P_SOL,P_L,P_gen,N,D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAIN LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_iter=100;
iter=0; 
conv_beh=zeros(max_iter,N+1); %%%Convergent behavior of best solution 
mean_beh=zeros(max_iter,1); %%%Convergent behavior of solution mean

%%%%%%%%%%%%%%%%%%Bee fitness calculation, sort and partition%%%%%%%%%%%%%%
while iter<max_iter
    iter=iter+1;
    
    %%%%Linear time varying weights 
    w_b=w_bmax-((w_bmax-w_bmin)/max_iter)*iter;
    
    w_g=w_gmin+((w_gmax-w_gmin)/max_iter)*iter;
    
    if iter==1
        mean_prev=mean(P_SOL);
    end
    
    %%%%Fitness calculation 
    for k=1:D
        OBJ_arr(k)=1/(1+a_P'*(P_SOL(k,:).^2)'+b_P'*P_SOL(k,:)'+c_P'+abs(sum(d_P'*sin(e_P'*(P_gen(:,4)'-P_SOL(k,:))')))); %%%(Fitness=1/(1+obj_F))
    end

    %%%Sorting
    [OBJ_arr,idx_or]=sort(OBJ_arr,'descend'); %%Fitness is sorted in an descending fashion

    %%%Partition 
    EFB=OBJ_arr(1:numel(EFB));
    OLK=OBJ_arr(numel(OLK)+1:D-numel(SCT));
    SCT=OBJ_arr(D-numel(SCT)+1:end);
    
    idx_EFB=idx_or(1:numel(EFB));
    idx_OLK=idx_or(numel(OLK)+1:D-numel(SCT));
    idx_SCT=idx_or(D-numel(SCT)+1:end);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%EXPERIENCED FORAGER BEES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if iter==1  
       for i=1:s_EFB
            best_food_EFB(i,:)=P_SOL(idx_or(i),:);
            Gbest(1,:)=best_food_EFB(1,:);

        end
    else
        fit_x= fit_EFB(iter,EFB,a_P,b_P,c_P,d_P,e_P,P_gen,P_SOL,idx_or);
        fit_b= fit_EFB(iter,EFB,a_P,b_P,c_P,d_P,e_P,P_gen,best_food_EFB);
        fit_g= fit_EFB(iter,EFB,a_P,b_P,c_P,d_P,e_P,P_gen,Gbest);

        for i=1:s_EFB
            %%%Cognitive knowledge
            if fit_x(i)> fit_b(i)
               best_food_EFB(i,:)=P_SOL(idx_or(i),:);
            end

            %%%Social knowledge
            if fit_b(i)> fit_g
                Gbest= best_food_EFB(i,:);
            end
        end
    end
    
    %%%Position update for the EFB
    for i=1:s_EFB
        r_e=randi(1);
        r_b=randi(1);
        for d=1:N
            P_SOL(idx_or(i),d)=P_SOL(idx_or(i),d)+w_b*r_b*(best_food_EFB(i,d)-P_SOL(idx_or(i),d))+w_g*r_e*(Gbest(d)-P_SOL(idx_or(i),d));
        end
    end
    
    [PV,P_SOL]= P_violation(P_SOL,P_L,P_gen,N,D); %#ok<ASGLU> %%%Checks for any violations and repair if necessary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%ONLOOKER BEES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:s_OLK
        %%%Onlooker elite bee selection by roulette wheel 
        prob_v=sort(EFB./sum(EFB),'ascend'); 
        idx_EFB_rev=flipud(idx_EFB);
        q_prob_v=zeros(s_EFB,1);

        for j=1:s_EFB
            if j==1
                q_prob_v(j)=prob_v(1);
            else
                q_prob_v(j)=sum(prob_v(1:j-1))+prob_v(j);
            end
        end
        r=rand; 
        q=find(q_prob_v>=r);
        q=q(1); 
        elite_bee=idx_EFB_rev(q); 

        %%Position update for the OLK
        for k=1:N
            P_SOL(idx_OLK(i),k)=P_SOL(idx_OLK(i),k)+w_g*rand*(P_SOL(idx_EFB(q),k)-P_SOL(idx_OLK(i),k));
        end
    end
    
    [PV,P_SOL]= P_violation(P_SOL,P_L,P_gen,N,D); %#ok<ASGLU> %%%Checks for any violations and repair if necessary
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SCOUT BEES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    l=round(1+rand); %%Obedience factor 
    
    %Position update for the SCT
    for i=1:s_SCT
        for k=1:N
            P_SOL(idx_SCT(i),k)=P_SOL(idx_SCT(i),k)+rand*(Gbest(k)-l*mean_prev(k));
        end
    end
    

   [PV,P_SOL]= P_violation(P_SOL,P_L,P_gen,N,D); %%%Checks for any violations and repair if necessary
   
   mean_prev=mean(P_SOL); %%%%Mean value of the population to be used in the
                           %%Next iteration 
    

   conv_beh(iter,:)=[1./fit_EFB(iter,EFB,a_P,b_P,c_P,d_P,e_P,P_gen,Gbest)-1,Gbest];
   mean_beh(iter,:)=1/mean(OBJ_arr);

end

disp('%%%%Classical Economic Dispatch problem solved by Bee Swarm Optimization(BSO)'); 
fprintf('\n')
toc
fprintf('\n')

%%%%Display resultsin workspace
disp('The optimal power output for each generator is:'); 
fprintf('\n')
disp(conv_beh(end,2:end)); 
fprintf('\n')
disp('The minimum production cost in $'); 
fprintf('\n')
disp(conv_beh(end,1)); 
fprintf('\n')
disp('Demand sufficed by optimal generation output?'); 
    Lcond=sum(conv_beh(end,2:end));
    if Lcond>=P_L
        fprintf('\n')
        disp('Yes, The current satisfied load in MW is'); 
        fprintf('\n')
        disp(Lcond); 
    else
        disp('The current demmand is not satisfied exacly'); 
        disp(Lcond); 
        eps=1e-3;
        fprintf('\n')
        if Lcond-P_L<eps
            disp('However MW Difference is neglible '); 
            disp(Lcond-P_L);
        end   
    end
    
%%%%Plot results%%%%
figure(1)
plot(conv_beh(:,1),'b-x');
title('Gbest vs iteration');
xlabel('Iteration');
ylabel('Gbest');
figure(2)
plot(mean_beh,'g-x');
hold off
end






































