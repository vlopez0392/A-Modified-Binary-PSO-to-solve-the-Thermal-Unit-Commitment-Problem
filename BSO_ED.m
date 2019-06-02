function [EDC_SOL,P_srv,P_COST,tot_gen_COST,itt] = BSO_ED(Init_SOL,T,I,P_D)
%Bee swarm optimization to solve the Economic Dispatch problem (EDP)
%Assuming a quadratic I/O cost function. F(Pgi)=ai*Pgi^2+bi*Pgi+ci
%May add valve point effects if necessary 
%No loss (Can be added)

    No=size(I,1);    
    ai_or=I(:,3);
    bi_or=I(:,4);
    ci_or=I(:,5);
    
    URR=I(:,12);
    DRR=I(:,12);
    EDC_SOL=zeros(No,T);
    P_srv=zeros(1,T);
    Pgi_max_or=I(:,1); %%%Upper generation limit for each generator
    Pgi_min_or=I(:,2); %%%Lower generation limit for each generator
    
    P_COST=zeros(1,T); %%%Total generation cost
    itt=zeros(1,T);
    
    for t=1:T
       u=find(Init_SOL(:,t)>0);
       sz_u=numel(u);
       
       %%%Generator limits 
       if t==1 %%%Ignore URR/DRR when considering gen limits
           Pgi_max=Pgi_max_or;
           Pgi_min=Pgi_min_or;
       
       else    %%%Include URR/DRR only for consecutive online units 
%            Pgi_max=Pgi_max_or;
%            Pgi_min=Pgi_min_or;
%            
%            u_t=Init_SOL(:,t);
%            u_t_1=Init_SOL(:,t-1);
%            EDC_SOL_1=EDC_SOL(:,t-1);
%            
%            for k=1:N
%                dec=u_t(k)-u_t_1(k);
%                switch dec
%                    case 0
%                        if u_t(k)==1
%                              Pgi_max(k)=min(Pgi_max(k),EDC_SOL_1(k)+URR(k));
%                              Pgi_min(k)=max(Pgi_min(k),EDC_SOL_1(k)-DRR(k));
%                        end
%                        
%                    case -1 %%%Unit turn on (May be set from 0 to Pgi_min)
%                        continue
%                end  
%            end
       end
       
       %%%Load demand & cost coefficients 
       P_L=P_D(t);
       a_P=ai_or(u);
       b_P=bi_or(u); 
       c_P=sum(ci_or(u));
       d_P=zeros(sz_u,1);
       e_P=zeros(sz_u,1);
       
       P_gen=[a_P,b_P,ci_or(u),Pgi_min(u),Pgi_max(u),d_P,e_P];

        
       %%%%BSO parameter tuning 
       N=size(P_gen,1);
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

           %%%Arrays and vectors for Solution and Type of bees 
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
           for k=1:N 
               P_SOL(:,k)= (P_gen(k,4)-P_gen(k,5))*rand(D,1)+ P_gen(k,5);
           end

           %%%Demand Violation correction (FUNCTION1) %%%Include losses with beta coef
           [~,P_SOL]= P_violation(P_SOL,P_L,P_gen,N,D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAIN LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            max_iter=10;
            iter=0;
            
            %%%%%%%%%%%%%%%%%%Bee fitness calculation, sort and partition%%%%%%%%%%%%%%
            while iter<max_iter
                iter=iter+1;
                
                if iter>max_iter
                    break
                end

                %%%%Linear time varying weights 
                w_b=w_bmax-((w_bmax-w_bmin)/max_iter)*iter;

                w_g=w_gmin+((w_gmax-w_gmin)/max_iter)*iter;

                if iter==1
                    mean_prev=mean(P_SOL);
                end

                %%%%Fitness calculation 
                for k=1:D
                    OBJ_arr(k)=1/(1+a_P'*(P_SOL(k,:).^2)'+b_P'*P_SOL(k,:)'+c_P+abs(sum(d_P'*sin(e_P'*(P_gen(:,4)'-P_SOL(k,:))')))); %%%(Fitness=1/(1+obj_F))
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
                   end
                        Gbest(1,:)=best_food_EFB(1,:);
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


               [~,P_SOL]= P_violation(P_SOL,P_L,P_gen,N,D); %%%Checks for any violations and repair if necessary
                
               dff=sum(Gbest)-P_L;
            end
       end
        
       EDC_SOL(u,t)=Gbest';
       
       P_srv(t)=sum(EDC_SOL(u,t));
       
       itt(t)=iter;
       %%%%Determining total cost 
       
       P_COST(t)=ai_or'*EDC_SOL(:,t).^2+bi_or'*EDC_SOL(:,t)+c_P;
    end
    
    tot_gen_COST=sum(P_COST);
end