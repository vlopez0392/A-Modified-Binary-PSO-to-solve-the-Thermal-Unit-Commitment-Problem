function [P_SOL,P_srv,P_COST,tot_gen_COST,itt,Init_SOL] = F_LIM_ED_RR(Init_SOL,T,I,P_D,IDX_INS,OPTS,POZ)
%Enhanced Lambda Iteration Method to solve the Economic Dispatch problem (EDP)
%Assuming a quadratic I/O cost function. F(Pgi)=ai*Pgi^2+bi*Pgi+ci
   
    %%%%Input data
    N=size(I,1);    
    ai_or=I(:,3);
    bi_or=I(:,4);
    ci_or=I(:,5);
    
    URR=I(:,12);
    DRR=URR;
    P_SOL=zeros(N,T);
    P_srv=zeros(1,T);
    PGI_MAX=I(:,1); %%%Upper generation limit for each generator
    PGI_MIN=I(:,2); %%%Lower generation limit for each generator
    
    P_COST=zeros(1,T); %%%Total generation cost
    itt=zeros(1,T);
    
    for t=1:T
       u=find(Init_SOL(:,t)>0);
       sz_u=numel(u);
       
       %%%Load demand & cost coefficients 
       PD=P_D(t);
       ai=ai_or(u);
       bi=bi_or(u); 
       ci=ci_or(u);
       
       %%%Generator limits 
       if t==1 %%%Ignore URR/DRR at the first hour 
           Pgi_max=PGI_MAX;
           Pgi_min=PGI_MIN;
       
       else    %%%Include URR/DRR only for consecutive online units 
           Pgi_max=PGI_MAX;
           Pgi_min=PGI_MIN;
           
           u_t=Init_SOL(:,t);
           u_t_1=Init_SOL(:,t-1);
           P_SOL_1=P_SOL(:,t-1);
           
           for k=1:N
               dec=u_t(k)-u_t_1(k);
               switch dec
                   case 0
                       if u_t(k)==1
                             Pgi_max(k)=min(Pgi_max(k),P_SOL_1(k)+URR(k));
                             Pgi_min(k)=max(Pgi_min(k),P_SOL_1(k)-DRR(k));
                       end
                       
                   case -1 %%%Unit turn on (May be set from 0 to Pgi_min)
               end  
           end
           
           if sum(Pgi_max(u))<PD
               type='under';
               
               [Init_SOL,u] = Recomm_swp(I,T,t,Init_SOL,PD,Pgi_max,IDX_INS,type);
                
               u=find(u>0);
               sz_u=numel(u);
               
               ai=ai_or(u);
               bi=bi_or(u); 
               ci=ci_or(u);
               
           elseif sum(Pgi_min(u))>PD
               type='over';
               
               [Init_SOL] = Recomm_swp(I,T,t,Init_SOL,PD,Pgi_min,IDX_INS,type);
               break
           end  
           
           Pgi_max=Pgi_max(u);
           Pgi_min=Pgi_min(u);
       end
       
       it=1;
       max_iter=300; 
       
       P_S=zeros(sz_u,1);
       tol=1e-2;
       
       %%%%%Defining the initial lambda 
       lmb_v=zeros(max_iter,1);
       lmb_v(1)=(PD+sum(bi./(2*ai)))/sum(1./(2*ai));
       
       lmb=lmb_v(1);
     
       del_P_v=zeros(max_iter,1);
       del_P_v(1)=inf;
       
       del_P=del_P_v(1);
       
       %%%%%Main loop
       while abs(del_P)>tol
           
           it=it+1;
           
           if it>max_iter
               P_S_0=P_S;
               
               P_S=fmincon(@(P_S) ai'*(P_S).^2+bi'*P_S+sum(ci),P_S_0,-ones(1,sz_u),-PD,[],[],Pgi_min,Pgi_max,[],OPTS);
               del_P=PD-sum(P_S);
               
               if abs(del_P)>tol
                   disp('oops');
               end
               
               break
           end
           
           status=zeros(sz_u,1);
           for j=1:sz_u
               
                P_S(j)=(lmb-bi(j))/(2*ai(j));

                if P_S(j)>Pgi_max(j)
                    P_S(j)=Pgi_max(j);
                    status(j)=1;
                    
                elseif P_S(j)<Pgi_min(j)
                    P_S(j)=Pgi_min(j);
                    status(j)=2;
                end     
           end
           
            del_P_v(it)=PD-sum(P_S);
            
            del_P=del_P_v(it);
           
            if it==2
                if del_P_v(it)>0
                    del_lmb=0.03*lmb;

                else
                    del_lmb=-0.03*lmb;   
                end
                
                lmb_v(it)=lmb_v(it-1)+del_lmb;
                lmb=lmb_v(it);
            else
                if del_P_v(it)==del_P_v(it-1)
                    
                    if del_P_v(it)>0
                        del_lmb=0.03*lmb;

                    else
                        del_lmb=-0.03*lmb;   
                    end
                
                    lmb_v(it)=lmb_v(it-1)+del_lmb;
                    lmb=lmb_v(it);
                    continue
                else
                    lmb_v(it)=abs((lmb_v(it-1)*del_P_v(it-1)-lmb_v(it-2)*del_P_v(it))/(del_P_v(it-1)-del_P_v(it)));
                    
                    lmb=lmb_v(it);
                end
            end
       end
       itt(t)=it;
       P_SOL(u,t)=P_S;
       P_srv(t)=sum(P_SOL(u,t));
       
       %%%%Compute cost per hour
       P_COST(t)=ai_or'*P_SOL(:,t).^2+bi_or'*P_SOL(:,t)+sum(ci_or(u));
    end
        
       %%%Compute total cost
       tot_gen_COST=sum(P_COST);
 end