function [P_SOL,P_srv,P_COST,tot_gen_COST,itt] = LIM_ED_RR(Init_SOL,T,I,P_D)
%Lambda Iteration Method to solve the Economic Dispatch problem (EDP)
%Assuming a quadratic I/O cost function. F(Pgi)=ai*Pgi^2+bi*Pgi+ci
%No loss (Can be added)

    %%%%Input data
    N=size(I,1);    
    ai_or=I(:,3);
    bi_or=I(:,4);
    ci_or=I(:,5);
    
    URR=I(:,12);
    DRR=I(:,12);
    P_SOL=zeros(N,T);
    P_srv=zeros(1,T);
    Pgi_max_or=I(:,1); %%%Upper generation limit for each generator
    Pgi_min_or=I(:,2); %%%Lower generation limit for each generator
    
    P_COST=zeros(1,T); %%%Total generation cost
    itt=zeros(1,T);
    
    for t=1:T

       u=find(Init_SOL(:,t)>0);
       sz_u=numel(u);
       
       %%%Load demand & cost coefficients 
       PD=P_D(t);
       ai=ai_or(u);
       bi=bi_or(u); 
       
       %%%Generator limits 
       if t==1 %%%Ignore URR/DRR at the first hour 
           Pgi_max=Pgi_max_or;
           Pgi_min=Pgi_min_or;
       
       else    %%%Include URR/DRR only for consecutive online units 
           Pgi_max=Pgi_max_or;
           Pgi_min=Pgi_min_or;
           
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
                       continue
               end  
           end
           
           if sum(Pgi_max(u))<PD
               disp('Undercommited at');
               disp(t);
               break
               
           elseif sum(Pgi_min(u))>PD
               disp('Overcommited at');
               disp(t);
           end  
       end
       
       %%%%%Defining the initial lambda 
       lmb=(PD+sum(bi./(2*ai)))/sum(1./(2*ai));
       P_S=zeros(sz_u,1);
       del_P=inf;
       tol=1e-3;
       
       iter=0;
       max_iter=100000; 
       status=zeros(sz_u,1);
       
       while abs(del_P)>tol
           iter=iter+1;
           
           if iter>max_iter
               disp('wow');
               break
           end
           
           for i=1:sz_u
                P_S(i)=(lmb-bi(i))/(2*ai(i));

                if P_S(i)>Pgi_max(u(i))
                    P_S(i)=Pgi_max(u(i));
                    status(i)=1;

                elseif P_S(i)<Pgi_min(u(i))
                    P_S(i)=Pgi_min(u(i));
                    status(i)=2;
                end
           end
          
           del_P=PD-sum(P_S);
           
           del_lmb=del_P/sum(1./(2*ai));
           lmb=lmb+del_lmb;
       end
       itt(t)=iter;
       
       P_SOL(u,t)=P_S;
       
       P_srv(t)=sum(P_SOL(u,t));
       
       %%%%Determining total cost 
       
       P_COST(t)=ai_or'*P_SOL(:,t).^2+bi_or'*P_SOL(:,t)+sum(ci_or(u));
    end
    
      tot_gen_COST=sum(P_COST);
end
