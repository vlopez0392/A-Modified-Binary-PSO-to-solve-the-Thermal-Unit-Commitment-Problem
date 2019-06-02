function [ P_SOL,P_srv,P_COST,tot_gen_COST,itt] = ED_fmincon(Init_SOL,T,I,P_D,OPTS)
%Find the minimum of a constrained non-linear optimization (fmincon)
%Assuming a quadratic I/O cost function. F(Pgi)=ai*Pgi^2+bi*Pgi+ci
    
    %%%%Input data
    N=size(I,1);    
    ai_or=I(:,3);
    bi_or=I(:,4);
    ci_or=I(:,5);
    URR=I(:,12);
    DRR=I(:,12);
    
    P_SOL=zeros(N,T);
    P_srv=zeros(1,T);
    
    %%%Generator limits 
    %%%Ignore URR/DRR when considering gen limits
    Pgi_max_or=I(:,1); %%%Upper generation limit for each generator
    Pgi_min_or=I(:,2); %%%Lower generation limit for each generator
    
    P_COST=zeros(1,T); %%%Total generation cost
    itt=zeros(1,T);
    
    for t=1:T
       u=find(Init_SOL(:,t)>0);
       sz_u=numel(u);
       
       P_S_0=zeros(sz_u,1);
       
       %%%Load demand & cost coefficients 
       PD=P_D(t);
       ai=ai_or(u);
       bi=bi_or(u); 
       ci=ci_or(u);
       
       %%%Generator limits 
       if t==1 %%%Ignore URR/DRR at the first hour 
           Pgi_max=Pgi_max_or(u);
           Pgi_min=Pgi_min_or(u);
       
       else    %%%Include URR/DRR only for consecutive online units 
%            Pgi_max=Pgi_max_or;
%            Pgi_min=Pgi_min_or;
%            
%            u_t=Init_SOL(:,t);
%            u_t_1=Init_SOL(:,t-1);
%            P_SOL_1=P_SOL(:,t-1);
%            
%            for k=1:N
%                dec=u_t(k)-u_t_1(k);
%                switch dec
%                    case 0
%                        if u_t(k)==1
%                              Pgi_max(k)=min(Pgi_max(k),P_SOL_1(k)+URR(k));
%                              Pgi_min(k)=max(Pgi_min(k),P_SOL_1(k)-DRR(k));
%                        end
%                        
%                    case -1 %%%Unit turn on (May be set from 0 to Pgi_min)
%                end  
%            end
%            
%            Pgi_min=Pgi_min(u);
%            Pgi_max=Pgi_max(u);
%            
%            if sum(Pgi_max)<PD
%                disp('Undercommited at');
%                disp(t);
%                break
%                
%            elseif sum(Pgi_min)>PD
%                disp('Overcommited at');
%                disp(t);
%            end  
       end
       Pgi_max=Pgi_max_or(u);
       Pgi_min=Pgi_min_or(u);
       
       %%%%%Defining the initial lambda & initial guess for P_S 
       lmb_xo=(PD+sum(bi./(2*ai)))/sum(1./(2*ai));
       
       for j=1:sz_u
            P_S_0(j)=(lmb_xo-bi(j))/(2*ai(j));
       end
       
       P_S=fmincon(@(P_S) ai'*(P_S).^2+bi'*P_S+sum(ci),P_S_0,-ones(1,sz_u),-PD,[],[],Pgi_min,Pgi_max,[],OPTS);
       
       P_SOL(u,t)=P_S;
       
       P_SOL(u,t)=P_S;
       P_srv(t)=sum(P_SOL(u,t));
       
       %%%%Compute total cost 
       
       P_COST(t)=ai_or'*P_SOL(:,t).^2+bi_or'*P_SOL(:,t)+sum(ci_or(u));
    end
    
      tot_gen_COST=sum(P_COST);
end

