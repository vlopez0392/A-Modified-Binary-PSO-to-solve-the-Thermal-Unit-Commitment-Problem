function [P_SOL,P_srv,P_COST,tot_gen_COST,itt] = F_LIM_ED(Init_SOL,T,I,P_D,OPTS)
%Enhanced Lambda Iteration Method to solve the Economic Dispatch problem (EDP)
%Assuming a quadratic I/O cost function. F(Pgi)=ai*Pgi^2+bi*Pgi+ci
%No loss (Can be added)

    %%%%Input data
    N=size(I,1);    
    ai_or=I(:,3);
    bi_or=I(:,4);
    ci_or=I(:,5);
   
    P_SOL=zeros(N,T);
    P_srv=zeros(1,T);
    
    %%%Generator limits 
    %%%Ignore URR/DRR when considering gen limits
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
      
       it=1;
       max_iter=3000; 
       
       P_S=zeros(sz_u,1);
       tol=1e-2;
       
       %%%Generation limits 
       Pgi_max=PGI_MAX(u);
       Pgi_min=PGI_MIN(u);
       
       %%%%%Defining the initial lambda 
       lmb_v=zeros(max_iter,1);
       lmb_v(1)=(PD+sum(bi./(2*ai)))/sum(1./(2*ai));
       
       lmb=lmb_v(1);
       
       del_P_v=zeros(max_iter,1);
       del_P_v(1)=inf;
       
       del_P=del_P_v(1);
       status=zeros(sz_u,1);
       
       while abs(del_P)>tol
           
           it=it+1;
           
           if it>max_iter
               P_S_0=P_S;
               
               P_S=fmincon(@(P_S) ai'*(P_S).^2+bi'*P_S+sum(ci),P_S_0,-ones(1,sz_u),-PD,[],[],PGI_MIN(u),PGI_MAX(u),[],OPTS);
               del_P=PD-sum(P_S);
               
               if abs(del_P)>tol
                   disp('oops');
               end
               
               break
           end
           
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
      
      %%%%Compute total cost 
      tot_gen_COST=sum(P_COST);
 end