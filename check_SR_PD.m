function [check_E,status_E] = check_SR_PD(T,P_D,SRREQ,Init_SOL,PGI_MAX)
%check_SR_PD
%Checks if the current schedule satisfies the SR constraint is satisfied in the UC problem 

check_E=zeros(4,T);
    for t=1:T
       u_jt=Init_SOL(:,t);
       
       check_E(2,t)= u_jt'*PGI_MAX;   
       check_E(3,t)=P_D(t)+SRREQ(t);
       
       if check_E(2,t)< check_E(3,t)
           check_E(1,t)=true;  %Violation at time t
       else
           check_E(1,t)=false; %%No violation at time t
           
           if check_E(2,t)== check_E(3,t)
               continue
           else
               check_E(4,t)=abs(check_E(2,t)-check_E(3,t));
           end    
       end
    end
    if check_E(1,:)==false
        status_E=true;    %Coupling constraints satisfied
    else 
        status_E=false;   %Violation of coupling constraints
    end
end

