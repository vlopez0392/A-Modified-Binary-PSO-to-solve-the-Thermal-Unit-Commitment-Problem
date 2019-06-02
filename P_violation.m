function [PV,P_SOL] = P_violation(P_SOL,P_D,P_gen,N,D)
%Computes power violation and adjusts power outputs of each generator 
PV=zeros(D,1);

%%%%Verify load can be satisfied by system 
system_limit=sum(P_gen(:,5));

if P_D>system_limit
    error('Load cannot be satisfied by current system configuration');
else
    %%%%%%Verify that each of the outputs in P_SOL are within limits and correct if else
    for i=1:N
        UPPER=P_SOL(:,i)<=P_gen(i,5);
        idx_up=find(UPPER~=1);

        if isempty(idx_up)==false
            for k=1:numel(idx_up)
               P_SOL(idx_up(k),i)=P_gen(i,5);
            end    
        end    

        LOWER=P_SOL(:,i)>=P_gen(i,4); 
        idx_low=find(LOWER~=1);

        if isempty(idx_low)==false
            for k=1:numel(idx_low)
               P_SOL(idx_low(k),i)=P_gen(i,4);
            end    
        end
    end  
    
    %%%%%%Demand correction taking only a feasible solution
    for k=1:D
         PV(k)=sum(P_SOL(k,:))-P_D; % Initial Violation of power demand 

         if PV(k)==0 %%%Very unlikely however possible
             continue 
         else
             unit_selec=1:N;
             while PV(k)~=0; 
                 randi_unit=unit_selec(randi(length(unit_selec)));            

                 P_gen_rand=P_SOL(k,randi_unit)-PV(k);

                 Pgmin=P_gen(randi_unit,4); %%%Get lower limit for the generator
                 Pgmax=P_gen(randi_unit,5); %%%Get upper limit for the generator

                 if P_gen_rand<Pgmin

                    P_gen_rand=Pgmin;

                 elseif P_gen_rand>Pgmax

                    P_gen_rand=Pgmax;  
                 end

                 P_SOL(k,randi_unit)=P_gen_rand; 

                 PV(k)=sum(P_SOL(k,:))-P_D; %%%%Recalculate violation 

                 if PV(k)==0 %%%If there is no violation continue 
                     continue 
                 else
                    idx=find(randi_unit);
                    unit_selec(idx)=[]; %#ok<FNDSB> %%Elimininate that TGU from next selection

                    if numel(unit_selec)==0 %%Retry if unit_selec becomes empty 
                        unit_selec=1:N;
                    end

                 end 
             end
         end
    end 
     
    check=(sum(sum(P_SOL))/P_D)-D;
    eps=1e-4;
    
    if check>eps
        error('Solution is unfeasible');
    else
    end
end

