function [Init_SOL] = sf_repair_MDT_MUT(N,T,MUT,MDT,Init_SOL,init_status,check_F_G)
%%%Straightforward repair of MDT/MUT constraint in the Thermal UC problem 

    for k=1:N
        n_sch=Init_SOL(k,:);   %%%Nth unit operation schedule

        if check_F_G(k,:)==0; %%%If all MDT/MUT requirements are met continue repairing next's unit schedule
            continue

        else     
            for t=1:T 
               if t==1 
                   if init_status(k)>0
                        u_t_1=1; %Unit was initially online u_t-1
                        tau_MUT=init_status(k);
                        tau_MDT=0;
                   else
                        u_t_1=0; %Unit was initially offline u_t-1
                        tau_MDT=init_status(k);
                        tau_MUT=0;
                   end
                        u_t=n_sch(1,t); %Status of current unit u_t
               else    
                  u_t_1=n_sch(1,t-1); %%idem.
                  u_t=n_sch(1,t);     %%idem.
               end   

               %%%%Straightforward repair

               %%%MUT
               if u_t_1*(1-u_t)==1 %%%Unit ON->OFF
                        if tau_MUT>=MUT(k)
                                tau_MUT=0;
                                tau_MDT=tau_MDT-1;   
                            continue
                        else   
                            n_sch(t)=1;
                            tau_MDT=0;
                            tau_MUT=tau_MUT+1;

                        end   

               %%%MDT
               elseif (1-u_t_1)*u_t==1 %%Unit OFF->ON

                        if -tau_MDT>=MDT(k)
                                tau_MDT=0;
                                tau_MUT=tau_MUT+1;  
                            continue 
                        else
                            n_sch(t)=0;
                            tau_MUT=0;
                            tau_MDT=tau_MDT-1;
                            continue
                        end                  

               else %%No change in unit status ON->ON or OFF->OFF
                   if u_t_1==1 && u_t==1
                        tau_MUT=tau_MUT+1;
                        continue
                   else
                        tau_MDT=tau_MDT-1;
                   end
               end 
            end
        end
        Init_SOL(k,:)=n_sch;
    end
end

