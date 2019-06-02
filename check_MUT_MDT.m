function [check_F_G,status_f,status_g]=check_MUT_MDT(N,T,Init_SOL,init_status,MUT,MDT)
%Check_MUT_MDT
%Checks if the current schedule satisfies the MUT/MDT constraints are satisfied in the UC problem 
%status_f=true when no violations to the MUT constraint occur
%status_g=true when no violations to the MDT constraint occur
COM=ones(1,T);
check_F_G=ones(N,2);
check_units=zeros(N,2); %%%Exclude fast starting units from being checked 
    for k=1:N
        if MUT(k)>1
            check_units(k,1)=1;
        else
            check_F_G(k,1)=false;
        end
        
        if MDT(k)>1
            check_units(k,2)=1;
        else
            check_F_G(k,2)=false;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%Check MUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:N
       n_sch=Init_SOL(k,:);
       if n_sch==COM 
           check_F_G(k,:)=[false,false]; %%%No violation if all units ON
           continue
       elseif n_sch==zeros(1,T)
           check_F_G(k,:)=[false,false]; %%%No violation if all units OFF
           continue
       else 
           if check_units(k,1)==0
               continue
           else
                MUT_check=find(n_sch);
                int_UP=diff(MUT_check);
                sz_int_UP=length(int_UP);
                sz_MUT_check=length(MUT_check);
               if all(int_UP==1) &&  sz_int_UP+1>=MUT(k) %%%One single interval in the schedule that meets the MUT requirements 
                    check_F_G(k,1)=false;
                    continue
               else                                      %%%Several intervals (MUT requirements must be checked)
                    counter_up=1;
                    intervals_UP=zeros(1,sz_MUT_check);
                    if sz_MUT_check>1;
                        for j=1:sz_MUT_check
                            if j<sz_MUT_check
                                el_t_1=MUT_check(j);
                                el_t=MUT_check(j+1);
                            else
                                el_t_1=MUT_check(j-1);
                                el_t=MUT_check(j);
                            end
                            if abs(el_t-el_t_1)==1
                                counter_up=counter_up+1;

                                if j==sz_MUT_check
                                    intervals_UP(j)=counter_up-1;
                                end     
                            else
                                intervals_UP(j)=counter_up;
                                counter_up=1;
                            end    
                        end
                    else  
                        if MUT_check==T
                            check_F_G(k,1)=false;

                        elseif MUT_check==1; %%%Only when init_status==MUT
                            check_F_G(k,1)=false;    

                        else
                            check_F_G(k,1)=true;   
                        end
                           continue
                    end   
               end  
                idx_UP_final_hr=MUT_check(intervals_UP>0);    %%%Last hour in the nth interval
                hours_UP=intervals_UP(intervals_UP>0);        %%%No of hrs UP in the nth interval
                index_interval=[(idx_UP_final_hr-hours_UP+1)',idx_UP_final_hr']; %%Hours UP index in n_sch             
                no_intervals=numel(idx_UP_final_hr); 
                holds_MUT=0;
                for n=1:no_intervals        %%%Check that interval length satisfy MUT constraint 
                    if hours_UP(n)>=MUT(k)  %%%Discard intervals that meet MUT requirement right away
                        holds_MUT=holds_MUT+1;
                    else                    %%%Check if hours UP occurs at beggining or end of current schedule
                        if index_interval(n,1)==1
                            if init_status(k)>0
                                holds_MUT=holds_MUT+1;
                            end    
                            
                        elseif index_interval(n,2)==24
                            holds_MUT=holds_MUT+1; 
                        end
                    end
                end
                if holds_MUT==no_intervals %%%All UP intervals meet the MUT
                     check_F_G(k,1)=false;
                else
                    continue
                end
           end
       end
    end 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%Check MDT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:N
       if  check_F_G(k,:)==false;
           continue
       end
        
       n_sch=Init_SOL(k,:);
        
       if check_units(k,2)==0
           continue
       else 
            MDT_check=find(n_sch==0);
            int_DW=diff(MDT_check);

            sz_int_DW=length(int_DW);
            sz_MDT_check=length(MDT_check);

                if all(int_DW==1) && sz_int_DW+1>=MDT(k) %%%One single interval that meets the MDT requirements 
                    check_F_G(k,2)=false;
                    continue
                else                                     %%%Several intervals (MDT requirements must be checked)
                    counter_dw=1;
                    intervals_DW=zeros(1,sz_MDT_check);

                    if sz_MDT_check>1
                        for j=1:sz_MDT_check
                            if j<sz_MDT_check
                                el_t_1=MDT_check(j);
                                el_t=MDT_check(j+1);

                            else
                                el_t_1=MDT_check(j-1);
                                el_t=MDT_check(j);
                            end

                            if abs(el_t-el_t_1)==1
                                counter_dw=counter_dw+1;

                                if j==sz_MDT_check
                                    intervals_DW(j)=counter_dw-1;
                                end     
                            else
                                intervals_DW(j)=counter_dw;
                                counter_dw=1;
                            end    
                        end
                    else 
                        intervals_DW=counter_dw;
                    end    
                end

                idx_DW_final_hr=MDT_check(intervals_DW>0);    %%%Last hour in the nth interval
                hours_DW=intervals_DW(intervals_DW>0);        %%%No of hrs UP in the nth interval   

                index_interval=[(idx_DW_final_hr-hours_DW+1)',idx_DW_final_hr']; %%Hours UP index in n_sch             
                no_intervals=numel(idx_DW_final_hr); 

                holds_MDT=0;

                for n=1:no_intervals        %%%Check that interval length satisfy MUT constraint 
                    if hours_DW(n)>=MDT(k)  %%%Discard intervals that meet MDT requirement right away
                        holds_MDT=holds_MDT+1;

                    else                    %%%Check if hours DW occurs at beggining or end of current schedule
                        if index_interval(n,1)==1
                            if init_status(k)<0
                                holds_MDT=holds_MDT+1;
                            end    

                        elseif index_interval(n,2)==24
                            holds_MDT=holds_MDT+1; 
                        end
                    end
                end

                if holds_MDT==no_intervals %%%All DW intervals meet the MUT
                     check_F_G(k,2)=false;
                else
                    continue
                end
        end   
    end      
    
    if check_F_G(:,1)==zeros(N,1)
        status_f=true;    %MUT constraint satisfied
        
    else 
        status_f=false;   %Violation of MUT
    end
    
    if check_F_G(:,2)==zeros(N,1)
        status_g=true;    %MDT constraint satisfied

    else 
        status_g=false;   %Violation of MDT
    end    
end


