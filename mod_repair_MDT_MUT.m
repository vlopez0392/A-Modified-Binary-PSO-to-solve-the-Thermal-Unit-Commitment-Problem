function [n_sch] = mod_repair_MDT_MUT(u,MUT,MDT,varargin)
%%%Modified straightforward repair for pivot heuristic
min_args=9;
max_args=9;
narginchk(min_args,max_args)
    
    n_sch=varargin{1};
    dec=varargin{3};
    pivot_no=varargin{6};

    switch pivot_no    
        case 1
            switch dec
            %%%%%%%Repair schedule to the left of idx_begin%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 'l' 
                init_status=varargin{2};
                idx_begin=varargin{4};

                if idx_begin>1
                    begin=idx_begin-1;
                else
                    begin=1;
                end    

                for t=begin:-1:1 
                   if t==begin 
                       if n_sch(t)==1 
                            u_t_1=1;   %Unit is initially online u_t-1
                            tau_MUT=1;
                            tau_MDT=0;
                       else
                            u_t_1=0;   %Unit was initially offline u_t-1
                            tau_MDT=-1;
                            tau_MUT=0;
                       end

                       if t>1
                            u_t=n_sch(t-1); %Status of current unit u_t
                       else
                            u_t=init_status(u)>0;
                       end

                   else
                      u_t_1=n_sch(t); %%idem.

                      if t>1
                        u_t=n_sch(t-1); %%idem.
                      else 
                        u_t=init_status(u)>0;
                      end
                   end   

                   %%%%Straightforward repair

                   %%%MUT
                   if u_t_1*(1-u_t)==1 %%%Unit ON->OFF
                            if tau_MUT>=MUT(u)
                                    tau_MUT=0;
                                    tau_MDT=tau_MDT-1;   
                                continue
                            else
                                if t>1
                                    n_sch(t-1)=1;
                                    tau_MDT=0;
                                    tau_MUT=tau_MUT+1;
                                else
                                    prev=n_sch(t+1:1:MUT(u)+1);

                                    if all(prev==1) && u_t>0
                                       continue
                                    else 
                                       n_sch(t:length(prev))=0;  
                                    end       
                                end

                                continue
                            end   

                   %%%MDT
                   elseif (1-u_t_1)*u_t==1 %%Unit OFF->ON

                            if -tau_MDT>=MDT(u)
                                    tau_MDT=0;
                                    tau_MUT=tau_MUT+1;  
                                continue 
                            else
                                if t>1
                                    n_sch(t-1)=0;
                                    tau_MUT=0;
                                    tau_MDT=tau_MDT-1;
                                else 
                                    prev=n_sch(t+1:1:MDT(u)+1);

                                    if all(prev==0) && u_t<1
                                       continue
                                    else 
                                        n_sch(t:length(prev))=1; 
                                    end
                                end

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

            %%%%%%%Repair schedule to the right of idx_end%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                case 'r' 

                idx_end=varargin{4};
                T=varargin{5};

                if idx_end<=23
                    n_end=idx_end+1;
                else
                    n_end=24;
                end

                for t=n_end:1:T
                   if t==n_end 
                       if n_sch(t)==1 
                            u_t_1=1;   %Unit is initially online u_t-1
                            tau_MUT=1;
                            tau_MDT=0;
                       else
                            u_t_1=0;   %Unit was initially offline u_t-1
                            tau_MDT=-1;
                            tau_MUT=0;
                       end

                       if t<T
                            u_t=n_sch(t+1); %Status of current unit u_t
                       else
                            u_t=n_sch(T);
                       end

                   else
                      u_t_1=n_sch(t); %%idem.

                      if t<T
                            u_t=n_sch(t+1); %%idem.
                      else 
                            u_t=n_sch(T);
                      end
                   end   

                   %%%%Right Straightforward repair

                   %%%MUT
                   if u_t_1*(1-u_t)==1 %%%Unit ON->OFF
                            if tau_MUT>=MUT(u)
                                    tau_MUT=0;
                                    tau_MDT=tau_MDT-1;   
                                continue
                            else
                                if t<T
                                    n_sch(t+1)=1;
                                    tau_MDT=0;
                                    tau_MUT=tau_MUT+1;
                                else
                                    continue
                                end

                                continue
                            end   

                   %%%MDT
                   elseif (1-u_t_1)*u_t==1 %%Unit OFF->ON

                            if -tau_MDT>=MDT(u)
                                    tau_MDT=0;
                                    tau_MUT=tau_MUT+1;  
                                continue 
                            else
                                if t<T
                                    n_sch(t+1)=0;
                                    tau_MUT=0;
                                    tau_MDT=tau_MDT-1;
                                else 
                                    continue
                                end

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
    end
end
