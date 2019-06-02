function [curr_SOL,status,trials] = constraint_repair(Init_SOL,I,Tdiff_UP,Tdiff_DW,P_D,SRREQ)
%Schedule repair function based on the pivot heuristic algorithm 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

N=size(I,1); %%%Number of TGU's
T=24;        %%%Study period duration

%%%%%%%%%%%%%%%%%%%%%%%%%%UNIT DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Generator limits
Pgi_max=I(:,1); %%%Upper generation limit for each generator

%%%I/O curve for each generator is modeled by a smooth quadratic function 
%%%i.e C(Pgi)=ai*Pgi^2+bi*Pgi+ci 
ai=I(:,3);
bi=I(:,4);
ci=I(:,5);

[~,~,~,idx_ins] = AFLC(ai,bi,ci,Pgi_max); %%%Unit priority list assuming same fuel cost for every unit

%%TIME DATA
MUT=I(:,6);     %%%Minimum up-time for each generator
MDT=I(:,7);     %%%Minimum down-time for each generator

%%START-UP AND SHUT-DOWN COSTS
init_status=I(:,11);

Ins_On_cap=zeros(3,T);     %%Insufficient online capacity
f_unit_ins_all=zeros(N,T); %%%This unit is not part of the exceptions MR,UN,Initialization
free_unit_ins=zeros(N,T);  %%%These units are offline and available to be turned on

complete=ones(1,T);        %%%Nth unit operation schedule is ON for all t

status_e=0;
status_f=0;
status_g=0;

status=[status_e,status_f,status_g];
trials=0;
max_trials=100;

    while ~all(status) 
    trials=trials+1;

    if trials>max_trials
        disp('top');
        break;
    end

    if trials==1
        %%%Insufficient online capacity initial schedule repair (RULE E)
        for t=1:T
            u_t=Init_SOL(:,t); %%Decompose into T single hour unit combinations

            Ins_On_cap(1,t)=u_t'*Pgi_max;
            Ins_On_cap(2,t)=P_D(t)+SRREQ(t);

            if Ins_On_cap(1,t)<P_D(t)+SRREQ(t)
               Ins_On_cap(3,t)=1;

               for k=1:N
                  %%%Check Exceptions 

                  if t<=Tdiff_UP(k,2)|| t<=Tdiff_DW(k,2)  %%%If the unit has been initialized (MUT/MDT) continue with next unit  
                        continue
                  else 
                      %%%Obtain available offline generators 
                      f_unit_ins_all(k,t)=k; %%Available units (Both offline and online)

                      if u_t(k)==0
                           free_unit_ins(k,t)=k; %%%Available offline units 
                      end  
                  end
               end

               %%%Repair initial schedule such that online capacity is sufficient
               free_unit_ins_tmp=find(free_unit_ins(:,t));
               [~, idx]=ismember(free_unit_ins_tmp,idx_ins);

               idx=sort(idx);
               s_idx_off=size(idx,1);

               while Ins_On_cap(3,t)==1
                    for i=1:s_idx_off
                         u_t(idx_ins(idx(i)))=1;

                         Ins_On_cap(1,t)=u_t'*Pgi_max;

                         if Ins_On_cap(1,t)<P_D(t)+SRREQ(t)
                             continue
                         else
                             Ins_On_cap(3,t)=0;
                             Init_SOL(:,t)=u_t;
                             break
                         end
                    end  

                    if Ins_On_cap(3,t)==1
                        format_ins='Insufficent capacity at hour %u. Check input data\n';
                        fprintf(format_ins,t);
                    end 

                    break
               end     
            end 
        end

        [check_F_G,~,~]=check_MUT_MDT(N,T,Init_SOL,init_status,MUT,MDT);

    else 
         prev_SOL=Init_SOL;
         v=find(check_E(1,:));
         sz_v=length(v);

         for t=1:sz_v
               v_t=v(t);
               for k=1:N
                   if check_E(1,v_t)==0
                      break
                   end

                   commit_first=idx_ins(k); %%%Commit the cheapeast generator 
                   u_t=Init_SOL(commit_first,:);

                   if u_t==complete
                      continue
                   else
                      if u_t(v_t)==0
                         u_t(v_t)=1;

                         Init_SOL(commit_first,:)=u_t;

                         check_E=check_SR_PD(T,P_D,SRREQ,Init_SOL,Pgi_max);  
                      end
                   end
               end
         end

         [check_F_G,~,~]=check_MUT_MDT(N,T,Init_SOL,init_status,MUT,MDT);   
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Check & Repair Time constraints %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if trials==1

       [Init_SOL] = sf_repair_MDT_MUT(N,T,MUT,MDT,Init_SOL,init_status,check_F_G); 

       [~,status_f,status_g]=check_MUT_MDT(N,T,Init_SOL,init_status,MUT,MDT);

       [check_E,status_e]=check_SR_PD(T,P_D,SRREQ,Init_SOL,Pgi_max);

    else %%%Pivot heuristic for MUT/MDT correction       
        
        if prev_SOL==Init_SOL
            break
            %%%%Repair directly since no pivots where found
            vio_u=zeros(N,1);
            
            for j=1:N
                  if ~all(check_F_G(j,:)==0)
                      vio_u(j)=j; 
                  end    
            end
            
            vio_u=vio_u((find(vio_u)));
            
        else
           pivot_idx=abs(prev_SOL-Init_SOL);

           vio_u=zeros(N,1); 

           for j=1:N
              if ~all(check_F_G(j,:)==0)
                  vio_u(j)=j; 
              end    
           end
           
           vio_u=vio_u((find(vio_u)));
           sz_vio_u=length(vio_u);
           
           pvt_u=zeros(sz_vio_u,1);
           
           for j=1:sz_vio_u
               k=vio_u(j);
               
               if all(pivot_idx(k,:)==0)
                   %%%%Repair directly schedule strings with no pivots
               else    
                  pvt_u(j)=k;  
               end    
           end
           
           pvt_u=pvt_u((find(pvt_u)));
           sz_pvt_u=length(pvt_u);

           %%%Initialize pivot structure
           p_s=struct('unit',zeros(sz_pvt_u,1),'pivot_no',zeros(sz_pvt_u,1),'begin_end',zeros(sz_pvt_u,6));
           p_s_type=cell(sz_pvt_u,3);

           %%%Count number of pivots and define type(Store in a structure)
           for j=1:sz_pvt_u
                k=pvt_u(j);
                p_s.unit(j)=k;

                pivot=find(pivot_idx(k,:));
                sz_pivot=length(pivot);

                diff_pivot=diff(pivot);

                if sz_pivot>1
                   if all(diff_pivot==1)
                        no_pivot=1;
                        p_s.begin_end(j,1:2)=[pivot(1),pivot(end)];
                        p_s_type{j,1}='cluster';
                   else
                       cluster=pivot(end)-pivot(1)+1;

                       if cluster<=MUT(k)
                           no_pivot=1;
                           p_s.begin_end(j,1:2)=[pivot(1),pivot(end)];
                           p_s_type{j,1}='cluster';

                       else %%%Multiple pivots in a string
                           if sz_pivot==2
                               no_pivot=2;
                               p_s.begin_end(j,1:4)=[pivot(1),pivot(1),pivot(end),pivot(end)];
                               p_s_type{j,1}='single';
                               p_s_type{j,2}='single';

                           else
                               no_pivot=0;
                               while ~isempty(pivot)
                                   no_pivot=no_pivot+1;

                                   add=pivot(1)+MUT(k);
                                   sub_pivot=pivot(pivot<=add);

                                   p_s.begin_end(j,2*no_pivot-1:2*no_pivot)=[sub_pivot(1),sub_pivot(end)];

                                   if sub_pivot(1)==sub_pivot(end)
                                      p_s_type{j,no_pivot}='single';
                                   else    
                                      p_s_type{j,no_pivot}='cluster'; 
                                   end    

                                   pivot(sub_pivot>0)=0;
                                   pivot=pivot(pivot>0);
                               end    
                           end   
                       end
                   end
                else
                    no_pivot=1;
                    p_s.begin_end(j,1)=pivot(1);
                    p_s_type{j,1}='single';
                end

                p_s.pivot_no(j)=no_pivot;
           end

           %%%%%Begin repair strategy for MUT/MDT violated hours 
           p_s_complete.pvt=p_s;
           p_s_complete.typ=p_s_type;
           search_l_r=zeros(sz_pvt_u,1);

           for k=1:sz_pvt_u
              u=p_s_complete.pvt.unit(k);
              n_sch=Init_SOL(u,:); 
              pivot_no=p_s_complete.pvt.pivot_no(k);

              %%%%Single pivot
              if pivot_no==1;
                   t=pivot_no;
                   typ=p_s_complete.typ{k,t};

                   switch typ
                       case 'single'
                          idx_begin=p_s_complete.pvt.begin_end(k);
                          idx_end=idx_begin;

                          search=abs(MUT(u)-1);

                       case 'cluster'
                          idx_begin=p_s_complete.pvt.begin_end(k,2*t-1);
                          idx_end=p_s_complete.pvt.begin_end(k,2*t);
                          D=abs(idx_end-idx_begin+1);

                          if ~all(n_sch(idx_begin:idx_end)==1)
                             n_sch(idx_begin:idx_end)=1;
                          end      

                          if D~=MUT(u) %%%%Checks left and right to make decision
                              search=abs(MUT(u)-D);
                          else
                              search=3;
                          end
                   end 
                   search_l_r(k,1)=search;

                   search_right=idx_end+search;
                   search_left= idx_begin-search;

                   if  search_right<=T &&  search_left>=1

                        idx_begin_mone=idx_begin-1;
                        idx_end_pone=idx_end+1;

                        all_right_on=all(n_sch(idx_end_pone:search_right)==1);

                        all_left_on=all(n_sch(idx_begin_mone:-1:search_left)==1);

                            if all_right_on==1 && all_left_on==1

                                dec='r';
                                n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,0,dec,idx_end,T,1);

                                dec='l';
                                n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,init_status,dec,idx_begin,0,1);

                            elseif all_right_on==1 && all_left_on==0
                                l=randi([0,1]); 

                                if l==1 
                                   z=until_zero(n_sch,idx_begin,0);
                                   n_sch(idx_begin_mone:-1:z)=1;
                                else
                                   dec='l'; 
                                   n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,init_status,dec,idx_begin,0,1);
                                end

                                dec='r';
                                n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,0,dec,idx_end,T,1); 

                            elseif all_right_on==0 && all_left_on==1     
                                r=randi([0,1]); 

                                if r==1 
                                   z=until_zero(n_sch,idx_end,1); 
                                   n_sch(idx_end_pone:1:z)=1;
                                else
                                   dec='r'; 
                                   n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,0,dec,idx_end,T,1);                            
                                end 

                                dec='l'; 
                                n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,init_status,dec,idx_begin,0,1);

                            else
                                l_r=randi([0,1]); 

                                if l_r==1
                                   z=until_zero(n_sch,idx_end,1); 
                                   n_sch(idx_end_pone:1:z)=1;

                                   dec='l'; 
                                   n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,init_status,dec,idx_begin,0,1);

                                else
                                   z=until_zero(n_sch,idx_begin,0);
                                   n_sch(idx_begin_mone:-1:z)=1;

                                   dec='r'; 
                                   n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,0,dec,idx_end,T,1);
                                end 
                            end    
                   else
                       if  search_right>T
                           n_sch(idx_end:1:T)=1;

                           dec='l';
                           n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,init_status,dec,idx_begin,0,1);  

                       elseif search_left<1
                           if init_status(u)<0 && idx_begin<MUT(u)
                              n_sch(1:MUT(u))=1;

                           elseif init_status(u)>0 && idx_begin<MDT(u)
                              n_sch(1:MUT(u))=1;

                           end

                           dec='r';
                           n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,0,dec,idx_end,T,1);

                       end
                   end

              %%%%%%Multiple pivots
              else
                   last=pivot_no;
                   n_sch=Init_SOL(u,:);
                   decomp=zeros(pivot_no,T);

                   for j=1:last %%%Decompose n_sch into no_pivot strings
                      typ=p_s_complete.typ{k,j};   
                      if j<last

                      else
                          idx_begin= p_s_complete.pvt.begin_end(k,2*j-1);
                          idx_end=   p_s_complete.pvt.begin_end(k,2*j);

                          decomp(j,idx_begin:T)=n_sch(idx_begin:T); 

                          if ~all(n_sch(idx_begin:idx_end)==1)
                             n_sch(idx_begin:idx_end)=1;
                          end    

                          if idx_end<T

                             idx_begin_mone=idx_begin-1; 

                             l_r=randi([0,1]); 

                             if l_r==1
                                z=until_zero(n_sch,idx_begin,0); 
                                n_sch(idx_begin_mone:-1:z)=1;

                                n_sch(idx_end+1:T)=0;
                             else   
                                n_sch(idx_end+1:T)=1;
                             end

                             dec='l'; 
                             n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,init_status,dec,idx_begin,0,1);
                          else

                             dec='l';  
                             n_sch=mod_repair_MDT_MUT(u,MUT,MDT,n_sch,init_status,dec,idx_begin,0,1);
                          end    
                      end    
                   end 
              end
              Init_SOL(u,:)=n_sch;
           end

           [~,status_f,status_g]=check_MUT_MDT(N,T,Init_SOL,init_status,MUT,MDT);

           [check_E,status_e]=check_SR_PD(T,P_D,SRREQ,Init_SOL,Pgi_max);
        end   
    end
       status=[status_e,status_f,status_g];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%De-committ excessive units%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%De-committment process
    for t=1:T
       SRR_PL=check_E(3,t);
       Rxs=check_E(2,t);

       uj_t=Init_SOL(:,t);
       u_on=find(uj_t);
       sz_u_on=numel(u_on);
       still_xs=true;

       while still_xs
           Rxs=0;

           for k=sz_u_on:-1:1
               uj_t(u_on(k))=0;

               Rxs=uj_t'*Pgi_max;

               if Rxs>SRR_PL
                   continue
               else
                  uj_t(u_on(k))=1; %%%When uj_t'*Pgi_max<P_L+SRR re-commit last unit
                  still_xs=false;
                  break
               end
           end    
       end

       Init_SOL(:,t)=uj_t;
    end    

    [check_F_G]=check_MUT_MDT(N,T,Init_SOL,init_status,MUT,MDT);

    %%%%%%%%%Repair the schedule in case of any further violations 

    temp_SOL=Init_SOL;
    trial=0;
    max_trial=10;

    while ~all(all(check_F_G==0))

        if trial>=max_trial
           break 
        end

        trial=trial+1;

        for k=1:N
            test_F_G=all(check_F_G(k,:)==0);

            if test_F_G==true;
                continue
            else
                if trial==1
                    n_sch=Init_SOL(k,:);
                else
                    n_sch=temp_SOL(k,:);
                end

                [intr_on,~] = count_intervals(n_sch);

                %%%%MUT violation correction (Correct first)
                if intr_on==0;
                   temp_SOL(k,:)=n_sch; 
                   continue
                else
                    int_vio_on=find(intr_on(:,3)<MUT(k));
                    sz_ivon=numel(int_vio_on);

                    for j=1:sz_ivon
                        u=int_vio_on(j);

                        idx_begin_on=intr_on(u,1);
                        idx_end_on=intr_on(u,2);

                        crr_on=MUT(k)-intr_on(u,3);

                        dec=randi([0,1]);

                        if dec==1

                            if idx_end_on+crr_on>=T
                                n_sch(idx_end_on:T)=1;
                            else
                                n_sch(idx_end_on:idx_end_on+crr_on)=1;
                            end
                        else

                            if idx_begin_on-crr_on<=1
                                n_sch(idx_begin_on:-1:1)=1;
                            else
                                n_sch(idx_begin_on:-1:idx_begin_on-crr_on)=1;
                            end                         
                        end
                    end
                end   

               [~,intr_off] = count_intervals(n_sch);

               %%%%MDT violation (Decommitt units only if possible)
               if intr_off==0;
                   temp_SOL(k,:)=n_sch; 
                   continue
               else
                   int_vio_off=find(intr_off(:,3)<MDT(k));
                   
                   if isempty(int_vio_off)
                     
                     temp_SOL(k,:)=n_sch; 
                     continue
                   else 
                       
                       sz_ivoff=numel(int_vio_off);

                       for j=1:sz_ivoff
                           u=int_vio_off(j);

                           idx_begin_off=intr_off(u,1);
                           idx_end_off=intr_off(u,2);

                           if init_status(k)<1 && idx_begin_off==1
                               continue

                           elseif idx_end_off==24
                               continue

                           else
                               if idx_end_off==idx_begin_off
                                    n_sch(idx_begin_off)=1;

                               else

                                   n_sch(idx_begin_off:idx_end_off)=1;
                               end
                           end   
                       end 
                   end
               end

                temp_SOL(k,:)=n_sch;
            end
        end

        [check_F_G,status_f,status_g]=check_MUT_MDT(N,T,temp_SOL,init_status,MUT,MDT);  
    end
    curr_SOL=temp_SOL;
    
    [~,status_e]=check_SR_PD(T,P_D,SRREQ,curr_SOL,Pgi_max);
    
    status=[status_e,status_f,status_g];
end

