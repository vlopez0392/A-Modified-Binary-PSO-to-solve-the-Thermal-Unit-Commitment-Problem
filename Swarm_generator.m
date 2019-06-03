%%%%%%Unit commitment schedule generator for swarm based metaheuristics%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Coded by Vick @NSYSU%%%
clc;
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST CASE I-II-V-VI%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%IEEE 10 UNIT-SYSTEM DATA

% %Pmax %%Pmin %%ai   %%bi   %%ci  %MUT %MDT %Shot %Scold %CSH %I.S %%D/URR  %%%MR %%%UN                         
I=[ 455  150 0.00048  16.19  1000  8    8    4500  9000   5     8   91        0       0    ;    %U1
    455  150 0.00031  17.26  970   8    8    5000  10000  5     8   91        0       0    ;    %U2
    130  20    0.002  16.60  700   5    5    550   1100   4    -5   26        0       0    ;    %U3
    130  20  0.00211  16.50  680   5    5    560   1120   4    -5   26        0       0    ;    %U4 
    162  25  0.00398  19.70  450   6    6    900   1800   4    -6   33        0       0    ;    %U5
     80  20  0.00712  22.26  370   3    3    170   340    2    -3   21        0       0    ;    %U6
     85  25  0.00079  27.74  480   3    3    260   520    2    -3   20        0       0    ;    %U7
     55  10  0.00413  25.92  660   1    1    30    60     0    -1   11        0       0    ;    %U8
     55  10  0.00222  27.27  665   1    1    30    60     0    -1   11        0       0    ;    %U9
     55  10  0.00173  27.79  670   1    1    30    60     0    -1   11        0       0   ];   %U10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST CASE III%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%IEEE 10 UNIT-SYSTEM DATA CONSIDERING EMISSION COEFFICIENTS (ton/h)
 
%    %Pmax   %Pmin  %alpha    %beta     %gamma   %MUT %MDT %Shot %Scold %CSH %I.S %%D/URR  %%%MR %%%UN                         
% I=[ 455    150    0.00312  -0.24444   10.33908   8    8    4500  9000   5     8   0      0       0    ;    %U1
%     455    150    0.00312  -0.24444   10.33908   8    8    5000  10000  5     8   0      0       0    ;    %U2
%     130    20     0.00509  -0.40695   30.03910   5    5    550   1100   4    -5   0      0       0    ;    %U3
%     130    20     0.00509  -0.40695   30.03910   5    5    560   1120   4    -5   0      0       0    ;    %U4 
%     162    25     0.00344  -0.38132   32.00006   6    6    900   1800   4    -6   0      0       0    ;    %U5
%      80    20     0.00344  -0.38132   32.00006   3    3    170   340    2    -3   0      0       0    ;    %U6
%      85    25     0.00465  -0.39023   33.00056   3    3    260   520    2    -3   0      0       0    ;    %U7
%      55    10     0.00465  -0.39023   33.00056   1    1    30    60     0    -1   0      0       0    ;    %U8
%      55    10     0.00465  -0.39524   35.00056   1    1    30    60     0    -1   0      0       0    ;    %U9
%      55    10     0.00470  -0.39864   36.00012   1    1    30    60     0    -1   0      0       0   ];   %U10

% % % %%%%%%%%%%%%%%%%%%%DEMAND DATA FOR 10 UNIT SYSTEM%%%%%%%%%%%%%%%%%%%%

% % Demand & Spinning reserve requirement per hour (CASES I,III,IV,V)
P_D=[700;750;850;950;1000;1100;1150;1200;1300;1400;1450;1500;1400;1300;1200;1050;1000;1100;1200;1400;1300;1100;900;800];
% % 
% % % Demand & Spinning reserve requirement per hour (CASE VI-DRUC)
% % P_D=[700;750;850;950;1000;1100;1150;1200;1040;1120;1160;1200;1120;1040;1200;1050;1000;1100;1200;1120;1040;1100;900;800];
% 
% %%%Common for all cases
SRREQ=0.10*P_D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEST CASE II%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%IEEE 26 UNIT-SYSTEM DATA (ROBUSTNESS TEST)

% %Pmax %%Pmin   %%ai   %%bi     %%ci     %MUT  %MDT %Shot  %Scold  %CSH  %I.S %%D/URR  %%%MR %%%UN                         
% I=[ 400  100   0.0019  7.5031   311.9102   8     5    500    500     10    10   0         0       0    ;   %U1
%     400  100   0.0019  7.4921   310.0021   8     5    500    500     10    10   0         0       0    ;   %U2
%     350  140   0.0015  10.8616  177.0575   8     5    300    200     8     10   0         0       0    ;   %U3
%     197  68.95 0.0026  23.2000  260.1760   5     4    200    200     8     -4   0         0       0    ;   %U4
%     197  68.95 0.0026  23.1000  259.6490   5     4    200    200     8     -4   0         0       0    ;   %U5
%     197  68.95 0.0026  23.0000  259.1310   5     4    200    200     8     -4   0         0       0    ;   %U6
%     155  54.25 0.0049  10.7583  143.5972   5     3    150    150     6      5   0         0       0    ;   %U7
%     155  54.25 0.0048  10.7367  143.3719   5     3    150    150     6      5   0         0       0    ;   %U8
%     155  54.25 0.0047  10.7154  143.0288   5     3    150    150     6      5   0         0       0    ;   %U9
%     155  54.25 0.0046  10.6940  142.7348   5     3    150    150     6      5   0         0       0    ;   %U10
%     100  25.00 0.0060  18.2000  218.7752   4     2     70     70     4     -3   0         0       0    ;   %U11
%     100  25.00 0.0061  18.1000  218.3350   4     2     70     70     4     -3   0         0       0    ;   %U12
%     100  25.00 0.0062  18.0000  217.8952   4     2     70     70     4     -3   0         0       0    ;   %U13
%      76  15.20 0.0093  13.4073   81.6259   3     2     50     50     3      3   0         0       0    ;   %U14
%      76  15.20 0.0091  13.3805   81.4641   3     2     50     50     3      3   0         0       0    ;   %U15
%      76  15.20 0.0089  13.3538   81.2980   3     2     50     50     3      3   0         0       0    ;   %U16
%      76  15.20 0.0088  13.3272   81.1364   3     2     50     50     3      3   0         0       0    ;   %U17
%      20  4.000 0.0143  37.8896  118.8206   1     1     20     20     2     -1   0         0       0    ;   %U18
%      20  4.000 0.0136  37.7770  118.4576   1     1     20     20     2     -1   0         0       0    ;   %U19
%      20  4.000 0.0126  37.6637  118.1083   1     1     20     20     2     -1   0         0       0    ;   %U20
%      20  4.000 0.0120  37.5510  117.7551   1     1     20     20     2     -1   0         0       0    ;   %U21
%      12  2.400 0.0285  26.0611   24.8882   1     1      0      0     1     -1   0         0       0    ;   %U22
%      12  2.400 0.0284  25.9318   24.7605   1     1      0      0     1     -1   0         0       0    ;   %U23
%      12  2.400 0.0280  25.8027   24.6382   1     1      0      0     1     -1   0         0       0    ;   %U24
%      12  2.400 0.0265  25.6753   24.4110   1     1      0      0     1     -1   0         0       0    ;   %U25
%      12  2.400 0.0253  25.5472   24.3891   1     1      0      0     1     -1   0         0       0    ];  %U26
% % 
% % % %%%%%%%%%%%%%%%DEMAND DATA FOR 26 UNIT SYSTEM%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % %Demand & Spinning reserve requirement per hour 
% P_D=[2223;2052;1938;1881;1824;1825.5;1881;1995;2280;2508;2565;2593.5;2565;2508;2479.5;2479.5;2593.5;2850;2821.5;2764.5;2679;2662;2479.5;2308.5];
% 
% SRREQ=0.05*P_D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=size(I,1); %%%Number of TGU's
T=24;        %%%Study period duration

%%%%fmincon.m solver options for F_LIM_ED
OPTS = optimoptions('fmincon','Algorithm','sqp','display','off','ConstraintTolerance',1e-2,'OptimalityTolerance',1.0000e-02,'MaxIterations',100, 'StepTolerance',1e-2,'MaxFunctionEvaluations',100);

%%%%%%%%%%%%%%%%%%%%%%%%%%UNIT DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Generator limits
PGI_MAX=I(:,1); %%%Upper generation limit for each generator
PGI_MIN=I(:,2); %%%Lower generation limit for each generator

%%%I/O curve for each generator is modeled by a smooth quadratic function 
%%%i.e C(Pgi)=ai*Pgi^2+bi*Pgi+ci 
ai=I(:,3);
bi=I(:,4);
ci=I(:,5);

[I_C_SORT_EXS,IDX_EXS,I_C_SORT_INS,IDX_INS] = AFLC(ai,bi,ci,PGI_MAX); %%%Unit priority list assuming same fuel cost for every unit

%%MINIMUM TIME DATA
MUT=I(:,6);     %%%Minimum up-time for each generator
MDT=I(:,7);     %%%Minimum down-time for each generator

%%START-UP AND SHUT-DOWN COSTS
SU_H=I(:,8);
SU_C=I(:,9);
CSH=I(:,10);
init_status=I(:,11);

%%%%OTHER CONSTRAINTS
M_R=I(:,13);    %%MR units 
U_N=I(:,14);    %%Unavailable units

%%%% MR & UN hours & units 
    
    %%Hr begin  Hr end
MR=[  10          15  ];
UN=[                  ];

%%%% Prohibitive operative zones (POZ)
     %%%Unit    %%%POZ   
POZ=[  [1  150 165; 1  448,453]
       [2  90  110; 2  240 250]
       [8  20   30; 8  40  45]
       [10 12   17; 10 35  45]];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%MAIN LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Initiaize solution (NxT) schedule solution as an empty array
Init_SOL=zeros(N,T);

complete=ones(1,T); %%%Nth unit operation schedule is ON for all t
sum_Pgi_max=sum(PGI_MAX);

%%%%%Define swarm size (Either particle or bee)
SWARM_SZ=20;
PART_BEE=0;
YES_CNT=0;
NO_CNT=0;

population=zeros(SWARM_SZ,N*T);
total_COST=zeros(SWARM_SZ,1);

while PART_BEE<SWARM_SZ
    PART_BEE=PART_BEE+1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Repairing initial schedule according to fixed constraints%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Determine the Must-Run and Unavailable units (Rules B&C)

for j=1:N
    if M_R(j)==0 && U_N(j)==0
        continue
    elseif M_R(j)==1 
        Init_SOL(j,MR(1):MR(end))=1;
        continue
    elseif U_N(j)==1
        Init_SOL(j,UN(1):UN(end))=0;
        continue
    end
end
    
MR_idx=find(M_R>0); %%%Find the exception indices for must run units
UN_idx=find(U_N>0); %%%Find the exception indices for unavailable units 

%%%%%%%%%%%%Determine initial status of units %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

Tdiff_UP=zeros(N,3);%%%Exception indices for coupling constraints Tdiff_UP
Tdiff_DW=zeros(N,3);%%%Exception indices for coupling constraints Tdiff_DW

for k=1:N
   Tdiff_UP(k,1)=k;  
   Tdiff_DW(k,1)=k;  
     
   if ismember(k,MR_idx) || ismember(k,UN_idx) %%%If must run or unavailable unit, continue schedule repair without modifying INIT_SOL
       continue

   else
       P=1+(PGI_MAX(k)/sum_Pgi_max);
       
       if init_status(k)>0 %%Unit Online for ti=init_status hours
           Tdiff_U=init_status(k)-MUT(k);
           
           if Tdiff_U<0
                Tdiff_UP(k,2)=abs(Tdiff_U);
                
                Tdiff_UP(k,3)=1;
               
                Init_SOL(k,1:Tdiff_UP(k,2))=1;                
           else 
                Tdiff_UP(k,2)=0;
                
                Tdiff_UP(k,3)=0;
           end

           for j=Tdiff_UP(k,2)+1:T
               
               X=-1+(P-(-1))*rand;
               
               if X>=0
                   Init_SOL(k,j)=1;
               else
                   Init_SOL(k,j)=0;
               end    
           end

       else             %%Unit is offline for ti=init_status hours
           Tdiff_D=MDT(k)+init_status(k);
           
           if Tdiff_D>0
                Tdiff_DW(k,2)=abs(Tdiff_D);
               
                Init_SOL(k,1:Tdiff_D)=0;
                
                Tdiff_DW(k,3)=1;  
           else 
                Tdiff_DW(k,2)=0;
                
                Tdiff_DW(k,3)=0;
           end

           for j=Tdiff_DW(k,2)+1:T
                
               X=-1+(P-(-1))*rand;
               
               if X>=0
                   Init_SOL(k,j)=1;
               else
                   Init_SOL(k,j)=0;
               end    
           end
       end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Check & Repair Coupling constraints %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%Decompose commitment schedule into T single-hour unit combinations 

Exs_On_cap=zeros(3,T); %%Excess online capacity
Ins_On_cap=zeros(3,T); %%Insufficient online capacity

f_unit_ins_all=zeros(N,T); %%%This unit is not part of the exceptions MR,UN,Initialization

free_unit_ins=zeros(N,T);  %%%These units are offline and available to be turned off 

%%Check Excess online capacity first for full commitment
check_full_exs=ones(1,N)*PGI_MIN;     %%%Check if system has excess online capacity for a complete commitment
check_1=find(check_full_exs>P_D);

    if isempty(check_1)
        clear excp_exs f_unit_exs_all free_unit_exs Exs_On_cap
    else 
    %%%Excess online capacity schedule repair
        for t=1:T
            u_t=Init_SOL(:,t); %%Decompose into T single hour unit combinations

            Exs_On_cap(1,t)=u_t'*PGI_MIN;
            Exs_On_cap(2,t)=P_D(t);

            if Exs_On_cap(1,t)>P_D(t)
                Exs_On_cap(3,t)=1;
                format_exs='Excess capacity at hour %u. Repairing schedule...\n';
                fprintf(format_ins,t);  
            end
        end 
    end     

status_e=0;
status_f=0;
status_g=0;

status=[status_e,status_f,status_g];
trials=0;
max_trials=15;

    while ~all(status) 
    trials=trials+1;
    
    if trials>max_trials
        break
    end
    
    if trials==1
        %%%Insufficient online capacity initial schedule repair
        for t=1:T
            u_t=Init_SOL(:,t); %%Decompose into T single hour unit combinations

            Ins_On_cap(1,t)=u_t'*PGI_MAX;
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
               [bool, idx]=ismember(free_unit_ins_tmp,IDX_INS);

               idx=sort(idx);
               s_idx_off=size(idx,1);

               while Ins_On_cap(3,t)==1
                    for i=1:s_idx_off
                         u_t(IDX_INS(idx(i)))=1;

                         Ins_On_cap(1,t)=u_t'*PGI_MAX;

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

                   commit_first=IDX_INS(k); %%%Commit the cheapeast generator 
                   u_t=Init_SOL(commit_first,:);

                   if u_t==complete
                      continue
                   else
                      if u_t(v_t)==0
                         u_t(v_t)=1;

                         Init_SOL(commit_first,:)=u_t;

                         check_E=check_SR_PD(T,P_D,SRREQ,Init_SOL,PGI_MAX);  
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

       [check_E,status_e]=check_SR_PD(T,P_D,SRREQ,Init_SOL,PGI_MAX);

    else %%%Pivot heuristic for MUT/MDT correction       
        
        if prev_SOL==Init_SOL
            
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
                   x=1;
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
                no_pivot=0;
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

           [check_E,status_e]=check_SR_PD(T,P_D,SRREQ,Init_SOL,PGI_MAX);
        end   
    end
       status=[status_e,status_f,status_g];
    end
    
    if ~all(status==1)
        PART_BEE=PART_BEE-1;
        continue
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

               Rxs=uj_t'*PGI_MAX;

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

               [intr_on,intr_off] = count_intervals(n_sch);

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
    
    Init_SOL=temp_SOL;
    
    [check_E,status_e]=check_SR_PD(T,P_D,SRREQ,Init_SOL,PGI_MAX);
    
    status=[status_e,status_f,status_g];
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%Economic-Dispatch%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%Ramping Rate constraints not included

[P_SOL,P_srv,P_COST,tot_gen_COST,itt] = F_LIM_ED(Init_SOL,T,I,P_D,OPTS);

%%%Ramping Rate constraints included

% [P_SOL,P_srv,P_COST,tot_gen_COST,itt,Init_SOL] = F_LIM_ED_RR(Init_SOL,T,I,P_D,IDX_INS,OPTS,POZ);

dff=abs(P_srv-P_D');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Start-up-cost computation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
[su_COST,tot_su_COST] = SU_COST(N,Init_SOL,init_status,SU_H,SU_C,CSH,MUT,MDT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%Display results%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if all((dff)<1e-2)
        YES_CNT=YES_CNT+1;

    else
        NO_CNT=NO_CNT+1;
        disp(PART_BEE);
    end
        
    population(PART_BEE,:)=reshape(Init_SOL,1,N*T);        

    total_COST(PART_BEE)=tot_gen_COST+tot_su_COST;    
    
end

disp('%%%%%%%%%%%%%%%%%%%%%SWARM GENERATOR MODULE%%%%%%%%%%%%%%%%%%%%%%%%');
disp('%%%%%%%%%%%%%%%%Coded by Vick @2018 NSYSU EMS-LAB%%%%%%%%%%%%%%%%%%');
fprintf('\n');

format_swarm='%u initial particles created for swarm.\n';
fprintf(format_swarm,SWARM_SZ);
fprintf('\n');

format_particle='%u particles succesfully repaired satisfying all constraints.\n';
fprintf(format_particle,YES_CNT);
fprintf('\n');

disp('Total cost for each particle (Generation+Startup):')
fprintf('\n');
disp(total_COST);
fprintf('\n');

toc