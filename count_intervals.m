function [intr_on,intr_off] = count_intervals(n_sch)
% Returns ON/OFF intervals & their respective duration in n_sch

    complete_ones=ones(1,numel(n_sch));
    complete_zeros=zeros(1,numel(n_sch));

    if n_sch==complete_ones
        intr_on(1:3)=[1,24,24];

        intr_off=0;

    elseif n_sch==complete_zeros
        intr_off(1:3)=[1,24,24];
        
        intr_on=0;
            
    else
        %%%Determine number of on/off intervals to begin counting process
        no_intr_on=numel(find(diff(find(n_sch>0))>1))+1;
        intr_on=zeros(no_intr_on,3);

        no_intr_off=numel(find(diff(find(n_sch<1))>1))+1;
        intr_off=zeros( no_intr_off,3);

        %%%Detetmine the duration of each of the intervals (Store idx_begin,
        %%%idx_end)

        %%%ON intervals
        look_on=find(n_sch>0);
        sz_look_on=numel(look_on)-1;
        int_UP=diff(look_on);

        if all(int_UP>1) %%%Individual ON One hour intervals 

            for j=1:no_intr_on
                intr_on(j,1:2)=look_on(j);
                intr_on(j,3)=1;
            end

        else %%%Several ON intervals of different duration
            on=1;
            bool=true;

            for k=1:sz_look_on

                 if bool==true
                     if on==1
                         intr_on(on,1)=look_on(1);
                     else
                         intr_on(on,1)=e_t;
                     end

                     bool=false;
                 end  

                 e_t_1=look_on(k);
                 e_t=look_on(k+1);

                 if e_t-e_t_1==1
                     continue
                 else
                     intr_on(on,2)=e_t_1;
                     intr_on(on,3)=e_t_1-intr_on(on,1)+1;

                     bool=true;
                     on=on+1;
                     continue
                 end
            end

            if intr_on(on,1)==0
                intr_on(on,2)=look_on(end);
                intr_on(on,1)=intr_on(on,2);
                intr_on(on,3)=1;
            else
                intr_on(on,2)=look_on(end);
                intr_on(on,3)=intr_on(on,2)-intr_on(on,1)+1; 
            end 
        end

        %%%OFF intervals
        look_off=find(n_sch<1);
        sz_look_off=numel(look_off)-1;
        int_DW=diff(look_off);

        if all(int_DW>1) %%%Individual one hour OFF intervals 

            for j=1:no_intr_off
                intr_off(j,1:2)=look_off(j);
                intr_off(j,3)=1;
            end

        else  %%%Several OFF intervals of different duration
            off=1;
            bool=true;

            for k=1:sz_look_off

                 if bool==true
                     if off==1
                         intr_off(off,1)=look_off(1);
                     else
                         intr_off(off,1)=e_t;
                     end

                     bool=false;
                 end  

                 e_t_1=look_off(k);
                 e_t=look_off(k+1);

                 if e_t-e_t_1==1
                     continue
                 else
                     intr_off(off,2)=e_t_1;
                     intr_off(off,3)=e_t_1-intr_off(off,1)+1;

                     bool=true;
                     off=off+1;
                     continue
                 end
            end

            if intr_off(off,1)==0
                intr_off(off,2)=look_off(end);
                intr_off(off,1)=intr_off(off,2);
                intr_off(off,3)=1;
            else
                intr_off(off,2)=look_off(end);
                intr_off(off,3)=intr_off(off,2)-intr_off(off,1)+1; 
            end
        end
    end
end

