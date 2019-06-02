function [su_COST,tot_su_COST,H_or_C] = SU_COST(N,Init_SOL,init_status,SU_H,SU_C,CSH,MUT,MDT)
%%Computes the start up and shut down costs (to expand) of the UC schedule 
    su_COST=zeros(N,2);
    
    H_or_C=zeros(N,2);
    
    for k=1:N
        n_sch=Init_SOL(k,:);

        %%Compute upper and lower limits to determine HSC or CSC
        ul=MDT(k)+CSH(k);
        ll=MUT(k);

        [intr_on,intr_off] = count_intervals(n_sch);

        if intr_off==0    %%%All units on during schedule
            if init_status(k)<0
                su_COST(k,1)=1;
                if abs(init_status(k))<=ul && abs(init_status(k))>=ll

                    su_COST(k,2)=SU_H(k);
                else

                    su_COST(k,2)=SU_C(k);
                end
            end

        elseif intr_on==0 %%%All units off during schedule (Expand if shut-down costs need to be included);


        else 
            no_intr_off=size(intr_off,1);

            for j=1:no_intr_off

                idx_off_e=intr_off(j,2);
                hours_off=intr_off(j,3);

                if j==1 
                    if idx_off_e==24
                         if init_status(k)<1
                            Ti_off=abs(init_status(k));
                            
                            if Ti_off>ul
                                su_COST(k,2)=SU_C(k);
                            else
                                su_COST(k,2)=SU_H(k);
                            end
                            
                            su_COST(k,1)=su_COST(k,1)+1;
                         end
                    else
                        if init_status(k)>1
                            Ti_off=intr_off(j,2);
                        else
                            Ti_off=abs(init_status(k))+hours_off;
                        end

                        if Ti_off>ul
                            su_COST(k,2)=SU_C(k);
                        else
                            su_COST(k,2)=SU_H(k);
                        end
                        
                        su_COST(k,1)=su_COST(k,1)+1;
                    end
                else

                    if idx_off_e==24
                        continue
                    else
                        Ti_off=hours_off;

                        if Ti_off>ul
                            su_COST(k,2)=su_COST(k,2)+SU_C(k);
                        else
                            su_COST(k,2)=su_COST(k,2)+SU_H(k);
                        end

                        su_COST(k,1)=su_COST(k,1)+1;
                    end
                end
            end
        end
    end
    
    tot_su_COST=sum(su_COST(:,2));
end
