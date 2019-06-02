function zero_idx = until_zero(n_sch,varargin)
%%%Returns the first or last occurring zero closest to the pivot 
%%%This function is part of the modified repair strategy by pivot heuristic

min_args=3;
max_args=3;
narginchk(min_args,max_args)

dec=varargin{2};
zeros=find(n_sch<1);

    switch dec
        case 1 %%%%Look to the right
            idx_end=varargin{1};

            zeros_right=zeros(zeros>idx_end);
            
            diff_z_idx=diff(zeros_right);

            if all(diff_z_idx==1)
                zero_idx=max(zeros_right);
            else
                sz_l=numel(zeros_right);
                
                for t=1:sz_l
                    if t<sz_l
                        e_t_1=zeros_right(t);
                        e_t=zeros_right(t+1);
                        
                        if e_t-e_t_1==1
                            continue
                        else
                           zero_idx=e_t_1;
                           break
                        end
                    else
                        continue
                    end
                end
            end
        case 0 %%%Look to the left
            idx_begin=varargin{1};

            zeros_left=zeros(zeros<idx_begin);
            
            diff_z_idx=diff(zeros_left);
            
            if all(diff_z_idx==1)
                zero_idx=min(zeros_left);
            else
                sz_l=numel(zeros_left);
                
                for t=sz_l:-1:1
                    if t>1
                        e_t=zeros_left(t);
                        e_t_1=zeros_left(t-1);
                        
                        if e_t-e_t_1==1
                            continue
                        else
                           zero_idx=e_t;
                           break
                        end
                    else
                        continue
                    end
                end        
            end
    end     
end

