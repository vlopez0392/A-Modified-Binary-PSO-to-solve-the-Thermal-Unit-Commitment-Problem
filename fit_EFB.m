function fit= fit_EFB(iter,EFB,a_P,b_P,c_P,d_P,e_P,P_gen,varargin )
%%%Computes the fitness for the experienced foragers (Either new or old)
%%%Computes the fitness of the current Gbest

min_args=9;
max_args=11;
narginchk(min_args,max_args)

s_EFB=size(EFB,1);
%%%%Checks for correct combination of inputs 
    if length(varargin)==2
        P_SOL=varargin{1};
        idx_or=varargin{2};
        
        if iter==2
            N=size(P_SOL,2);
            D=size(P_SOL,1);
            
            if N< size(idx_or,2);%%%N of units=N; N>1 for all iter
                [P_SOL,idx_or]=swap(P_SOL,idx_or);
                N=size(P_SOL,2);  %#ok<NASGU>

            elseif D~=size(idx_or,1)
                e=1; %%%Error 1: Input of incorrect combination of arguments
                
                if e==1
                disp('Incorrect combination of arguments');
                end
            end
        end
        
        %%%Computes the fitness values of the new experienced forager bees in each iter
        %(fit_x)
                fit=zeros(s_EFB,1);
                for i=1:s_EFB
                    fit(i)=1./(1+a_P'*(P_SOL(idx_or(i),:)).^2'+b_P'*P_SOL(idx_or(i),:)'+c_P'+abs(sum(d_P'*sin(e_P'*(P_gen(:,4)'-P_SOL(i,:))'))));
                end
     else   
            if size(varargin{1},1)<s_EFB %%%N>1 for all iter
                c='g';
                Gbest=varargin{1};
            else 
                c='b';
                best_food_EFB=varargin{1}; 
            end
       
        switch c
            case 'b'
            %%%Computes the f.values of the bees for the previously found best food sources
            %(fit_b)

                fit=zeros(s_EFB,1);

                for j=1:s_EFB
                    fit(j)=1./(1+a_P'*(best_food_EFB(j,:)).^2'+b_P'*best_food_EFB(j,:)'+c_P'+abs(sum(d_P'*sin(e_P'*(P_gen(:,4)'-best_food_EFB(j,:))'))));
                end

            case 'g' 
            %%%Computes the f value of the current Gbest food source %fit_g 

                fit=1./(1+a_P'*(Gbest).^2'+b_P'*(Gbest)'+c_P'+abs(sum(d_P'*sin(e_P'*(P_gen(:,4)'-Gbest)'))));
        end
    end
end


