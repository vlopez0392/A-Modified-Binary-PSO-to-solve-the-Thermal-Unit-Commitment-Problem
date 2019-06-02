function [Init_SOL,u] = Recomm_swp(I,T,t,Init_SOL,PD,Pgi_max,IDX_INS,type)
%Re-Commit or swap units to satisfy  System Power Demand at time t
    
    %%%%Input data
    MUT=I(:,6); 
    MDT=I(:,7);
    init_status=I(:,11);
    
    switch type
        case 'under'
        u=Init_SOL(:,t);
        
        %%%Determine the zero positions that satisfy the MUT/MDT constraints
        z_p=find(u<1);
        sz_zp=numel(z_p);
        
        temp_SOL=Init_SOL;
        temp_SOL(u<1,t)=1;
        
        [avlb]=check_MUT_MDT(sz_zp,T,temp_SOL(z_p,:),init_status(z_p),MUT(z_p),MDT(z_p));
        
        for j=1:sz_zp 
           if avlb(j,1)==1 || avlb(j,2)==1
               z_p(j)=0;
           else
               continue
           end    
        end
        
        z_p=intersect(IDX_INS,z_p(z_p>0));
        sz_zp=numel(z_p);
        
        %%%%Re-commit/swap units to satisfy demand 
        iter=0;
        
        while sum(Pgi_max(u>0))<PD
            iter=iter+1;
            
            if iter> sz_zp
                break
            end
            u(z_p(iter))=1;
        end
        Init_SOL(:,t)=u;
    end
end














