cvx_begin 
    variables p1 p2 p3 p4;
        %%%Generator parameters             
        
        %%%SAMPLE TEST I unit characteristics                   %V.P
        %            %ai    bi      ci      %PGi min  %Pgi max  %di    %ei
%         P_gen= [0.001562    7.92    561     150        600      0      0     ;
%                 0.001940    7.85    310     100        400      0      0     ;
%                 0.004820    7.97    78      50         200      0      0     ];

        %%%SAMPLE TEST II unit characteristics           %V.P
        %           %ai   bi      ci  %PGi min  %Pgi max  %di    %ei
%         P_gen=[0.0063   5.32    500   100       600       0      0   ; 
%                0.0413   5.13    400   150       700       0      0   ;
%                0.0093   4.83    450   50        400       0      0   ;  
%                0.0073   4.61    550   50        200       0      0  ]; 


        %%%SAMPLE TEST III unit characteristics                   %V.P
                   %ai    bi      ci      %PGi min  %Pgi max  %di    %ei
%         P_gen= [0.001280    6.48    459     150        600      0      0     ;
%                 0.001940    7.85    310     100        400      0      0     ;
%                 0.004820    7.97    78      50         200      0      0     ];    

P_gen=  [0.00048  16.19  1000	 150	    455    0    0;                
         0.00031  17.26  970	 150	    455    0    0]; 	
     
        a_P=P_gen(:,1);
        b_P=P_gen(:,2);
        c_P=sum(P_gen(:,3));
        
        %%%%Solution vector
        PSOL=[p1;p2];
        
        %%%I/O of generator modeled by quadratic function
        F=a_P'*PSOL.^2+b_P'*PSOL+c_P';
   
        minimize F
            subject to %%Constraints 
                PSOL>=P_gen(:,4); %#ok<VUNUS>
                PSOL<=P_gen(:,5); %#ok<VUNUS>
                sum(PSOL)==700; %#ok<EQEFF>
cvx_end
