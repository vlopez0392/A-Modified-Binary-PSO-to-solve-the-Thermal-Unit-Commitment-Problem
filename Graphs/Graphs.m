%%%%%%%%%%%%Graphing best results by VS-BPSO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%CASE I
%%%%Convergence behavior per velocity mapping function
load('S1.mat')
load('S2.mat')
load('S3.mat')
load('S4.mat')
load('V1.mat')
load('V2.mat')
load('V3.mat')
load('V4.mat')
load('LOAD_DEMAND.mat')
load('P_SOL_stacked.mat');

%%%%Plotting the results 
%%%%Convergence results 
figure(1)
plot(S1(:),'k','LineStyle','--','LineWidth',1);
hold on 
plot(S2(:),'k','LineStyle','-','LineWidth',1);
plot(S3(:),'k','LineStyle',':','LineWidth',1);
plot(S4(:),'k','LineStyle','-.','LineWidth',1);
plot(V1(:),'b','LineStyle','--','LineWidth',1);
plot(V2(:),'b','LineStyle','-','LineWidth',1);
plot(V3(:),'b','LineStyle',':','LineWidth',1);
plot(V4(:),'b','LineStyle','-.','LineWidth',1);
hold off

ylim([563000,570000])
xlabel('Iteration');
ylabel('Total cost in $');

legend('S-BPSO1','S-BPSO2', 'S-BPSO3','S-BPSO4','V-BPSO1','V-BPSO2','V-BPSO3','V-BPSO4')

%%%%Output power results 
figure(2)
bar(P_SOL_stacked','stacked')
colormap parula
hold on
plot(P_D(:),'k','LineStyle','-','LineWidth',1.5);
hold off
set(gca,'XTick',1:24);
ylim([0,1800])
xlabel('Hour');
ylabel('Output power (MW)');
legend('Unit 1','Unit 2', 'Unit 3','Unit 4','Unit 5','Unit 6','Unit 7','Unit 8','Unit 9','Unit 10', 'Load demand','Location','eastoutside','Orientation','Vertical')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%CASE II
%%%%Convergence behavior for VBPSO3 
load('V3_RR.mat')
load('P_SOL_OPT_RR');

U5_BC=P_SOL_stacked(5,:);
U5_RR=P_SOL_OPT_RR(5,:);
B_PLT=[U5_BC;U5_RR];

sz_BLT=numel(B_PLT);

for i=1:sz_BLT
   if B_PLT(i)==0
       B_PLT(i)=1;
   end
end

%%%%Convergence results 
figure(3)
plot(V3_RR(:),'b','LineStyle','-.','LineWidth',3);
ylim([564500,570000])
xlabel('Iteration');
ylabel('Total cost in $');

%%%%Output power results 
figure(4)
c = categorical({'Case I','Case II'});
bar(c,B_PLT)
ylabel('Output power (MW)');

xlabel('Cases');

legend('H1','H2', 'H3','H4','H5','H6','H7','H8','H9','H10', 'H11','H12','H13', 'H14','H15','H16','H17','H18','H19','H20','H21', 'H22','H23', 'H24','Location','eastoutside','Orientation','Vertical')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%CASE III
load('V3_EM.mat')
load('P_SOL_OPT_EM');
load('EM.mat')

%%%%Convergence behavior for VBPSO3 
figure(5) 
plot(V3_EM(:),'b','LineStyle','--','LineWidth',2);
xlabel('Iteration');
ylabel('Total Emission in TON');
legend('Emissions for case III')
title('Case III');

%%%%Output power & emission results 
figure(6) 
plot(EM(:),'r','LineStyle','-.','LineWidth',1);
xlabel('Iteration');
ylabel('Total Emission in TON');
title('Case I');
legend('Emissions for case I','Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%CASE IV
%%%%Convergence behavior for VBPSO3 
load('V3_DR.mat');
load('PDR_OPT.mat');
load('P_D_DR.mat');

figure(7)
plot(V3_DR(:),'m','LineStyle','-.','LineWidth',3);
ylim([502000,520000])
xlabel('Iteration');
ylabel('Total cost in $');

%%%%Load Demand
figure(8)
% plot(P_D_DR(:),'k','LineStyle','-','LineWidth',2);
hold on
plot(P_D(:),'r','LineStyle','-','LineWidth',1.5);
SR=P_D(:)+0.10*P_D(:);
plot(SR,'m','LineStyle','-.','LineWidth',1.5);
hold off
set(gca,'XTick',1:24);
ylim([600,1750])
xlabel('Hour');
ylabel('Load Demand (MW)');
% legend('DR-Load demand','Original Load Demand','Location','south','Orientation','Horizontal')
legend('Load Demand','Spinning Reserve requirement','Location','southoutside','Orientation','Horizontal')

%%%%Output power results 
figure(9)
bar(PDR_OPT','stacked')
colormap parula
hold on
plot(P_D_DR(:),'k','LineStyle','-','LineWidth',2);
plot(P_D(:),'r','LineStyle','-','LineWidth',1.5);
hold off
set(gca,'XTick',1:24);
ylim([0,1800])
xlabel('Hour');
ylabel('Output power (MW)');
legend('Unit 1','Unit 2', 'Unit 3','Unit 4','Unit 5','Unit 6','Unit 7','Unit 8','Unit 9','Unit 10', 'DR-Load demand','Original Load Demand','Location','eastoutside','Orientation','Vertical')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%CASE V
load('C_26_OPTTTT.mat');
load('P_SOL_OPT_26_OPTTTTT.mat');
load('P_D_26.mat')

%%%%Convergence behavior for VBPSO3 
figure(10)
plot(C26_OPTTTT(:),'b','LineStyle','-.','LineWidth',3);
ylim([760000,790000])
xlabel('Iteration');
ylabel('Total cost in $');

%%%Output power results 
figure(11)
bar(P_SOL_26_OPTTTTT,'stacked')
colormap parula
hold on
plot(P_D_26(:),'k','LineStyle','-','LineWidth',1.5);
hold off
set(gca,'XTick',1:24);
ylim([0,2900])
xlabel('Hour');
ylabel('Output power (MW)');

legend('U1','U2','U3','U4','U5','U6','U7','U8','U9','U10', 'U11','U12','U13', 'U14','U15','U16','U17','U18','U19','U20','U21','U22','U23','U24','U25','U26','Load Demand','Location','eastoutside','Orientation','Vertical')


%%%%Load Demand
figure(12)
hold on
plot(P_D_26(:),'r','LineStyle','-','LineWidth',1.5);
SR=P_D_26(:)+0.05*P_D_26(:);
plot(SR,'m','LineStyle','-.','LineWidth',1.5);
hold off
set(gca,'XTick',1:24);
ylim([1700,3100])
xlabel('Hour');
ylabel('Load Demand (MW)');
% legend('DR-Load demand','Original Load Demand','Location','south','Orientation','Horizontal')
legend('Load Demand','Spinning Reserve requirement','Location','southoutside','Orientation','Horizontal')