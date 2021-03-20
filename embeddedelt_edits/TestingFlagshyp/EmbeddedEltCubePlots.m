clear; clc; close all;

% fid0=fopen("EmbeddedElt_Mod_Load0.txt",'r');
% fid1=fopen("EmbeddedElt_Mod_Load102.txt",'r');
% fid4=fopen("EmbeddedElt_Mod_Load4.txt",'r');



fid0=fopen("Flagshyp_1h_1e_en.txt",'r');
fid1=fopen("Flagshyp_1h_1e-HalfStiff.txt",'r');
fid4=fopen("Flagshyp_1h_1e-HalfDense.txt",'r');


Legend3=["Original","HalfStiff","HalfDens"];
Color0 = '#00008B';
Color1 = '#006400';
Color4 = '#B22222';
set(0,'defaultfigurecolor',[1 1 1]);

for i=1:4
    tline = fgetl(fid0);
    tline = fgetl(fid1);
    tline = fgetl(fid4);
end
Energy0 = fscanf(fid0, '%f ',[10,201]);
Energy1 = fscanf(fid1, '%f ',[10,201]);
Energy4 = fscanf(fid4, '%f ',[10,201]);

for i=1:4
    tline = fgetl(fid0);
    tline = fgetl(fid1);
    tline = fgetl(fid4);
end

Elt0 = fscanf(fid0, '%f ',[7,101]);
Elt1 = fscanf(fid1, '%f ',[7,101]);
Elt4 = fscanf(fid4, '%f ',[7,101]);

for i=1:3
    tline = fgetl(fid0);
    tline = fgetl(fid1);
    tline = fgetl(fid4);
end

Nodes0 = fscanf(fid0, '%f ',[17,101]);
Nodes1 = fscanf(fid1, '%f ',[17,101]);
Nodes4 = fscanf(fid4, '%f ',[17,101]);

fclose('all');
%% Energy Plots

EnergyLeg=["Time","Artificial Strain Energy", "Creep Dissapation","Total Strain Energy","Kinetic Energy","Plastic Dissipation", "Elastic Strain","Viscous Disipation", "External Work","Total Energy"];

%Elt0
figure(); hold on; grid on;
Legend0=[];
plot(Energy0(1,:),Energy0(2,:),'b','Linewidth', 2);
Legend0=[Legend0 EnergyLeg(2)];
% plot(Energy0(1,:),Energy0(3,:),'Linewidth', 1);
% Legend0=[Legend0 EnergyLeg(3)];
plot(Energy0(1,:),Energy0(4,:),'r','Linewidth', 2);
Legend0=[Legend0 EnergyLeg(4)];
plot(Energy0(1,:),Energy0(5,:),'k','Linewidth', 2);
Legend0=[Legend0 EnergyLeg(5)];
% plot(Energy0(1,:),Energy0(6,:),'Linewidth', 1);
% Legend0=[Legend0 EnergyLeg(6)];
plot(Energy0(1,:),Energy0(7,:),'m','Linewidth', 2);
Legend0=[Legend0 EnergyLeg(7)];
plot(Energy0(1,:),-Energy0(8,:),'g','Linewidth', 2);
Legend0=[Legend0 EnergyLeg(8)];
plot(Energy0(1,:),Energy0(9,:),'c','Linewidth', 2);
Legend0=[Legend0 EnergyLeg(9)];
xlabel("Time (s)");
ylabel("Energy (J)");
title("No Fibers");
legend(Legend0);
hold off;

%Elt1
figure(); hold on; grid on;
Legend1=[];
plot(Energy1(1,:),Energy1(2,:),'b','Linewidth', 2);
Legend1=[Legend1 EnergyLeg(2)];
% plot(Energy1(1,:),Energy1(3,:),'Linewidth', 1);
% Legend1=[Legend1 EnergyLeg(3)];
plot(Energy1(1,:),Energy1(4,:),'r','Linewidth', 2);
Legend1=[Legend1 EnergyLeg(4)];
plot(Energy1(1,:),Energy1(5,:),'k','Linewidth', 2);
Legend1=[Legend1 EnergyLeg(5)];
% plot(Energy1(1,:),Energy1(6,:),'Linewidth', 1);
% Legend1=[Legend1 EnergyLeg(6)];
plot(Energy1(1,:),Energy1(7,:),'m','Linewidth', 2);
Legend1=[Legend1 EnergyLeg(7)];
plot(Energy1(1,:),-Energy1(8,:),'g','Linewidth', 2);
Legend1=[Legend1 EnergyLeg(8)];
plot(Energy1(1,:),Energy1(9,:),'c','Linewidth', 2);
Legend1=[Legend1 EnergyLeg(9)];
xlabel("Time (s)");
ylabel("Energy (J)");
title("102 Fibers");
legend(Legend1);
hold off;

%Elt4
figure(); hold on; grid on;
Legend4=[];
plot(Energy4(1,:),Energy4(2,:),'b','Linewidth', 2);
Legend4=[Legend4 EnergyLeg(2)];
% plot(Energy4(1,:),Energy4(3,:),'Linewidth', 1);
% Legend4=[Legend4 EnergyLeg(3)];
plot(Energy4(1,:),Energy4(4,:),'r','Linewidth', 2);
Legend4=[Legend4 EnergyLeg(4)];
plot(Energy4(1,:),Energy4(5,:),'k','Linewidth', 2);
Legend4=[Legend4 EnergyLeg(5)];
% plot(Energy4(1,:),Energy4(6,:),'Linewidth', 1);
% Legend4=[Legend4 EnergyLeg(6)];
plot(Energy4(1,:),Energy4(7,:),'m','Linewidth', 2);
Legend4=[Legend4 EnergyLeg(7)];
plot(Energy4(1,:),-Energy4(8,:),'g','Linewidth', 2);
Legend4=[Legend4 EnergyLeg(8)];
plot(Energy4(1,:),Energy4(9,:),'c','Linewidth', 2);
Legend4=[Legend4 EnergyLeg(9)];
xlabel("Time (s)");
ylabel("Energy (J)");
title("4 Fibers");
legend(Legend4);
hold off;

for i=2:10
    figure; hold on; grid on;
    plot(Energy0(1,:),Energy0(i,:),'Color', Color0,'Linewidth', 2);
    plot(Energy1(1,:),Energy1(i,:),'Color', Color1,'Linewidth', 2);
    plot(Energy4(1,:),Energy4(i,:),'Color', Color4,'Linewidth', 2);
    xlabel("Time (s)");
    ylabel("Energy (J)");
    title(EnergyLeg(i));
    legend(Legend3);
%     legend("0Fiber","102Fiber","4Fiber");
end
%% Energy Summations 

%["Time","Artificial Strain Energy", "Creep Dissapation","Total Strain Energy","Kinetic Energy","Plastic Dissipation", "Elastic Strain","Viscous Disipation", "External Work","Total Energy"];

    %Energy Balance All Sims
    figure; hold on; grid on;
    plot(Energy0(1,:),Energy0(4,:),'Color', '#00008B' ,'Linewidth', 1,'MarkerIndices',1:5:100); 
    plot(Energy0(1,:),Energy0(4,:)-Energy0(9,:)+Energy0(5,:)+Energy0(2,:),'Color','#4169E1','Linewidth', 1,'MarkerIndices',1:5:100);
    plot(Energy0(1,:),-Energy0(9,:),'Color', '#1E90FF','Linewidth', 1,'MarkerIndices',1:5:100);
    
    plot(Energy1(1,:),Energy1(4,:),'Color','#006400','Linewidth', 1,'MarkerIndices',1:5:100);
    plot(Energy1(1,:),Energy1(4,:)-Energy1(9,:)+Energy1(5,:)+Energy1(2,:),'Color',  '#228B22','Linewidth', 1,'MarkerIndices',1:5:100);
    plot(Energy1(1,:),-Energy1(9,:),'Color','#3CB371' ,'Linewidth', 1,'MarkerIndices',1:5:100);
    
    plot(Energy4(1,:),Energy4(4,:),'Color','#B22222' ,'Linewidth', 1,'MarkerIndices',1:5:100);
    plot(Energy4(1,:),Energy4(4,:)-Energy4(9,:)+Energy4(5,:)+Energy4(2,:),'Color','#FF4500' ,'Linewidth', 1,'MarkerIndices',1:5:100);
    plot(Energy4(1,:),-Energy4(9,:), 'Color','#FF6347','Linewidth', 1,'MarkerIndices',1:5:100);
    
    
    xlabel("Time (s)");
    ylabel("Energy (J)");
    title("Energy Balance");   
    legend("0Fiber Strain Energy","0Fiber Total Energy","0Fiber External Work","102Fiber Strain Energy","102Fiber Total Energy","102Fiber External Work","4Fiber Strain Energy","4Fiber Total Energy","4Fiber External Work");
    
    %Summation of Energy Terms
%     figure; hold on; grid on;
%     plot(Energy0(1,:),Energy0(4,:)-Energy0(8,:)+Energy0(5,:)+Energy0(2,:),'Linewidth', 1);
%     plot(Energy1(1,:),Energy1(4,:)-Energy1(7,:)+Energy1(5,:)+Energy1(2,:),'Linewidth', 1);
%     plot(Energy4(1,:),Energy4(4,:)-Energy4(7,:)+Energy4(5,:)+Energy4(2,:),'Linewidth', 1);
%     xlabel("Time (s)");
%     ylabel("Energy (J)");
%     title("Summation of All Energy");
%     legend("Elt0","Elt1","Elt4");
%     
     
    %Conservation of Energy 0
    figure; hold on; grid on;
    plot(Energy0(1,:),Energy0(4,:),'Color', '#00008B','Linewidth', 1);
    plot(Energy0(1,:),-Energy0(9,:),'Color','#4169E1','Linewidth', 1);
    plot(Energy0(1,:),Energy0(4,:)-Energy0(9,:),'Color', '#1E90FF','Linewidth', 1);
    xlabel("Time (s)");
    ylabel("Energy (J)");
    title("Elt0");
    legend("Strain Energy","External Work","Sum");
    
    %Conservation of Energy 1
    figure; hold on; grid on;
    plot(Energy1(1,:),Energy1(4,:),'Color','#006400','Linewidth', 1);
    plot(Energy1(1,:),-Energy1(9,:),'Color',  '#228B22','Linewidth', 1);
    plot(Energy1(1,:),Energy1(4,:)-Energy1(9,:),'Color','#3CB371','Linewidth', 1);
    xlabel("Time (s)");
    ylabel("Energy (J)");
    title("Elt102");
    legend("Strain Energy","External Work","Sum");
    
    %Conservaion of Energy 4
    figure; hold on; grid on;
    plot(Energy4(1,:),Energy4(4,:),'Color','#B22222' ,'Linewidth', 1);
    plot(Energy4(1,:),-Energy4(9,:),'Color','#FF4500','Linewidth', 1);
    plot(Energy4(1,:),Energy4(4,:)-Energy4(9,:),'Color','#FF6347','Linewidth', 1);
    xlabel("Time (s)");
    ylabel("Energy (J)");
    title("Elt4");
    legend("Strain Energy","External Work","Sum");

%% Element Outputs

EltLeg=["Time (s)","ER33","ER Max Principal","LE33","LE Max Principal","S Max Principal","S33"];

%ER33
% figure; hold on; grid on;
% plot(Elt0(1,:),Elt0(2,:),'Linewidth', 2);
% plot(Elt1(1,:),Elt1(2,:),'Linewidth', 2);
% plot(Elt4(1,:),Elt4(2,:),'Linewidth', 2);
% xlabel("Time (s)");
% ylabel("Strain Rate (1/s)");
% title(EltLeg(2));
% legend(Legend3);

%ER Max Principle
figure; hold on; grid on;
plot(Elt0(1,:),Elt0(3,:),'Linewidth', 2);
plot(Elt1(1,:),Elt1(3,:),'Linewidth', 2);
plot(Elt4(1,:),Elt4(3,:),'Linewidth', 2);
xlabel("Time (s)");
ylabel("Strain Rate (1/s)");
title(EltLeg(3));
legend(Legend3);

%LE33
% figure; hold on; grid on;
% plot(Elt0(1,:),Elt0(4,:),'Linewidth', 2);
% plot(Elt1(1,:),Elt1(4,:),'Linewidth', 2);
% plot(Elt4(1,:),Elt4(4,:),'Linewidth', 2);
% xlabel("Time (s)");
% ylabel("Strain");
% title(EltLeg(4));
% legend(Legend3);
% LE33=[Elt0(4,:);Elt1(4,:);Elt4(4,:)];

%LE Max Principle
figure; hold on; grid on;
plot(Elt0(1,:),Elt0(5,:),'Color', Color0,'Linewidth', 2);
plot(Elt1(1,:),Elt1(5,:),'Color', Color1,'Linewidth', 2);
plot(Elt4(1,:),Elt4(5,:),'Color', Color4,'Linewidth', 2);
xlabel("Time (s)");
ylabel("Strain");
title(EltLeg(5));
legend(Legend3);

%S Max Principle
% figure; hold on; grid on;
% plot(Elt0(1,:),Elt0(6,:),'Linewidth', 2);
% plot(Elt1(1,:),Elt1(6,:),'Linewidth', 2);
% plot(Elt4(1,:),Elt4(6,:),'Linewidth', 2);
% xlabel("Time (s)");
% ylabel("Stress (Pa)");
% title(EltLeg(6));
% legend(Legend3);

%S33
% figure; hold on; grid on;
% plot(Elt0(1,:),Elt0(7,:),'Linewidth', 2);
% plot(Elt1(1,:),Elt1(7,:),'Linewidth', 2);
% plot(Elt4(1,:),Elt4(7,:),'Linewidth', 2);
% xlabel("Time (s)");
% ylabel("Stress (Pa)");
% title(EltLeg(7));
% legend(Legend3);

%Stress vs Strain
% figure; hold on; grid on;
% plot(Elt0(4,:),Elt0(7,:),'Linewidth', 2);
% plot(Elt1(4,:),Elt1(7,:),'Linewidth', 2);
% plot(Elt4(4,:),Elt4(7,:),'Linewidth', 2);
% xlabel("Strain");
% ylabel("Stress (Pa)");
% title("Stress Strain");
% legend(Legend3);

%% Node Outputs

Force0=Nodes0(6,:)+Nodes0(7,:)+Nodes0(8,:)+Nodes0(9,:);
NodesTot0=[Nodes0(1,:);Nodes0(10,:);Nodes0(14,:);Nodes0(2,:);Force0];
Force1=Nodes1(6,:)+Nodes1(7,:)+Nodes1(8,:)+Nodes1(9,:);
NodesTot1=[Nodes1(1,:);Nodes1(10,:);Nodes1(14,:);Nodes1(2,:);Force1];
Force4=Nodes4(6,:)+Nodes4(7,:)+Nodes4(8,:)+Nodes4(9,:);
NodesTot4=[Nodes4(1,:);Nodes4(10,:);Nodes4(14,:);Nodes4(2,:);Force4];

NodesLeg=["Time (s)", "Displacement", "Velocity","Acceleration","Applied Force"];

%Front Face Displacement
figure; hold on; grid on;
plot(NodesTot0(1,:),NodesTot0(2,:),'Linewidth', 2);
plot(NodesTot1(1,:),NodesTot1(2,:),'Linewidth', 2);
plot(NodesTot4(1,:),NodesTot4(2,:),'Linewidth', 2);
xlabel("Time (s)");
ylabel("Dispacement (m)");
title(NodesLeg(2));
legend(Legend3);

%Front Face Velocity
figure; hold on; grid on;
plot(NodesTot0(1,:),NodesTot0(3,:),'Linewidth', 2);
plot(NodesTot1(1,:),NodesTot1(3,:),'Linewidth', 2);
plot(NodesTot4(1,:),NodesTot4(3,:),'Linewidth', 2);
xlabel("Time (s)");
ylabel("Velocity (m/s)");
title(NodesLeg(3));
legend(Legend3);

%Front Face Acceleration 
figure; hold on; grid on;
plot(NodesTot0(1,:),NodesTot0(4,:),'Linewidth', 2);
plot(NodesTot1(1,:),NodesTot1(4,:),'Linewidth', 2);
plot(NodesTot4(1,:),NodesTot4(4,:),'Linewidth', 2);
xlabel("Time (s)");
ylabel("Acceleration (m/s^2)");
title(NodesLeg(4));
legend(Legend3);

% Applied Force (- of Reaction Force)
figure; hold on; grid on;
plot(NodesTot0(1,:),-NodesTot0(5,:),'Linewidth', 2);
plot(NodesTot1(1,:),-NodesTot1(5,:),'Linewidth', 2);
plot(NodesTot4(1,:),-NodesTot4(5,:),'Linewidth', 2);
xlabel("Time (s)");
ylabel("Applied Force (N)");
title(NodesLeg(5));
legend(Legend3);


%Force vs Strain
% figure; hold on; grid on;
% plot(Elt0(4,:),-NodesTot0(5,:),'Linewidth', 2);
% % plot(Elt0(4,:),-NodesTot0(5,:),'bo');
% plot(Elt1(4,:),-NodesTot1(5,:),'Linewidth', 2);
% % plot(Elt1(4,:),-NodesTot1(5,:),'bo');
% plot(Elt4(4,:),-NodesTot4(5,:),'Linewidth', 2);
% % plot(Elt4(4,:),-NodesTot4(5,:),'bo');
% xlabel("Log Strain");
% ylabel("Applied Force (N)");
% title("Reaction Force vs Log Strain");
% legend(Legend3);

%%

% Applied Force (- of Reaction Force) vs Displacment
figure; hold on; grid on;
plot(NodesTot0(2,:),-NodesTot0(5,:),'Linewidth', 2);
plot(NodesTot1(2,:),-NodesTot1(5,:),'Linewidth', 2);
plot(NodesTot4(2,:),-NodesTot4(5,:),'Linewidth', 2);
xlabel("Distance (m)");
ylabel("Applied Force (N)");
title("Applied Force vs Displacement");
legend(Legend3);

%Integrate Fdt and plot results
%Use Trapezodal Rule

LinearF=NodesTot0(1,:)*4;
ApE0 = zeros(1,size(NodesTot0,2));
ApE1 = zeros(1,size(NodesTot1,2));
ApE4 = zeros(1,size(NodesTot4,2));
for i=2:size(NodesTot0,2)
    ApE0(i)= ApE0(i-1)+(NodesTot0(2,i)-NodesTot0(2,i-1))*0.5*(-NodesTot0(5,i-1)-NodesTot0(5,i));
    ApE1(i)= ApE1(i-1)+(NodesTot1(2,i)-NodesTot1(2,i-1))*0.5*(-NodesTot1(5,i-1)-NodesTot1(5,i));
    ApE4(i)= ApE4(i-1)+(NodesTot4(2,i)-NodesTot4(2,i-1))*0.5*(-NodesTot4(5,i-1)-NodesTot4(5,i));
end
% Applied (Linear) Force vs Displacment
figure; hold on; grid on;
LinearF=NodesTot0(1,:)*4;
plot(NodesTot0(2,:),LinearF,'Color',Color0,'Linewidth', 2);
plot(NodesTot1(2,:),LinearF,'Color',Color1,'Linewidth', 2);
plot(NodesTot4(2,:),LinearF,'Color',Color4,'Linewidth', 2);
xlabel("Distance (m)");
ylabel("Applied (Linearized) Force (N)");
title("Applied (Linearized) Force vs Displacement");
legend(Legend3);

% Added Energy vs Displacement
figure; hold on; grid on;
plot(NodesTot0(2,:),ApE0,'Linewidth', 2);
plot(NodesTot1(2,:),ApE1,'Linewidth', 2);
plot(NodesTot4(2,:),ApE4,'Linewidth', 2);
xlabel("Distance (m)");
ylabel("Applied Energy (Nm)");
title("Applied Energy vs Displacement");
legend(Legend3);

% x1=find(ApE1==max(ApE1),1);
% plot(NodesTot0(2,:),ones(1,101)*max(ApE1),'k','Linewidth', 1);
% x0=find(ApE0>=max(ApE1),1);
% x4=find(ApE4>=max(ApE1),1);

%Strain Energy vs Displacement 
% figure; hold on; grid on;
% plot(NodesTot0(2,:),Energy0(4,:),'Linewidth', 2);
% plot(NodesTot1(2,:),Energy1(4,:),'Linewidth', 2);
% plot(NodesTot4(2,:),Energy4(4,:),'Linewidth', 2);
% xlabel("Distance (m)");
% ylabel("Strain Energy (J)");
% title("Strain Energy vs Displacement");
% legend(Legend3);;
%%
    %Energy Balance All Sims Truncated Displacement
% 
%     x1=find(max(Energy1(4,:))==(max(Energy1(4,:))),1);
%     x0=find(Energy0(4,:)>=max(Energy1(4,:)),1);
%     x4=find(Energy4(4,:)>=max(Energy1(4,:)),1);
% 
%     
%     figure; hold on; grid on;
%     plot(NodesTot0(2,1:x0),Energy0(4,1:x0),'Color', '#00008B' ,'Linewidth', 1,'MarkerIndices',1:5:100); 
%     plot(NodesTot0(2,1:x0),Energy0(4,1:x0)-Energy0(8,1:x0)+Energy0(5,1:x0)+Energy0(2,1:x0),'Color','#4169E1','Linewidth', 1,'MarkerIndices',1:5:100);
%     plot(NodesTot0(2,1:x0),-Energy0(8,1:x0),'Color', '#1E90FF','Linewidth', 1,'MarkerIndices',1:5:100);
%     
%     plot(NodesTot1(2,1:x1),Energy1(4,1:x1),'Color','#006400','Linewidth', 1,'MarkerIndices',1:5:100);
%     plot(NodesTot1(2,1:x1),Energy1(4,1:x1)-Energy1(8,1:x1)+Energy1(5,1:x1)+Energy1(2,1:x1),'Color',  '#228B22','Linewidth', 1,'MarkerIndices',1:5:100);
%     plot(NodesTot1(2,1:x1),-Energy1(8,1:x1),'Color','#3CB371' ,'Linewidth', 1,'MarkerIndices',1:5:100);
%     
%     plot(NodesTot4(2,1:x4),Energy4(4,1:x4),'Color','#B22222' ,'Linewidth', 1,'MarkerIndices',1:5:100);
%     plot(NodesTot4(2,1:x4),Energy4(4,1:x4)-Energy4(8,1:x4)+Energy4(5,1:x4)+Energy4(2,1:x4),'Color','#FF4500' ,'Linewidth', 1,'MarkerIndices',1:5:100);
%     plot(NodesTot4(2,1:x4),-Energy4(8,1:x4), 'Color','#FF6347','Linewidth', 1,'MarkerIndices',1:5:100);
%     
%     plot(NodesTot0(2,:),ones(1,101)*max(Energy1(4,1:x1)),'k','Linewidth', 1);
%     plot(NodesTot0(2,:),ones(1,101)*-max(Energy1(8,1:x1)),'k','Linewidth', 1);
%     
%     xlabel("Displacment (m)");
%     ylabel("Energy (J)");
%     title("Energy Balance");   
%     legend("0Fiber Strain Energy","0Fiber Total Energy","0Fiber External Work","102Fiber Strain Energy","102Fiber Total Energy","102Fiber External Work","4Fiber Strain Energy","4Fiber Total Energy","4Fiber External Work");
   
    
    %Strain Energy vs Displa
    figure; hold on; grid on;
    plot(NodesTot0(2,:),Energy0(4,:),'Color', '#00008B' ,'Linewidth', 2);
    plot(NodesTot1(2,:),Energy1(4,:),'Color','#006400','Linewidth', 2);
    plot(NodesTot4(2,:),Energy4(4,:),'Color','#B22222' ,'Linewidth', 2);
%     plot(NodesTot4(2,:),Energy4(4,1:2:201),'Color','#B22222' ,'Linewidth', 2);
    xlabel("Distance (m)");
    ylabel("Strain Energy (J)");
    title("Strain Energy vs Displacement");
    legend(Legend3);
    