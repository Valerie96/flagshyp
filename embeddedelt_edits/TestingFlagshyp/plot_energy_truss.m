
% Read Flagshyp file

file=fopen('C:/Users/Valerie/Documents/GitHub/flagshyp/embeddedelt_edits/job_folder/truss_only/energy.dat','r');
formatSpec = '%e %e %e %e';
sizeA = [4 inf ];
A = fscanf(file,formatSpec,sizeA);
fclose(file);

fid = fopen('AbaqusTrussEnergy.dat');
    for i=1:4
        tline=fgetl(fid);
    end
    size =[5 Inf];
    info = fscanf(fid, '%f %f %f %f %f\n', size);
fclose(fid);


figure();
hold on; grid on;
plot(info(1,:),info(2,:),'DisplayName','Abaqus Internal Work','LineWidth',2)
plot(info(1,:),info(3,:),'DisplayName','Abaqus Kinetic Energy','LineWidth',2)
plot(info(1,:),-info(4,:),'DisplayName','Abaqus External Work','LineWidth',2)
plot(info(1,:),info(5,:),'DisplayName','Abaqus Total Energy','LineWidth',2)

plot(A(1,:),A(2,:),'DisplayName','Flagshyp Kinetic Energy','LineWidth',2)
plot(A(1,:),A(3,:),'DisplayName','Flagshyp Internal Work','LineWidth',2)
plot(A(1,:),-A(4,:),'DisplayName','Flagshyp External Work','LineWidth',2)
plot(A(1,:),A(2,:) + A(3,:)-A(4,:),'DisplayName','Flagshyp Total Energy','LineWidth',2)
legend('show')
ylabel('Energy(J)')
xlabel('Time (s)')
% 
figure();
hold on; grid on;
plot(A(1,:),A(2,:) + A(3,:)-A(4,:),'DisplayName','Total Energy','LineWidth',2)
plot(info(1,:),info(5,:),'DisplayName','Abaqus Total Energy','LineWidth',2)
legend('show')
ylabel('Energy (J)')
xlabel('Time (s)')
% 
% 
figure();
hold on; grid on;
plot(info(1,:),info(3,:),'DisplayName','Abaqus Kinetic Energy','LineWidth',2)
plot(A(1,:),A(2,:),'DisplayName','Kinetic Energy','LineWidth',2)
legend('show')
ylabel('Energy (J)')
xlabel('Time (s)')

figure();
hold on; grid on;
plot(info(1,:),info(2,:),'DisplayName','Abaqus Internal Work','LineWidth',2)
plot(A(1,:),A(3,:),'DisplayName','Flagshyp Internal Work','LineWidth',1)
legend('show')
ylabel('Energy (J)')
xlabel('Time (s)')



%% Plot Force Displacement
%Run Flagshyp first

fid=fopen('C:/Users/Valerie/Documents/GitHub/flagshyp/OneTruss.txt','r');
    for i=1:3
        tline=fgetl(fid);
    end
    size =[5 Inf];
    info = fscanf(fid, '%f %f %f %f %f\n', size);
    Atime   = info(1,:);
    AStrain = info(2,:);
    AForce  = info(3,:);
    AStress = info(4,:);
    ADisp  = info(5,:);    
fclose(fid);

ALamb = -ADisp + ones(1,length(ADisp)); 
A_J   = ALamb.^(1-2*0.3);

file=fopen('C:/Users/Valerie/Documents/GitHub/flagshyp/embeddedelt_edits/job_folder/truss_only/InternalForce.txt','r');
formatSpec = '%e %e %e';
sizeA = [7 inf ];
A = fscanf(file,formatSpec,sizeA);
fclose(file);

FForce = A(2,:);
FStress = A(4,:);
FJ = A(5,:);
Fa = A(6,:);
Fl = A(7,:);

FDisp = Fl-ones(1,length(Fl));

Area = 0.2;

figure();
hold on; grid on;
plot(AStress,-AForce,'DisplayName','Abaqus','LineWidth',2)
plot(FStress,FForce,'DisplayName','Flagshyp','LineWidth',2)
plot(AStress,AStress.*Area.*A_J./ALamb,'DisplayName','Abaqus Eq','LineWidth',2)
plot(FStress,FStress.*Area.*FJ./Fl,'DisplayName','Flagshyp Eq','LineWidth',1)
legend('show')
ylabel('Force (N)')
xlabel('Stress (Pa)')
hold off;

figure();
hold on; grid on;
plot(ALamb,-AForce,'o','DisplayName','Abaqus')
plot(Fl,FForce,'.','DisplayName','Flagshyp')
plot(ALamb,AStress.*Area.*A_J./ALamb,'DisplayName','Abaqus Eq','LineWidth',2)
plot(Fl,FStress.*Area.*FJ./Fl,'DisplayName','Flagshyp Eq','LineWidth',1)
plot([1 1.05], [0 -AForce(end)],'DisplayName','SomethingLinear', 'LineWidth',2)
legend('show')
ylabel('Force (N)')
xlabel('lambda')
hold off;

figure();
hold on; grid on;
plot(-ADisp,-AForce,'DisplayName','Abaqus','LineWidth',2)
plot(FDisp,FForce,'DisplayName','Flagshyp','LineWidth',2)
legend('show')
ylabel('Force (N)')
xlabel('Disp (m)')
hold off;


figure();
hold on; grid on;
plot(-ADisp,AStress,'DisplayName','Abaqus','LineWidth',2)
plot(FDisp,FStress,'DisplayName','Flagshyp','LineWidth',1)
legend('show')
ylabel('Stress')
xlabel('Disp (m)')
hold off;

FStress(end) - AStress(end)

% 
% f = fit(AStress',AForce1', 'poly1')
% fo = fitoptions('Method','NonlinearLeastSquares',...
%                'StartPoint',[0 0]);
% ft = fittype({'x', '0'},'options',fo);
% ff = fit(AStress',AForce1', ft)
% 
% plot( f, AStress,AForce1);
% f( 600 )