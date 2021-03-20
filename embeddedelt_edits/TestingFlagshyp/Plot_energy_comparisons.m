
% Read Flagshyp file

AbqOneHost = ReadHost;
AbqOneTruss = ReadTruss;
[AbqEHost, AbqETruss, AbqE] = ReadHostTruss();

%explicit_embedded_4elt_new
file=fopen('C:/Users/Valerie/Documents/GitHub/flagshyp/embeddedelt_edits/job_folder/explicit_embedded_truss/energy.dat','r');
formatSpec = '%e %e %e %e';
sizeA = [4 inf ];
TKIE = fscanf(file,formatSpec,sizeA);
fclose(file);

figure();
plot(TKIE(1,:),TKIE(2,:),'DisplayName','Kinetic Energy','LineWidth',2)
hold on; grid on;
plot(TKIE(1,:),TKIE(3,:),'DisplayName','Internal Work','LineWidth',2)
plot(TKIE(1,:),-TKIE(4,:),'DisplayName','External Work','LineWidth',2)
plot(TKIE(1,:),TKIE(2,:) + TKIE(3,:)-TKIE(4,:),'DisplayName','Total Energy','LineWidth',2)
legend('show')
ylabel('Energy(J)')
xlabel('Time (s)')

figure();
plot(AbqEHost.time, AbqE.KE,'DisplayName','Abaqus Kinetic Energy','LineWidth',2)
hold on; grid on;
plot(AbqEHost.time, AbqE.IE,'DisplayName','Abaqus Internal Work','LineWidth',2)
plot(AbqEHost.time, -AbqE.WK,'DisplayName','Abaqus External Work','LineWidth',2)
plot(AbqEHost.time, AbqE.ETOTAL,'DisplayName','Abaqus Total Energy','LineWidth',2)
legend('show')
ylabel('Energy(J)')
xlabel('Time (s)')

figure();
hold on; grid on;
plot(TKIE(1,:),TKIE(2,:) + TKIE(3,:)-TKIE(4,:),'DisplayName','Total Energy','LineWidth',2)
plot(AbqEHost.time, AbqE.ETOTAL,'DisplayName','Abaqus','LineWidth',1);
legend('show')
ylabel('Energy (J)')
xlabel('Time (s)')


figure();
hold on; grid on;
plot(TKIE(1,:),TKIE(2,:),'DisplayName','Kinetic Energy','LineWidth',2)
plot(AbqEHost.time, AbqE.KE,'DisplayName','Abaqus','LineWidth',1);
legend('show')
ylabel('Energy (J)')
xlabel('Time (s)')

figure();
hold on; grid on;
plot(TKIE(1,:),-TKIE(4,:),'DisplayName','External Work','LineWidth',2)
plot(AbqEHost.time, -AbqE.WK,'DisplayName','Abaqus','LineWidth',1);
legend('show')
ylabel('Energy (J)')
xlabel('Time (s)')

figure();
hold on; grid on;
plot(TKIE(1,:),TKIE(3,:),'DisplayName','Internal Work','LineWidth',2)
plot(AbqEHost.time, AbqE.IE,'DisplayName','Abaqus','LineWidth',1);
legend('show')
ylabel('Energy (J)')
xlabel('Time (s)')