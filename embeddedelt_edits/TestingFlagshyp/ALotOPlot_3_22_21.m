%Plotting Flagshyp vs Abaqus, Flagshyp vs Flagshyp
clear; clc; close all;
file1="explicit_embedded_truss";
name1 = "Flagshyp Uncorrected";

file2="embedded_truss_redundant_fixed";
name2 = "Flagshyp Corrected";

steps = 83;
nplot = 1;
numsteps = floor(steps/nplot);

graphsize=[100 100 800 400];

FLAG_1 = ReadFlagshypOutputFile(file1, 83,1);
FLAG_2 = ReadFlagshypOutputFile(file2, 82,1);

FLAG_o = ReadFlagshypOutputFile("explicit_3D", 82,1);

AbqOneHost = ReadHost;
AbqOneTruss = ReadTruss;
[AbqEHost, AbqETruss, AbqE] = ReadHostTruss('OneHostOneTrussResults');



%% Plot Embedded Flagshyp 1 vs Embedded Abaqus: Energy
file1="explicit_embedded_truss";
name1 = "Flagshyp Uncorrected";

FLAG_1 = ReadFlagshypOutputFile(file1, 83,1);
[AbqEHost, AbqETruss, AbqE] = ReadHostTruss('OneHostOneTrussResults'); 

figure();
plot(FLAG_1.Etime, FLAG_1.KE,'DisplayName','Flagshyp Kinetic Energy','LineWidth',2)
hold on; grid on;
plot(FLAG_1.Etime, FLAG_1.IE,'DisplayName','Flagshyp Internal Work','LineWidth',2)
plot(FLAG_1.Etime, FLAG_1.WK,'DisplayName','Flagshyp External Work','LineWidth',2)
plot(FLAG_1.Etime, FLAG_1.ET,'DisplayName','Flagshyp Total Energy','LineWidth',2)
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


PlotEnergy([FLAG_1.Etime, FLAG_1.KE], [AbqEHost.time, AbqE.KE], name1, 'Abaqus','Kinetic Energy')
PlotEnergy([FLAG_1.Etime, FLAG_1.IE], [AbqEHost.time, AbqE.IE], name1, 'Abaqus','Internal Energy')
PlotEnergy([FLAG_1.Etime, FLAG_1.WK], [AbqEHost.time, -AbqE.WK], name1, 'Abaqus','External Work')
PlotEnergy([FLAG_1.Etime, FLAG_1.ET], [AbqEHost.time, AbqE.ETOTAL], name1, 'Abaqus','Total Energy')

%% Plot Flagshyp 1 vs Abaqus: Feild Output

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.HostS(:,4),'DisplayName','Host YY Stress elt1','LineWidth',2);
plot(AbqEHost.time, AbqEHost.Stress(:,2),'DisplayName','Host Abaqus','LineWidth',2);
plot(FLAG_1.time,FLAG_1.TrussS(:,1),'DisplayName','Embedded Stress','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Stress(:,1),'DisplayName','Embedded Abaqus','LineWidth',1);
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.HostS(:,1),'DisplayName','Flagshyp','LineWidth',2);
plot(AbqEHost.time, AbqEHost.Stress(:,1),'DisplayName','Abaqus','LineWidth',2);
title("Host XX Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.TrussS(:,1),'DisplayName','Flagshyp','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Stress(:,1),'DisplayName','Abaqus','LineWidth',2);
title("Embedded XX Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.Disp(:,1,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(AbqEHost.time, AbqEHost.Displacement(:,1),'DisplayName','Host Abaqus','LineWidth',2);
plot(FLAG_1.time,FLAG_1.Disp(:,1,10),'DisplayName','Embedded Flagshyp','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Displacement(:,1),'DisplayName','Embedded Abaqus','LineWidth',2);
title("X Displacement");
xlabel("Time (s)");
ylabel("Displacement (m) ");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
% plot(FLAG_1.time,FLAG_1.HostLE(:,4,1),'DisplayName','Host Flagshyp','LineWidth',2);
% plot(AbqEHost.time, AbqEHost.Strain(:,2),'DisplayName','Host Abaqus','LineWidth',1);
plot(FLAG_1.time,FLAG_1.TrussLE(:,1,1),'DisplayName','Embedded Flagshyp','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Strain(:,1),'DisplayName','Embedded Abaqus','LineWidth',1);
title("XX Log Strain");
xlabel("Time (s)");
ylabel("Strain");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.HostLE(:,1,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(AbqEHost.time, AbqEHost.Strain(:,1),'DisplayName','Host Abaqus','LineWidth',1);
title("XX Log Strain: Host");
xlabel("Time (s)");
ylabel("Strain ");
legend('show');


figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.TrussLE(:,1,1),'DisplayName','Embedded Flagshyp','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Strain(:,1),'DisplayName','Embedded Abaqus','LineWidth',1);
title("XX Log Strain: Embedded Truss");
xlabel("Time (s)");
ylabel("Strain ");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
% plot(time,Displacements(:,1,1),'DisplayName','ux 1','LineWidth',2);
plot(FLAG_1.time,FLAG_1.Acc(:,1,1),'DisplayName','Flagshyp','LineWidth',2);

% plot(AbqEHost.time,AbqEHost.Displacement(:,1),'DisplayName','Abaqus ux 1','LineWidth',1);
plot(AbqEHost.time,AbqEHost.Acceleration(:,1),'DisplayName','Abaqus','LineWidth',1);
title("X Acceleration");
xlabel("Time (s)");
ylabel("Acceleration (m/s^2) ");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.RF(:,2,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(AbqEHost.time, AbqEHost.Force(:),'DisplayName','Host Abaqus','LineWidth',2);
plot(FLAG_1.time,FLAG_1.RF(:,2,10),'DisplayName','Embedded Flagshyp','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Force(:,1),'DisplayName','Embedded Abaqus','LineWidth',1);
title("Y Reaction Force");
xlabel("Time (s)");
ylabel("Reaction Force (N)");
legend('show');

%% Flagshyp 1 vs Flagshyp 2: Energy

file1="explicit_3D";
name1 = "Tension";

file2="explicit_wShear";
name2 = "Shear";

% file2="embedded_truss_redundant_fixed";
% name2 = "Flagshyp Corrected";

FLAG_1 = ReadFlagshypOutputFile(file1, 83,1); 
FLAG_2 = ReadFlagshypOutputFile(file2, 82,1);

PlotEnergy([FLAG_1.Etime, FLAG_1.KE], [FLAG_2.Etime, FLAG_2.KE], name1, name2,'Kinetic Energy')
PlotEnergy([FLAG_1.Etime, FLAG_1.IE], [FLAG_2.Etime, FLAG_2.IE], name1, name2,'Internal Energy')
PlotEnergy([FLAG_1.Etime, FLAG_1.WK], [FLAG_2.Etime, FLAG_2.WK], name1, name2,'External Work')
PlotEnergy([FLAG_1.Etime, FLAG_1.ET], [FLAG_2.Etime, FLAG_2.ET], name1, name2,'Total Energy')

%% Flagshyp 1 vs Flagshyp 2: Field Output

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.HostS(:,4),'DisplayName','Host YY Stress elt1','LineWidth',2);
plot(FLAG_2.time,FLAG_2.HostS(:,4),'DisplayName','Host YY Stress elt1','LineWidth',2);
% plot(FLAG_1.time,FLAG_1.TrussS(:,1),'DisplayName','Embedded Stress','LineWidth',2);
% plot(FLAG_2.time,FLAG_2.TrussS(:,1),'DisplayName','Embedded Stress','LineWidth',2);
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.HostS(:,1),'DisplayName','Flagshyp','LineWidth',2);
plot(FLAG_2.time,FLAG_2.HostS(:,1),'DisplayName','Flagshyp','LineWidth',2);
title("Host XX Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

% figure();
% hold on; grid on;
% % fig=gcf; fig.Position=graphsize;
% plot(FLAG_1.time,FLAG_1.TrussS(:,1),'DisplayName','Flagshyp','LineWidth',2);
% % plot(FLAG_2.time,FLAG_2.TrussS(:,1),'DisplayName','Flagshyp','LineWidth',2);
% title("Embedded XX Stress");
% xlabel("Time (s)");
% ylabel("Stress (Pa)");
% legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.Disp(:,1,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(FLAG_2.time,FLAG_2.Disp(:,1,1),'DisplayName','Host Flagshyp','LineWidth',2);
% plot(FLAG_1.time,FLAG_1.Disp(:,1,10),'DisplayName','Embedded Flagshyp','LineWidth',2);
% plot(FLAG_2.time,FLAG_2.Disp(:,1,10),'DisplayName','Embedded Flagshyp','LineWidth',2);
title("X Displacement");
xlabel("Time (s)");
ylabel("Displacement (m) ");
legend('show');

% figure();
% hold on; grid on;
% plot(FLAG_1.time,FLAG_1.TrussLE(:,1,1),'DisplayName','Embedded Flagshyp','LineWidth',2);
% % plot(FLAG_2.time,FLAG_2.TrussLE(:,1,1),'DisplayName','Embedded Flagshyp','LineWidth',2);
% title("XX Log Strain");
% xlabel("Time (s)");
% ylabel("Strain");
% legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.HostLE(:,1,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(FLAG_2.time,FLAG_2.HostLE(:,1,1),'DisplayName','Host Flagshyp','LineWidth',2);
title("XX Log Strain: Host");
xlabel("Time (s)");
ylabel("Strain ");
legend('show');


% figure();
% hold on; grid on;
% % fig=gcf; fig.Position=graphsize;
% plot(FLAG_1.time,FLAG_1.TrussLE(:,1,1),'DisplayName','Embedded Flagshyp','LineWidth',2);
% % plot(FLAG_2.time,FLAG_2.TrussLE(:,1,1),'DisplayName','Embedded Flagshyp','LineWidth',2);
% title("XX Log Strain: Embedded Truss");
% xlabel("Time (s)");
% ylabel("Strain ");
% legend('show');

figure();
hold on; grid on;
plot(FLAG_1.time,FLAG_1.Acc(:,1,1),'DisplayName','Flagshyp','LineWidth',2);
plot(FLAG_2.time,FLAG_2.Acc(:,1,1),'DisplayName','Flagshyp','LineWidth',2);
title("X Acceleration");
xlabel("Time (s)");
ylabel("Acceleration (m/s^2) ");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.RF(:,2,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(FLAG_2.time,FLAG_2.RF(:,2,1),'DisplayName','Host Flagshyp','LineWidth',2);
% plot(FLAG_1.time,FLAG_1.RF(:,2,10),'DisplayName','Embedded Flagshyp','LineWidth',2);
% plot(FLAG_2.time,FLAG_2.RF(:,2,10),'DisplayName','Embedded Flagshyp','LineWidth',2);
title("Y Reaction Force");
xlabel("Time (s)");
ylabel("Reaction Force (N)");
legend('show');

%% Abaqus vs Embedded Abaqus: Energy

[AbqEHost, AbqETruss, AbqE] = ReadHostTruss('OneHostOneTrussResults');
AbqOneHost = ReadHost;
graphsize=[100 100 800 400];
name3 = "Abaqus Embedded";
name4 = "Abaqus Solid";

PlotEnergy([AbqEHost.time, AbqE.KE], [AbqOneHost.time, AbqOneHost.KE], name3, name4,'Kinetic Energy')
PlotEnergy([AbqEHost.time, AbqE.IE], [AbqOneHost.time, AbqOneHost.IE], name3, name4,'Internal Energy')
PlotEnergy([AbqEHost.time, -AbqE.WK], [AbqOneHost.time, -AbqOneHost.WK], name3, name4,'External Work')
PlotEnergy([AbqEHost.time, AbqE.ETOTAL], [AbqOneHost.time, AbqOneHost.ETOTAL], name3, name4,'Total Energy')

%% One Host Elt Tests Flagshyp vs Abaqus: Energy
FLAG_o = ReadFlagshypOutputFile("explicit_3D", 82,1);
FLAG_o = ReadFlagshypOutputFile("embedded_truss_redundant_fixed", 83,1);
AbqOneHost = ReadHost;
graphsize=[100 100 800 400];
name3 = "Flagshyp Corrected";

PlotEnergy([FLAG_o.Etime, FLAG_o.KE], [AbqOneHost.time, AbqOneHost.KE], name3, 'Abaqus','Kinetic Energy')
PlotEnergy([FLAG_o.Etime, FLAG_o.IE], [AbqOneHost.time, AbqOneHost.IE], name3, 'Abaqus','Internal Energy')
PlotEnergy([FLAG_o.Etime, FLAG_o.WK], [AbqOneHost.time, -AbqOneHost.WK], name3, 'Abaqus','External Work')
PlotEnergy([FLAG_o.Etime, FLAG_o.ET], [AbqOneHost.time, AbqOneHost.ETOTAL], name3, 'Abaqus','Total Energy')

%% One Host Elt Tests Flagshyp vs Abaqus: Field Output
figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.HostS(:,4),'DisplayName','Host YY Stress elt1','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Stress(:,2),'DisplayName','Host Abaqus','LineWidth',2);
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.HostS(:,1),'DisplayName','Flagshyp','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Stress(:,1),'DisplayName','Abaqus','LineWidth',2);
title("Host XX Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.Disp(:,1,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Displacement(:,1),'DisplayName','Host Abaqus','LineWidth',2);
title("X Displacement");
xlabel("Time (s)");
ylabel("Displacement (m) ");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.HostLE(:,1,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Strain(:,1),'DisplayName','Host Abaqus','LineWidth',1);
title("XX Log Strain: Host");
xlabel("Time (s)");
ylabel("Strain ");
legend('show');


figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.Acc(:,1,1),'DisplayName','Flagshyp','LineWidth',2);
plot(AbqOneHost.time,AbqOneHost.Acceleration(:,1),'DisplayName','Abaqus','LineWidth',1);
title("X Acceleration");
xlabel("Time (s)");
ylabel("Acceleration (m/s^2) ");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.RF(:,2,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Force(:),'DisplayName','Host Abaqus','LineWidth',2);
title("Y Reaction Force");
xlabel("Time (s)");
ylabel("Reaction Force (N)");
legend('show');




%%
figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.HostS(:,4),'DisplayName','Host YY Stress elt1','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Stress(:,2),'DisplayName','Host Abaqus','LineWidth',2);
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.HostS(:,1),'DisplayName','Flagshyp','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Stress(:,1),'DisplayName','Abaqus','LineWidth',2);
title("Host XX Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.Disp(:,1,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Displacement(:,1),'DisplayName','Host Abaqus','LineWidth',2);
title("X Displacement");
xlabel("Time (s)");
ylabel("Displacement (m) ");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.HostLE(:,1,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Strain(:,1),'DisplayName','Host Abaqus','LineWidth',1);
title("XX Log Strain: Host");
xlabel("Time (s)");
ylabel("Strain ");
legend('show');


figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.Acc(:,1,1),'DisplayName','Flagshyp','LineWidth',2);
plot(AbqOneHost.time,AbqOneHost.Acceleration(:,1),'DisplayName','Abaqus','LineWidth',1);
title("X Acceleration");
xlabel("Time (s)");
ylabel("Acceleration (m/s^2) ");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_o.time,FLAG_o.RF(:,2,1),'DisplayName','Host Flagshyp','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Force(:),'DisplayName','Host Abaqus','LineWidth',2);
title("Y Reaction Force");
xlabel("Time (s)");
ylabel("Reaction Force (N)");
legend('show');


%% Flagshyp 1 vs Flagshyp 2 vs Flagshyp 3 Energy

file1="explicit_3D";
name1 = "Flagshyp No Truss";

file2="explicit_embedded_truss";
name2 = "Flagshyp Uncorrected";

file3="embedded_truss_redundant_fixed";
name3 = "Flagshyp Corrected";

FLAG_1 = ReadFlagshypOutputFile(file1, 83,1); 
FLAG_2 = ReadFlagshypOutputFile(file2, 82,1);
FLAG_3 = ReadFlagshypOutputFile(file3, 83,1);

PlotEnergy3([FLAG_1.Etime, FLAG_1.KE], [FLAG_2.Etime, FLAG_2.KE],[FLAG_3.Etime, FLAG_3.KE], name1, name2,name3,'Kinetic Energy')
PlotEnergy3([FLAG_1.Etime, FLAG_1.IE], [FLAG_2.Etime, FLAG_2.IE],[FLAG_3.Etime, FLAG_3.IE], name1, name2,name3,'Internal Energy')
PlotEnergy3([FLAG_1.Etime, FLAG_1.WK], [FLAG_2.Etime, FLAG_2.WK],[FLAG_3.Etime, FLAG_3.WK], name1, name2,name3,'External Work')
PlotEnergy3([FLAG_1.Etime, FLAG_1.ET], [FLAG_2.Etime, FLAG_2.ET],[FLAG_3.Etime, FLAG_3.ET], name1, name2,name3,'Total Energy')

%% Flagshyp 1 vs Flagshyp 2 vs Flagshyp 3: Field Output

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.HostS(:,4),'DisplayName',name1,'LineWidth',2);
plot(FLAG_2.time,FLAG_2.HostS(:,4),'DisplayName',name2,'LineWidth',2);
plot(FLAG_2.time,FLAG_2.HostS(:,4),'DisplayName',name3,'LineWidth',1);
% plot(FLAG_1.time,FLAG_1.TrussS(:,1),'DisplayName','Embedded Stress','LineWidth',2);
% plot(FLAG_2.time,FLAG_2.TrussS(:,1),'DisplayName','Embedded Stress','LineWidth',2);
xlabel("Time (s)");
ylabel("Host YY Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.HostS(:,1),'DisplayName',name1,'LineWidth',2);
plot(FLAG_2.time,FLAG_2.HostS(:,1),'DisplayName',name2,'LineWidth',2);
plot(FLAG_3.time,FLAG_3.HostS(:,1),'DisplayName',name3,'LineWidth',1);
title("Host XX Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_2.time,FLAG_2.TrussS(:,1),'DisplayName',name2,'LineWidth',2);
plot(FLAG_3.time,FLAG_3.TrussS(:,1),'DisplayName',name3,'LineWidth',1);
title("Embedded XX Stress");
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.Disp(:,1,1),'DisplayName',name1,'LineWidth',2);
plot(FLAG_2.time,FLAG_2.Disp(:,1,1),'DisplayName',name2,'LineWidth',2);
plot(FLAG_3.time,FLAG_3.Disp(:,1,1),'DisplayName',name3,'LineWidth',1);
% plot(FLAG_1.time,FLAG_1.Disp(:,1,10),'DisplayName','Embedded Flagshyp','LineWidth',2);
% plot(FLAG_2.time,FLAG_2.Disp(:,1,10),'DisplayName','Embedded Flagshyp','LineWidth',2);
title("Host X Displacement");
xlabel("Time (s)");
ylabel("Displacement (m) ");
legend('show');

figure();
hold on; grid on;
plot(FLAG_2.time,FLAG_2.TrussLE(:,1,1),'DisplayName',name2,'LineWidth',2);
plot(FLAG_3.time,FLAG_3.TrussLE(:,1,1),'DisplayName',name3,'LineWidth',1);
title("XX Log Strain: Embedded");
xlabel("Time (s)");
ylabel("Strain");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.HostLE(:,1,1),'DisplayName',name1,'LineWidth',2);
plot(FLAG_2.time,FLAG_2.HostLE(:,1,1),'DisplayName',name2,'LineWidth',2);
plot(FLAG_3.time,FLAG_3.HostLE(:,1,1),'DisplayName',name3,'LineWidth',1);
title("XX Log Strain: Host");
xlabel("Time (s)");
ylabel("Strain ");
legend('show');


figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_2.time,FLAG_2.TrussLE(:,1,1),'DisplayName',name2,'LineWidth',2);
plot(FLAG_3.time,FLAG_3.TrussLE(:,1,1),'DisplayName',name3,'LineWidth',1);
title("XX Log Strain: Embedded Truss");
xlabel("Time (s)");
ylabel("Strain ");
legend('show');

figure();
hold on; grid on;
plot(FLAG_1.time,FLAG_1.Acc(:,1,1),'DisplayName',name1,'LineWidth',2);
plot(FLAG_2.time,FLAG_2.Acc(:,1,1),'DisplayName',name2,'LineWidth',2);
plot(FLAG_3.time,FLAG_3.Acc(:,1,1),'DisplayName',name3,'LineWidth',1);
title("X Acceleration");
xlabel("Time (s)");
ylabel("Acceleration (m/s^2) ");
legend('show');

figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(FLAG_1.time,FLAG_1.RF(:,2,1),'DisplayName',name1,'LineWidth',2);
plot(FLAG_2.time,FLAG_2.RF(:,2,1),'DisplayName',name2,'LineWidth',2);
plot(FLAG_3.time,FLAG_3.RF(:,2,1),'DisplayName',name3,'LineWidth',1);
plot(FLAG_2.time,FLAG_2.RF(:,2,10),'DisplayName',name2,'LineWidth',2);
plot(FLAG_3.time,FLAG_3.RF(:,2,10),'DisplayName',name3,'LineWidth',1);
title("Y Reaction Force");
xlabel("Time (s)");
ylabel("Reaction Force (N)");
legend('show');


%% Abaqus vs Embedded Abaqus vs   Flagshyp 1 vs Flagshyp 2 vs Flagshyp 3 Energy

[AbqEHost, AbqETruss, AbqE] = ReadHostTruss('OneHostOneTrussResults');
AbqOneHost = ReadHost;
graphsize=[100 100 800 400];
name1a = "Abaqus Solid";
name2a = "Abaqus Embedded";

file1="explicit_3D";
name1 = "Flagshyp No Truss";

file2="explicit_embedded_truss";
name2 = "Flagshyp Uncorrected";

file3="embedded_truss_redundant_fixed";
name3 = "Flagshyp Corrected";

FLAG_1 = ReadFlagshypOutputFile(file1, 83,1); 
FLAG_2 = ReadFlagshypOutputFile(file2, 82,1);
FLAG_3 = ReadFlagshypOutputFile(file3, 83,1);

PlotEnergy5([AbqOneHost.time, AbqOneHost.KE],[AbqEHost.time, AbqE.KE], [FLAG_1.Etime, FLAG_1.KE], [FLAG_2.Etime, FLAG_2.KE],[FLAG_3.Etime, FLAG_3.KE], name1a, name2a,name1, name2,name3,'Kinetic Energy')
PlotEnergy5([AbqOneHost.time, AbqOneHost.IE],[AbqEHost.time, AbqE.IE], [FLAG_1.Etime, FLAG_1.IE], [FLAG_2.Etime, FLAG_2.IE],[FLAG_3.Etime, FLAG_3.IE], name1a, name2a,name1, name2,name3,'Internal Energy')
PlotEnergy5([AbqOneHost.time, -AbqOneHost.WK],[AbqEHost.time, -AbqE.WK],[FLAG_1.Etime, FLAG_1.WK], [FLAG_2.Etime, FLAG_2.WK],[FLAG_3.Etime, FLAG_3.WK],  name1a, name2a,name1, name2,name3,'External Work')
PlotEnergy5([AbqOneHost.time, AbqOneHost.ETOTAL],[AbqEHost.time, AbqE.ETOTAL], [FLAG_1.Etime, FLAG_1.ET], [FLAG_2.Etime, FLAG_2.ET],[FLAG_3.Etime, FLAG_3.ET], name1a, name2a,name1, name2,name3,'Total Energy')


figure();
hold on; grid on;
% fig=gcf; fig.Position=graphsize;
plot(AbqOneHost.time,AbqOneHost.Force,'bo','DisplayName',name1a);
plot(AbqEHost.time,AbqEHost.Force,'ro' ,'DisplayName',name2a);
plot(FLAG_1.time,FLAG_1.RF(:,2,1),'b','DisplayName',name1,'LineWidth',3);
plot(FLAG_2.time,FLAG_2.RF(:,2,1),'r','DisplayName',name2,'LineWidth',3);
plot(FLAG_3.time,FLAG_3.RF(:,2,1),'g','DisplayName',name3,'LineWidth',2);
title("External Reaction Force");
xlabel("Time (s)");
ylabel("Force (N)");
% legend('show');

%% Function Defs

function PlotEnergy(Data1, Data2, Name1, Name2,Title)
    figure();
    hold on; grid on;
    plot(Data1(:,1), Data1(:,2),'DisplayName',Name1,'LineWidth',2)
    plot(Data2(:,1), Data2(:,2),'DisplayName',Name2,'LineWidth',2)
    legend('show')
    title(Title);
    ylabel('Energy(J)')
    xlabel('Time (s)')
end

function PlotEnergy3(Data1, Data2, Data3, Name1, Name2,Name3,Title)
    figure();
    hold on; grid on;
    plot(Data1(:,1), Data1(:,2),'DisplayName',Name1,'LineWidth',2)
    plot(Data2(:,1), Data2(:,2),'DisplayName',Name2,'LineWidth',2)
    plot(Data3(:,1), Data3(:,2),'DisplayName',Name3,'LineWidth',1)
    legend('show')
    title(Title);
    ylabel('Energy(J)')
    xlabel('Time (s)')
end

function PlotEnergy5(Data1, Data2, Data3, Data4,Data5, Name1, Name2,Name3,Name4,Name5,Title)
    figure();
    hold on; grid on;
    plot(Data1(:,1), Data1(:,2),'bo','DisplayName',Name1)
    plot(Data2(:,1), Data2(:,2),'ro','DisplayName',Name2)
    plot(Data3(:,1), Data3(:,2),'b','DisplayName',Name3,'LineWidth',3)
    plot(Data4(:,1), Data4(:,2),'r','DisplayName',Name4,'LineWidth',2)
    plot(Data5(:,1), Data5(:,2),'g','DisplayName',Name5,'LineWidth',2)
    legend('show')
    title(Title);
    ylabel('Energy(J)')
    xlabel('Time (s)')
end