%Plot Stress
clear; clc; close all;
name="explicit_embedded_truss";
steps = 83;
nplot = 1;

basedir=strcat('C:/Users/Valerie/Documents/GitHub/flagshyp/embeddedelt_edits/job_folder/',name);
set(0,'defaultfigurecolor',[1 1 1]);
file = strcat(basedir,'/', name, '-OUTPUT.txt');
file =fopen(file, 'r');
numsteps = floor(steps/nplot);
graphsize=[100 100 800 400];


%loop once to get initial values, num elt, num nodes, elt type. 
%loop though rest of steps with tha data 

[numbers,n_nodes,nodeinfo1,nodeinfo2,nelttype,n_elt,npts,...
    ncomp,eltinfo1,eltinfo2] = getInitialStep(file);

time = zeros(numsteps,1);
dt = zeros(numsteps,1);

Displacements = zeros(numsteps, 3, n_nodes);
Forces = zeros(numsteps, 3, n_nodes);
Velocity = zeros(numsteps, 3, n_nodes);
Acceleration = zeros(numsteps, 3, n_nodes);

Stress1 = zeros(numsteps, ncomp(1), n_elt(1));
Strain1 = zeros(numsteps, ncomp(1), n_elt(1));
Stress2 = zeros(numsteps, ncomp(2), n_elt(2));
Strain2 = zeros(numsteps, ncomp(2), n_elt(2));

time(1) = numbers(3);
dt(1) = numbers(4);

for i = 1:n_nodes
    Forces(1,:,i) = nodeinfo1(i,6:8);
    Displacements(1,:,i) = nodeinfo2(i,2:4);
    Velocity(1,:,i) = nodeinfo2(i,5:7);
    Acceleration(1,:,i) = nodeinfo2(i,8:10);
end


for i = 1:n_elt(1)
    Stress1(1, :, i) = mean(eltinfo1(:, 1:6),1);
    Strain1(1, :, i) = mean(eltinfo1(:, 7:12),1);
end
for i = 1:n_elt(2)
    Stress2(1, :, i) = eltinfo2(i,1);
    Strain2(1, :, i) = eltinfo2(i,2);
end

%Read the rest of the steps. This could be simplified by not looking for
%the number of nodes/elements etc every time, but we can do that later

for k = 2:numsteps
    [numbers,n_nodes,nodeinfo1,nodeinfo2,nelttype,n_elt,npts,...
    ncomp,eltinfo1,eltinfo2] = getInitialStep(file);
    

    time(k) = numbers(3);
    dt(k) = numbers(4);

    for i = 1:n_nodes
        Forces(k,:,i) = nodeinfo1(i,6:8);
        Displacements(k,:,i) = nodeinfo2(i,2:4);
        Velocity(k,:,i) = nodeinfo2(i,5:7);
        Acceleration(k,:,i) = nodeinfo2(i,8:10);
    end


    for i = 1:n_elt(1)
        Stress1(k, :, i) = mean(eltinfo1(:, 1:6),1);
        Strain1(k, :, i) = mean(eltinfo1(:, 7:12),1);
    end
    for i = 1:n_elt(2)
        Stress2(k, :, i) = eltinfo2(i,1);
        Strain2(k, :, i) = eltinfo2(i,2);
    end
end

fclose(file);

AbqOneHost = ReadHost;
AbqOneTruss = ReadTruss;
[AbqEHost, AbqETruss, AbqE] = ReadHostTruss();
%% Plots for Host and Truss

figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
plot(time,Stress1(:,4),'DisplayName','Host YY Stress elt1','LineWidth',2);
plot(AbqEHost.time, AbqEHost.Stress(:,2),'DisplayName','Host Abaqus','LineWidth',2);
plot(time,Stress2(:,1),'DisplayName','Embedded Stress','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Stress(:,1),'DisplayName','Embedded Abaqus','LineWidth',1);
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
plot(time,Stress1(:,1),'DisplayName','Host XX Stress elt1','LineWidth',2);
plot(AbqEHost.time, AbqEHost.Stress(:,1),'DisplayName','Abaqus','LineWidth',2);
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
plot(time,Stress2(:,1),'DisplayName','Embedded XX Stress','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Stress(:,1),'DisplayName','Abaqus','LineWidth',2);
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
plot(time,Displacements(:,1,1),'DisplayName','Host ux 1','LineWidth',2);
plot(AbqEHost.time, AbqEHost.Displacement(:,1),'DisplayName','Host Abaqus','LineWidth',2);
plot(time,Displacements(:,1,10),'DisplayName','Embedded ux 1','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Displacement(:,1),'DisplayName','Embedded Abaqus','LineWidth',2);
xlabel("Time (s)");
ylabel("Displacement (m) ");
legend('show');

figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
% plot(time,Strain1(:,4,1),'DisplayName','Host LEy 1','LineWidth',2);
% plot(AbqEHost.time, AbqEHost.Strain(:,2),'DisplayName','Host Abaqus','LineWidth',1);
plot(time,Strain2(:,1,1),'DisplayName','Embedded LEx 1','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Strain(:,1),'DisplayName','Embedded Abaqus','LineWidth',1);
xlabel("Time (s)");
ylabel("Strain");
legend('show');

figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
plot(time,Strain1(:,1,1),'DisplayName','Host LEx 1','LineWidth',2);
plot(AbqEHost.time, AbqEHost.Strain(:,1),'DisplayName','Host Abaqus','LineWidth',1);
xlabel("Time (s)");
ylabel("Strain ");
legend('show');


figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
plot(time,Strain2(:,1,1),'DisplayName','Embedded LEx 1','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Strain(:,1),'DisplayName','Embedded Abaqus','LineWidth',1);
xlabel("Time (s)");
ylabel("Strain ");
legend('show');

figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
% plot(time,Displacements(:,1,1),'DisplayName','ux 1','LineWidth',2);
plot(time,Acceleration(:,1,1),'DisplayName','ax 1','LineWidth',2);

% plot(AbqEHost.time,AbqEHost.Displacement(:,1),'DisplayName','Abaqus ux 1','LineWidth',1);
plot(AbqEHost.time,AbqEHost.Acceleration(:,1),'DisplayName','Abaqus ax 1','LineWidth',1);
xlabel("Time (s)");
ylabel("Acceleration (m/s^2) ");
legend('show');

figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
plot(time,Forces(:,2,1),'DisplayName','Host RFy elt1','LineWidth',2);
plot(AbqEHost.time, AbqEHost.Force(:),'DisplayName','Host Abaqus','LineWidth',2);
plot(time,Forces(:,2,10),'DisplayName','Embedded RFy','LineWidth',2);
plot(AbqEHost.time, AbqETruss.Force(:,1),'DisplayName','Embedded Abaqus','LineWidth',1);
xlabel("Time (s)");
ylabel("Reaction Force (N)");
legend('show');



%% Plots for just host
figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
plot(time,Stress1(:,1),'DisplayName','XX Stress elt1','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Stress(:,1),'DisplayName','Abaqus','LineWidth',2);
xlabel("Time (s)");
ylabel("Stress (Pa)");
legend('show');

figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
plot(time,Displacements(:,1,1),'DisplayName','ux 1','LineWidth',2);
plot(AbqOneHost.time, AbqOneHost.Displacement(:,1),'DisplayName','Abaqus','LineWidth',2);
xlabel("Time (s)");
ylabel("stuff ");
legend('show');

figure();
hold on; grid on;
fig=gcf; fig.Position=graphsize;
plot(time,Displacements(:,1,1),'DisplayName','ux 1','LineWidth',2);
% plot(time,Velocity(:,1,1),'DisplayName','vx 1','LineWidth',2);
plot(time,Acceleration(:,1,1),'DisplayName','ax 1','LineWidth',2);

plot(AbqOneHost.time,AbqOneHost.Displacement(:,1),'DisplayName','Abaqus ux 1','LineWidth',1);
% plot(AbqHost.time,AbqHost.Velocity(:,1),'DisplayName','Abaqus vx 1','LineWidth',1);
plot(AbqOneHost.time,AbqOneHost.Acceleration(:,1),'DisplayName','Abaqus ax 1','LineWidth',1);
xlabel("Time (s)");
ylabel("stuff ");
legend('show');
%%

function [numbers,n_nodes,nodeinfo1,nodeinfo2,nelttype,n_elt,npts,...
    ncomp,eltinfo1,eltinfo2] = getInitialStep(fid)

    %Read time step info
    text = fgetl(fid);
    numbers = sscanf(text, "3-D EmbeddedElt       at increment:       %d,        load:  %d       time:  %f       dt:  %f");
    text = fgetl(fid);
    n_nodes = sscanf(text, "%d");

    trash = fgetl(fid);
    trash = fgetl(fid);

    nodeinfo1 = zeros(n_nodes, 8);
    for i = 1:n_nodes
        text = fgetl(fid);
        nodeinfo1(i,:) = sscanf(text, "%d %d %f %f %f %f %f %f");
    end

    trash = fgetl(fid);
    trash = fgetl(fid);
    trash = fgetl(fid);

    nodeinfo2 = zeros(n_nodes, 10);
    for i = 1:n_nodes
        text = fgetl(fid);
        nodeinfo2(i,:) = sscanf(text, "%d %f %f %f %*c %*c %f %f %f %*c %*c %f %f %f");
    end

    trash = fgetl(fid);

    text = fgetl(fid);
    nelttype = sscanf(text,"Element Types: %d");
    trash = fgetl(fid);

    n_elt = zeros(1,nelttype+1);
    npts = zeros(1,nelttype+1);
    ncomp = zeros(1,nelttype+1);
    eltinfo1 = zeros(2); %garbage values until element type is defined
    eltinfo2 = zeros(2);

    for xx = 1:nelttype
        eltype = fgetl(fid);
        switch eltype
            case 'hexa8'
                npts(xx) = 8;
                ncomp(xx) = 6;
                text = fgetl(fid);
                n_elt(xx) = sscanf(text,"Elements: %d");
                eltinfo1 = zeros(8*n_elt(xx),12);
            case 'truss2'
                npts(xx) = 1;
                ncomp(xx) = 1;
                text = fgetl(fid);
                n_elt(xx) = sscanf(text,"Elements: %d");
                eltinfo2 = zeros(1*n_elt(xx),2);
        end


        trash = fgetl(fid);
        trash = fgetl(fid);


        switch eltype
            case 'hexa8'
                for i = 1:n_elt(xx)
                   for j = 1:npts(xx)
                    text = fgetl(fid);
                    eltinfo1(j,:) = sscanf(text, "%f %f %f %f %f %f %*c %*c %f %f %f %f %f %f"); 
                   end 
                   trash = fgetl(fid);
                   trash = fgetl(fid);
                end
            case 'truss2'
                 for i = 1:n_elt(xx)
                   for j = 1:npts(xx)
                    text = fgetl(fid);
                    eltinfo2(j,:) = sscanf(text, "%f %*c %*c %f"); 
                   end 
                end
        end
    end

    for i=1:4
        trash = fgetl(fid);
    end

end


