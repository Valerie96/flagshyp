%--------------------------------------------------------------------------
% Check for equilibrium convergence.
%--------------------------------------------------------------------------
function [energy_value, max_energy] = check_energy_explicit(PRO,FEM,CON,BC,GLOBAL,...
       disp_n, disp_prev,fi_n,fi_prev,fe_n,fe_prev,time)

global Wint_n;
global Wext_n;

string='a';
%PRO.rest
%CON.incrm

if (~PRO.rest && CON.incrm==1)
    string='w';
    system('rm energy.dat');
    Wint_n = 0.0;
    Wext_n = 0.0;
end

string3=sprintf('energy.dat');

fid5= fopen(string3,string);

%internal work
sum = (disp_n(BC.hostdof) - disp_prev(BC.hostdof))' * (fi_n(BC.hostdof) + fi_prev(BC.hostdof));
Wint_n = Wint_n + 0.5 * sum;

%external work
sum = (disp_n(BC.hostdof) - disp_prev(BC.hostdof))' * (fe_n(BC.hostdof) + fe_prev(BC.hostdof));
Wext_n = Wext_n + 0.5 * sum;

% kinetic energy
WKE = 0;
% GLOBAL.velocities(BC.tiedof) = 0;
for i=1:FEM(1).mesh.n_dofs
    WKE = WKE + GLOBAL.M_KE(i,i) *( GLOBAL.velocities(i) * GLOBAL.velocities(i));
end
WKE = WKE * 0.5;

% for i=1:3:FEM(1).mesh.n_dofs
%     v2 = (GLOBAL.velocities(i)^2 + GLOBAL.velocities(i+1)^2 + GLOBAL.velocities(i+2)^2);
%     WKE = WKE + GLOBAL.M(i,i) *v2;
% end
% WKE = WKE * 0.5;

energy_value = abs(WKE + Wint_n - Wext_n);
numbers = [WKE, Wint_n, Wext_n];
max_energy = max(numbers);
energy_tolerance = 0.01;
if(energy_value > (energy_tolerance * max_energy))
    disp('Energy Violation!')
end

%--------------------------------------------------------------------------
% Print energy information.
%--------------------------------------------------------------------------
fprintf(fid5,'%5.5e %5.5e  %5.5e %5.5e \n',time, WKE, Wint_n, Wext_n);

fclose(fid5);

