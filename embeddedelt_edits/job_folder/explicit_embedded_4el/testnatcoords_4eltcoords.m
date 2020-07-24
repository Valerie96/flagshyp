clc; close all; 
% h_elets = FEM.mesh.connectivity(:,FEM.mesh.host);
% 
% e_elets = BC.tienodes;
% e_nodes = e_elets(:);
% xi = GEOM.x0(:,e_nodes(1));
% xn = GEOM.x0(:,h_elets(:,2));

xn = [1,1,0,0,1,1,0,0;0,1,1,0,0,1,1,0;0,0,0,0,1,1,1,1];


xn1c = mean(xn(1,:)); xn2c = mean(xn(2,:)); xn3c = mean(xn(3,:));
xi = x(:,3);
% xi(2)=xi(2)-0.5;
xi
% xi=[xi(1); -xi(2); xn3c]
% xi(3) = xi(3)-0.1;


inn =  point_in_hexahedron_test(xi',xn);


Zeta = find_natural_coords(xi, xn, 'hex')
X=find_xyz_in_host(Zeta,xn)