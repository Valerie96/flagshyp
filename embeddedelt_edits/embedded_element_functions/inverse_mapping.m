function GEOM = inverse_mapping(GEOM,FEM,tienodes)
% Calculates the host element coordinates of embedded elements and creates
% NodeHost: list of nodes and their associated host, and ElementHost: list
% of elements and their associated host. 

% clear;clc;
% load('InverseMappingTest.mat');
% tienodes = BC.tienodes;
%     e_elts = FEM.mesh.embedded;
%     h_elts = FEM.mesh.host;
    e_nodes = tienodes(:);
%-----------------------------------------------------------------------
%!!!!!!!!!!!
%For now, assume all elements of the first type are hosts and all of the
%second type are embedded
%!!!!!!!!!!!
%-----------------------------------------------------------------------
    h_elts=FEM(1).mesh.host; 
    e_elts=FEM(2).mesh.embedded;
    Global_nums = 1:GEOM.total_n_elets;



    NodeHost = zeros(GEOM.npoin,1);
    ElementHost = zeros(GEOM.total_n_elets,8);
    HostTotals = zeros(GEOM.total_n_elets,2);
    Zeta = zeros(3,GEOM.npoin);
    
    %Loop throgh host elements
    for j = 1:length(h_elts)
        h = h_elts(j);
        h_connectivity = FEM(1).mesh.connectivity(:,h);
        x_h = GEOM.x0(:,h_connectivity);  %Host node global coordinates

        %Loop over all embedded nodes
        for i = 1:length(e_nodes)
            ne = e_nodes(i);

            %Check if we already did this node
            if NodeHost(ne) == 0
                x_ne = GEOM.x0(:,ne);
                inel = point_in_hexahedron(x_ne', x_h);

                %Check if it's in this host elt
                if inel
                    NodeHost(ne) = h;
                    HostTotals(h,1) = HostTotals(h,1) + 1;

                    %Find natrual coordinates of the node
                    Zeta(:,ne) = find_natural_coords(x_ne, x_h, FEM(1).mesh.element_type);
                end
            end
        end    
    end

    %Assign hosts to embedded elements based on first node
    for i=1:length(e_elts)
        e = e_elts(i); 
        e_connectivity = FEM(2).mesh.connectivity(:,e);
        n1 = e_connectivity(1); %Choose first node
        host = NodeHost(n1);
        ElementHost(e,:) = host;
        HostTotals(host,2) = HostTotals(host,2) + 1;
    end

    GEOM.embedded.NodeHost = NodeHost;
    GEOM.embedded.ElementHost = ElementHost;
    GEOM.embedded.HostTotals = HostTotals;
    GEOM.embedded.Embed_Zeta = Zeta;

end