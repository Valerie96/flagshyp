function GEOM = inverse_mapping(GEOM,FEM,tienodes)
% clear;clc;
% load('InverseMappingTest.mat');
% tienodes = BC.tienodes;
    e_elts = FEM.mesh.embedded;
    h_elts = FEM.mesh.host;
    e_nodes = tienodes(:);

    num_e = length(FEM.mesh.embedded);
    num_h = length(FEM.mesh.host);

%     NodeHost = zeros(length(e_nodes),2); NodeHost(:,1) = e_nodes;
%     ElementHost = zeros(num_e, 9); ElementHost(:,1) = FEM.mesh.embedded;
%     HostTotals = zeros(length(h_elts), 3); 
%     Zeta = zeros(4,length(e_nodes)); Zeta(1,:) = e_nodes;


    NodeHost = zeros(GEOM.npoin,1);
    ElementHost = zeros(FEM.mesh.nelem,8);
    HostTotals = zeros(FEM.mesh.nelem,2);
    Zeta = zeros(3,GEOM.npoin);
    
    %Loop throgh host elements
    for j = 1:length(h_elts)
        h = h_elts(j);
        h_connectivity = FEM.mesh.connectivity(:,h);
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
                    Zeta(:,ne) = find_natural_coords(x_ne, x_h, FEM.mesh.element_type);
                end
            end
        end    
    end

    %Assign hosts to embedded elements based on first node
    for i=1:length(e_elts)
        e = e_elts(i); 
        e_connectivity = FEM.mesh.connectivity(:,e);
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