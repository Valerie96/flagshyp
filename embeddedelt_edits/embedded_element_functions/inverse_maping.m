% function GEOM = inverse_maping(GEOM,FEM,tienodes)
load('InverseMappingTest.mat');
tienodes = BC.tienodes;
    e_elts = FEM.mesh.embedded;
    h_elts = FEM.mesh.host;
    e_nodes = tienodes(:);

    num_e = length(FEM.mesh.embedded);
    num_h = length(FEM.mesh.host);

%     NodeHost = zeros(length(e_nodes),2); NodeHost(:,1) = e_nodes;
%     ElementHost = zeros(num_e, 9); ElementHost(:,1) = FEM.mesh.embedded;
%     HostTotals = zeros(length(h_elts), 3); HostTotals(:,1) = FEM.mesh.host;
%     Zeta = zeros(4,length(e_nodes)); Zeta(1,:) = e_nodes;


    NodeHost = zeros(GEOM.npoin,2);
    ElementHost = zeros(FEM.mesh.nelem,9);
    HostTotals = zeros(GEOM.npoin,3);
    Zeta = zeros(4, GEOM.npoin);
    
    %Loop throgh host elements
    for j = 1:length(h_elts)
        h = h_elts(j);
        h_connectivity = FEM.mesh.connectivity(:,h);
        x_h = GEOM.x0(:,h_connectivity);  %Host node global coordinates

        %Loop over all embedded nodes
        for i = 1:length(e_nodes)
            ne = e_nodes(i);

            %Check if we already did this node
            if NodeHost(i,2) == 0
                x_ne = GEOM.x0(:,ne);
                inel = point_in_hexahedron(x_ne', x_h);

                %Check if it's in this host elt
                if inel
                    NodeHost(i,2) = h;
                    HostTotals(j,2) = HostTotals(j,2) + 1;

                    %Find natrual coordinates of the node
                    z = find(Zeta(1,:)==ne);
                    Zeta(1,z) = ne;
                    Zeta(2:4,z) = find_natural_coords(x_ne, x_h, FEM.mesh.element_type);
                end
            end
        end    
    end

    %Assign hosts to embedded elements based on first node
    for i=1:length(e_elts)
        e_connectivity = FEM.mesh.connectivity(:,e_elts(i));
        n1 = e_connectivity(1); %Choose first node
        nodeIndex=find(NodeHost==n1);
        host = NodeHost(nodeIndex,2);
        ElementHost(i,2:9) = host;

        hostIndex = (find(h_elts==host));
        HostTotals(hostIndex,3) = HostTotals(hostIndex,3) + 1;
    end

    GEOM.embedded.NodeHost = NodeHost;
    GEOM.embedded.ElementHost = ElementHost;
    GEOM.embedded.HostTotals = HostTotals;
    GEOM.embedded.Embed_Zeta = Zeta;

% end