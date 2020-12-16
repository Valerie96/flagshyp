%number of dimensions
ndims = GEOM.ndime; 
massSize = ndims*GEOM.npoin;
LumpedMass = zeros(massSize,massSize);
n_nodes_elem = FEM.mesh.n_nodes_elem;

switch FEM.mesh.element_type
    case 'hexa8'
        Chi = [-1 1 1 -1 -1 1 1 -1; -1 -1 1 1 -1 -1 1 1; -1 -1 -1 -1 1 1 1 1];
        N = zeros(n_nodes_elem, n_nodes_elem);
        DN_chi = zeros(ndims,n_nodes_elem);

        for k = 1:n_nodes_elem
            interpolation = shape_functions_library(Chi(:,k),FEM.mesh.element_type);
            N(:,k)        = interpolation.N;
            DN_chi(:,k) = interpolation.DN_chi(:,k);
        end   
end

M_total = 0;
%Loop over all (host) elements
 ii=1;
    ielement=FEM.mesh.host(ii);
    %----------------------------------------------------------------------
    % Temporary variables associated with a particular element.
    %----------------------------------------------------------------------
    global_nodes    = FEM.mesh.connectivity(:,ielement);   
    material_number = MAT.matno(ielement);     
    matyp           = MAT.matyp(material_number);        
    properties      = MAT.props(:,material_number); % needs density included.
    xlocal          = GEOM.x(:,global_nodes);                     
    x0local         = GEOM.x0(:,global_nodes); 
    Ve              = GEOM.Ve(ielement);
   

    T = DN_chi'*x0local;
    
    M = eye(8)*2
    
    T*M*inv(T')