%Calculates the embedded element effect on host element during
%InternalForce_explicit


function [T_internal] = DistributeCorrectInternalForce_explicit2(ielement,...
          T_internal,FEM,xlocal,x0local,QUADRATURE,CONSTANT,GEOM,GlobT_int,...
          PLAST,KINEMATICS,MAT,DAMPING,eelt)


global explicit
dim=GEOM.ndime;

global VolumeCorrect;

% *** For now we are assuming that host elements are of FEM type 1 and
% embedded are FEM type 2
    %----------------------------------------------------------------------
    % GATHER material properties of the host and embedded elt
    %----------------------------------------------------------------------   
    material_number   = MAT(1).matno(ielement);     
    matyp_h           = MAT(1).matyp(material_number);        
    properties_h      = MAT(1).props(:,material_number);                       
    Ve_h              = GEOM.Ve(ielement);  
    
    material_number   = MAT(2).matno(eelt);     
    matyp_e           = MAT(2).matyp(material_number);        
    properties_e      = MAT(2).props(:,material_number);
    Ve_e              = GEOM.Ve(eelt);  
    
    
    switch matyp_h
     case {5,7,17}
          [pressure_h,kappa_bar_h,DN_x_mean_h,ve_h] = ...
           mean_dilatation_pressure(FEM,dim,matyp_h,properties_h,Ve_h,...
                                    QUADRATURE,KINEMATICS);
     otherwise
          pressure_h = 0;
    end
    
    switch matyp_e
     case {5,7,17}
          [pressure_e,kappa_bar_e,DN_x_mean_e,ve_e] = ...
           mean_dilatation_pressure(FEM,dim,matyp_e,properties_e,Ve_e,...
                                    QUADRATURE_eh,KINEMATICS_eh);
     otherwise
          pressure_e = 0;
    end



    %--------------------------------------------------------------------------
    %Step 1 
    %Get the embedded element quadrature points in the host element domain
    %---------------------------------------------------------------------------
    
    %Get the number of Gauss points in the embedded element
    nGp = GEOM.embedded.HostTotals(ielement,1)/GEOM.embedded.HostTotals(ielement,2);
    
    %Get host element nodes
    h_connectivity = FEM(1).mesh.connectivity(:,ielement);
    x_h = GEOM.x0(:,h_connectivity);  %Host node global coordinates
    
    %Get embedded element information
    e_connectivity = FEM(2).mesh.connectivity(:,eelt);
    x_e = GEOM.x0(:,e_connectivity);
    xelocal  = GEOM.x(:,e_connectivity);                     
    e_nodes_zeta = GEOM.embedded.Embed_Zeta(:,e_connectivity);
    
    QUADRATURE_E = QUADRATURE(2); %Quadrature of embedded elt
    QUADRATURE_EH = QUADRATURE(2); %Quadrature of embedded elt in the host domain 
    
    %Initialize KINEMATICS_EH to use the (8) shape functions of the host
    %and the (2) gauss points of the embedded element
        % Spatial gradient of the shape functions.
        KINEMATICS_EH.DN_x   = zeros(GEOM.ndime,FEM(1).mesh.n_nodes_elem,nGp);  
        % Jacobian of the mapping between spatial and isoparametric domains.
        KINEMATICS_EH.Jx_chi = zeros(nGp,1);  
        % Deformation gradient.
        KINEMATICS_EH.F      = zeros(GEOM.ndime,GEOM.ndime, nGp);                 
        % Jacobian of the deformation gradient.
        KINEMATICS_EH.J      = zeros(nGp ,1);            
        % Left Cauchy-Green strain tensor (b).
        KINEMATICS_EH.b      = zeros(GEOM.ndime,GEOM.ndime, nGp);   
        % First invariant of b.
        KINEMATICS_EH.Ib     = zeros(nGp,1);      
        % Principal stretches.
        KINEMATICS_EH.lambda = zeros(GEOM.ndime,nGp);                
        % Spatial principal directions.
        KINEMATICS_EH.n      = zeros(GEOM.ndime,GEOM.ndime,nGp);
        
    KINEMATICS_E = KINEMATICS(2);

    gp_x = zeros(3,2); 
    gp_zeta = zeros(3,2);  
    for gp = 1: nGp
              
        %Get the coordinates of the gauss point in the embedded space (LE)
        gp_nu = QUADRATURE_E.element.Chi(gp,:);

        %Find Gauss points in x
         N_nu = shape_function_values_at(gp_nu, 'truss2'); %mapping with shape functions
         for i=1:length(N_nu)
           gp_x(:,gp) = gp_x(:,gp) + N_nu(i,1)*x_e(:,i); 
         end
        
         %Find Gauss points in zeta
         gp_zeta(:,gp) = find_natural_coords(gp_x(:,gp), x_h, FEM(1).mesh.element_type);
        
        
        
        %Test accuracy of that conversion 
            %Find gp_x from gp_nu and gp_zeta
            gp_xz = find_xyz_in_host(gp_zeta(:,gp), x_h);
            gp_xn = find_xyz_in_truss(gp_nu, x_e);
            check = gp_xz - gp_xn;
            if abs(check(1))>1E-10 || abs(check(2))>1E-10 || abs(check(3)) >1E-10
                fprintf("Space conversion failure: host %u, guest %u, gp %u\n",ielement, eelt, gp);
                fprintf("     Error amount: %d %d %d\n", abs(check(1)), abs(check(2)), abs(check(3)));
            end  
            

 
    end
 
    
    %----------------------------------------------------------------------
    %Step B 
    %Compute deformation measures at embedded element quad points in host
    %element space
    %----------------------------------------------------------------------
    
    QUADRATURE_EH.element.Chi = gp_zeta;
    
Ne      = zeros(FEM(1).mesh.n_nodes_elem, nGp);
DNe_chi = zeros(GEOM.ndime,FEM(1).mesh.n_nodes_elem, nGp);

for igauss = 1:nGp 
    interpolation       = shape_functions_library(gp_zeta(:,igauss),...
                                                   FEM(1).mesh.element_type);
    Ne(:,igauss)        = interpolation.N;
    DNe_chi(:,:,igauss) = interpolation.DN_chi;
end
    
INTERPOLATION_EH.N = Ne;
INTERPOLATION_EH.DN_chi = DNe_chi; 


    KINEMATICS_EH = gradients(xlocal,x0local,INTERPOLATION_EH.DN_chi,...
             QUADRATURE_EH.element,KINEMATICS_EH); 
         
    KINEMATICS_A = gradients(xlocal,x0local,FEM(1).interpolation.element.DN_chi,...
             QUADRATURE(1).element,KINEMATICS(1)); %This is just a reference to what kinematics for the host element are normaly
         
    area_t = properties_e(4);   
    
    %If Host Element is Neo Hookean (need to change if otherwise)
    lam=properties_h(3);
    mu=properties_h(2);
    K=lam+(2*mu/3);
    nu_h   = (3*K-2*mu)/(2*(3*K+mu));
    
    KINEMATICS_E = gradientsTruss(xelocal,x_e,FEM(2).interpolation.element.DN_chi,...
             QUADRATURE_E.element,KINEMATICS_E,nu_h,area_t);      
%--------------------------------------------------------------------------
%Find Embedded element orientation wrt the global csys
L       = norm(x_e(:,2) - x_e(:,1));  
dx      = xelocal(:,2) - xelocal(:,1);        
l       = norm(dx);                           
n       = dx/l;  

%--------------------------------------------------------------------------
% Gauss quadrature integration loop. Loop of embedded element Gauss Points
%--------------------------------------------------------------------------    
% T_C = zeros(FEM(1).mesh.n_dofs_elem, 1);
TC = zeros(FEM(2).mesh.n_dofs_elem, 1);
    for igauss = 1:nGp
        
        %Step B
        %----------------------------------------------------------------------
        % Extract kinematics at the particular Gauss point.
        %----------------------------------------------------------------------
        kinematics_gauss = kinematics_gauss_point(KINEMATICS_E,igauss);   
        

        
        %Step D
        %----------------------------------------------------------------------
        % Calculate embeddded element stress measure in the host elt system
        % using the host elt material model (ie correction stress)
        %----------------------------------------------------------------------    
        [Cauchy_C,PLAST_h,...
         plast_gauss] = Cauchy_type_selection(kinematics_gauss,properties_h,...
                                              CONSTANT,dim,matyp_h,PLAST,igauss);
%         %----------------------------------------------------------------------
%         % Obtain elasticity tensor (for incompressible or nearly incompressible, 
%         % only deviatoric component).
%         %----------------------------------------------------------------------   
%         if(explicit==0)
%             c = elasticity_modulus_selection(kinematics_gauss,properties_h,CONSTANT,...
%                                           dim,matyp_h,PLAST_h,plast_gauss,igauss);
%         else
%             c = 0;
%         end
        %----------------------------------------------------------------------
        % Add pressure contribution to stresses and elasticity tensor.
        %----------------------------------------------------------------------    
%         [Cauchy_C,c] = mean_dilatation_pressure_addition(Cauchy_C,c,CONSTANT,pressure_h,matyp_h);    
        %----------------------------------------------------------------------
        
        
        %----------------------------------------------------------------------
        % Compute numerical integration multipliers.
        %----------------------------------------------------------------------
        JW = kinematics_gauss.Jx_chi*QUADRATURE_E.element.W(igauss)*...
             thickness_plane_stress(properties_h,kinematics_gauss.J,matyp_h);
        %----------------------------------------------------------------------
        % Compute contribution to (internal) force vector.
        %----------------------------------------------------------------------
        T = Cauchy_C*kinematics_gauss.DN_x;
%         T_C = T_C + T(:)*JW;
        TC = TC + T(:)*JW;
        

        
    end
    
                
        % Step C
        %----------------------------------------------------------------------
        % Get embeddded element internal force and convert to force on host
        % nodes
        %---------------------------------------------------------------------- 
        edof = FEM(1).mesh.dof_nodes(:,e_connectivity);
        Tint_e = GlobT_int(edof);
        N_node1 = shape_function_values_at(e_nodes_zeta(:,1), FEM(1).mesh.element_type);
        N_node2 = shape_function_values_at(e_nodes_zeta(:,2), FEM(1).mesh.element_type);
        
        %Force from embedded nodes distribted over host nodes
        T_e1 = zeros(GEOM.ndime*8, 1);
        T_e2 = zeros(GEOM.ndime*8, 1);
        for i = 1:3:24
           T_e1(i:i+2) = Tint_e(:,1)*N_node1((i-1)/3 + 1); 
           T_e2(i:i+2) = Tint_e(:,2)*N_node2((i-1)/3 + 1); 
        end
        T_e = T_e1 + T_e2;
    
        
       %Do the same for the correction force
        T_C1 = zeros(GEOM.ndime*8, 1);
        T_C2 = zeros(GEOM.ndime*8, 1);
        for i = 1:3:24
           T_C1(i:i+2) = TC(1:3)*N_node1((i-1)/3 + 1); 
           T_C2(i:i+2) = TC(4:6)*N_node2((i-1)/3 + 1); 
        end
        T_C = T_C1 + T_C2;
        
        
        %Step E
        %----------------------------------------------------------------------
        % Compute equivilant (internal) force vector of the host element.
        %----------------------------------------------------------------------
        if VolumeCorrect
            T_internal = T_internal + (T_e - T_C);
        else
            T_internal = T_internal + (T_e);
        end
 
    
    
%     Step_globalT_int = force_vectors_assembly(T_internal,global_nodes,...
%                    zeros(FEM.mesh.n_dofs,1),FEM.mesh.dof_nodes);
      
end