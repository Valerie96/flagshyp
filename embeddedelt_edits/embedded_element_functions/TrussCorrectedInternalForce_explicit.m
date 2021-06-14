%Calculates the embedded element effect on host element during
%InternalForce_explicit


function [T_internal] = TrussCorrectedInternalForce_explicit(ielement,...
          T_internal,FEM,QUADRATURE,GEOM,GlobT_int,...
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
    % Get Correction force (Embedded Truss internal forces using host element
    % material properties)
    %--------------------------------------------------------------------------    
    %Assuming both elements are neo hookean materials
    properties_eh = properties_e; %eh is the same as e, except that nu and E are replaced bu nu and E of the host element
    properties_eh(1) = properties_h(1);
    properties_eh(2) = E_h;
    properties_eh(3) = nu_h;

    [TC,~,PLAST,~,~,Cauchy,~,CauchyTensor] = element_force_truss(...
      properties_eh,xelocal,x_e,FEM(2),PLAST,1,GEOM,DAMPING,1);

    
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

    %----------------------------------------------------------------------
    % Compute equivilant (internal) force vector of the host element.
    %----------------------------------------------------------------------
    if VolumeCorrect
        T_internal = T_internal + (T_e - T_C);
    else
        T_internal = T_internal + (T_e);
    end
 
        
end