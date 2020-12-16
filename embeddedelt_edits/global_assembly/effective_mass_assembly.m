%--------------------------------------------------------------------------
% Computes and assemble lumped mass matrix
%--------------------------------------------------------------------------
function [GLOBAL] = effective_mass_assembly(GEOM,MAT,FEM,GLOBAL,QUADRATURE)   

MassCorrect = false;

%number of dimensions
ndims = GEOM.ndime; 
massSize = ndims*GEOM.npoin;
LumpedMass = zeros(massSize,massSize);
n_nodes_elem = FEM.mesh.n_nodes_elem;

M_total = 0;
%Loop over all (host) elements
for ii=1:length(FEM.mesh.host)
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
   

    
    % assign density
    rho = properties(1);
    
    
    % quadrature locations
    QUADRATURE.Chi;
    % quadrature weights
    QUADRATURE.W;
   

    %b: calculate element mass
    me = rho * Ve;
    
    %c: Loop over embedded elements
    mf = 0; mc =0;
    
        for jj=1:length(FEM.mesh.embedded)
           jelement = FEM.mesh.embedded(jj); 
           if GEOM.embedded.ElementHost(jelement,2) == ielement
               %i: get element vol and other properties
               %later for "truss" elements calculate volume
                    global_nodes_f    = FEM.mesh.connectivity(:,jelement);   
                    material_numberf = MAT.matno(jelement);     
                    matyp_f           = MAT.matyp(material_numberf);        
                    properties_f      = MAT.props(:,material_numberf); % needs density included.
                    xflocal          = GEOM.x(:,global_nodes_f);                     
                    xf0local         = GEOM.x0(:,global_nodes_f); 
                    Vf              = GEOM.Ve(jelement);
                    rho_f            = properties_f(1);
               
                %ii: Calculate embedded element mass
                    mf = rho_f * Vf;
                    
                %iii: Calculate correction mass 
                    mc = rho * Vf;
                    
           %***Until I work out removing the embedded nodes as dof, I need to
           %     add them to the mass matrix
                    %e: divide mass equally among nodes
                    Meff = (mf/n_nodes_elem);

                    %f: form diagonal mass matrix
                    Me = eye(n_nodes_elem* ndims,n_nodes_elem* ndims) * Meff;

        
                    %g: scatter mass matrix to global 
                    global_dof = FEM.mesh.dof_nodes(:,global_nodes_f);
                    LumpedMass(global_dof,global_dof) = LumpedMass(global_dof,global_dof) + Me;

                    
           end
            
        end
    
        %d: calculate effective mass
        if MassCorrect
            meff = me + mf - mc;
        else
            meff = me;
        end
        M_total = M_total + meff + mf;
        
        %e: divide mass equally among nodes
        Meff = (meff/n_nodes_elem);
        
        %f: form diagonal mass matrix
        Me = eye(n_nodes_elem* ndims,n_nodes_elem* ndims) * Meff;

        
        %g: transform and scatter scatter mass matrix to global 
        global_dof = FEM.mesh.dof_nodes(:,global_nodes);
        LumpedMass(global_dof,global_dof) = LumpedMass(global_dof,global_dof) + Me;
        
        
end % loop on elements

GLOBAL.M=LumpedMass;
fprintf('Total mesh mass is: %15.5f \n', M_total);

M_total = 0 ;
for i = 1:massSize
    M_total = M_total + LumpedMass(i,i);
end
fprintf('Sum Mass: %15.5f \n', M_total/3);

fprintf('\n');

end

 


 
