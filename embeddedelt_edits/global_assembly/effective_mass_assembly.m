%--------------------------------------------------------------------------
% Computes and assemble lumped mass matrix. Mass Correct gives the option
% to remove redundant mass from the system when using embedded elements.
%--------------------------------------------------------------------------
function [GLOBAL] = effective_mass_assembly(GEOM,MAT,FEM,GLOBAL,QUADRATURE)   


global VolumeCorrect;

%number of dimensions
ndims = GEOM.ndime; 
massSize = ndims*GEOM.npoin;
LumpedMass_eff = zeros(massSize,massSize); %Embedded element mass included in host element mass
LumpedMass_act = zeros(massSize,massSize); %Host elements without embedded element mass


M_total = 0;
%Loop over all (host) elements
for ii=1:length(FEM(1).mesh.nelem)
    ielement=ii;
    %----------------------------------------------------------------------
    % Temporary variables associated with a particular element.
    %----------------------------------------------------------------------
    n_nodes_elem    = FEM(1).mesh.n_nodes_elem;
    global_nodes    = FEM(1).mesh.connectivity(:,ielement);   
    material_number = MAT(1).matno(ielement);     
    matyp           = MAT(1).matyp(material_number);        
    properties      = MAT(1).props(:,material_number); % needs density included.
    xlocal          = GEOM.x(:,global_nodes);                     
    x0local         = GEOM.x0(:,global_nodes); 
    Ve              = GEOM.Ve(ielement,1);
   

    
    % assign density
    rho = properties(1);
    
    
    % quadrature locations
    QUADRATURE(1).element.Chi;
    % quadrature weights
    QUADRATURE(1).element.W;
   

    %b: calculate element mass
    me = rho * Ve;
    
    %c: Loop over embedded elements
    mf_tot = 0; mc =0;
    
        for jj=1:GEOM.embedded.HostTotals(ielement,2)
           jelement = jj; 
           if GEOM.embedded.ElementHost(jelement,1) == ielement
               %i: get element vol and other properties
               %later for "truss" elements calculate volume
                    n_nodes_elem_f    = FEM(2).mesh.n_nodes_elem;
                    global_nodes_f    = FEM(2).mesh.connectivity(:,jelement);   
                    material_numberf = MAT(2).matno(jelement);     
                    matyp_f           = MAT(2).matyp(material_numberf);        
                    properties_f      = MAT(2).props(:,material_numberf); % needs density included.
                    xflocal          = GEOM.x(:,global_nodes_f);                     
                    xf0local         = GEOM.x0(:,global_nodes_f); 
                    Vf               = GEOM.Ve(jelement,2);
                    rho_f            = properties_f(1);
               
                %ii: Calculate embedded element mass
                    mf = rho_f * Vf;
                    mf_tot = mf_tot + mf;
                %iii: Calculate correction mass 
                    mc = mc + rho * Vf;
                    
           %***Until I work out removing the embedded nodes as dof, I need to
           %     add them to the mass matrix
                    %e: divide mass equally among embedded nodes
                    Meff = (mf/n_nodes_elem_f);

                    %f: form diagonal mass matrix
                    Identity = eye(n_nodes_elem_f* ndims,n_nodes_elem_f* ndims);
                    Me = Identity * Meff;
                    
                    %g: scatter mass matrix to global 
                    global_dof = FEM(2).mesh.dof_nodes(:,global_nodes_f);
                    LumpedMass_eff(global_dof,global_dof) = LumpedMass_eff(global_dof,global_dof) + Me;
                    LumpedMass_act(global_dof,global_dof) = LumpedMass_act(global_dof,global_dof) + Me;

                    
           end
            
        end
    
        %d: calculate effective mass
        if VolumeCorrect
            meff = me + mf_tot - mc;
            mact = me - mc;
        else
            meff = me + mf_tot;
            mact = me;
        end
        M_total = M_total + mact + mf_tot;
        
        %e: divide mass equally among nodes
        Meff = (meff/n_nodes_elem);
        Mact = (mact/n_nodes_elem);
        
        %f: form diagonal mass matrix
        Me_eff = eye(n_nodes_elem* ndims,n_nodes_elem* ndims) * Meff;
        Me_act = eye(n_nodes_elem* ndims,n_nodes_elem* ndims) * Mact;

        
        %g: transform and scatter scatter mass matrix to global 
        global_dof = FEM(1).mesh.dof_nodes(:,global_nodes);
        LumpedMass_eff(global_dof,global_dof) = LumpedMass_eff(global_dof,global_dof) + Me_eff;
        LumpedMass_act(global_dof,global_dof) = LumpedMass_act(global_dof,global_dof) + Me_act;
        
end % loop on elements

GLOBAL.M    = LumpedMass_eff;
GLOBAL.M_KE = LumpedMass_act;
fprintf('Total mesh mass is: %15.5f \n', M_total);

M_total = 0 ;
for i = 1:massSize
    M_total = M_total + LumpedMass_act(i,i);
end
fprintf('Sum Mass: %15.5f \n', M_total/3);

fprintf('\n');

end

 


 
