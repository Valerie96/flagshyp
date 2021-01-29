%--------------------------------------------------------------------------
% Computes and assemble residual force vector and global tangent stiffness
% matrix except surface (line) element pressure contributions.
%--------------------------------------------------------------------------
function [GLOBAL,updated_PLAST,Jn_1_vec,VolRate_vec] = getForce_explicit(xlamb,...
          GEOM,MAT,FEM,GLOBAL,CONSTANT,QUADRATURE,PLAST,KINEMATICS,BC,DAMPING,dt)    
%--------------------------------------------------------------------------
% Initialisation of the updated value of the internal variables.
%--------------------------------------------------------------------------
updated_PLAST = PLAST;
%--------------------------------------------------------------------------
% Initialises the total external load vector except pressure contributions.
%--------------------------------------------------------------------------
GLOBAL.external_load = xlamb*GLOBAL.nominal_external_load;

Jn_1_vec=ones(FEM(1).mesh.nelem,1);
VolRate_vec=zeros(FEM(1).mesh.nelem,QUADRATURE(1).element.ngauss);



Step_globalT_int = zeros(size(GLOBAL.T_int,1),1);

counter          = 1;                                        

% Main element loop.
%--------------------------------------------------------------------------
% for ielement=1:FEM.mesh.nelem
for ii=1:FEM(1).mesh.nelem
    ielement=ii;
    %----------------------------------------------------------------------
    % GATHER Temporary variables associated with a particular element.
    %----------------------------------------------------------------------
    global_nodes    = FEM(1).mesh.connectivity(:,ielement);   
    material_number = MAT(1).matno(ielement);     
    matyp           = MAT(1).matyp(material_number);        
    properties      = MAT(1).props(:,material_number); 
    xlocal          = GEOM.x(:,global_nodes);                     
    x0local         = GEOM.x0(:,global_nodes);                       
    Ve              = GEOM.Ve(ielement);      
%     v_e             = GLOBAL.velocities(global_nodes);
%     a_e             = GLOBAL.accelerations(global_nodes);


    %----------------------------------------------------------------------
    % Select internal variables within the element (plasticity).
    %----------------------------------------------------------------------
    PLAST_element = selecting_internal_variables_element(PLAST,matyp,ielement);    
    %----------------------------------------------------------------------
    % Compute internal force and stiffness matrix for an element.
    %----------------------------------------------------------------------    
    switch FEM(1).mesh.element_type
      case 'truss2'
       [T_internal,indexi,indexj,global_stiffness,counter,PLAST_element] = ...
        element_force_and_stiffness_truss(properties,xlocal,x0local,...
        global_nodes,FEM(1),PLAST_element,counter,indexi,indexj,...
        global_stiffness,GEOM);
      otherwise
       [T_internal,counter,PLAST_element,Jn_1,VolRate] = ...
        InternalForce_explicit(ielement,FEM,xlocal,x0local,global_nodes,...
        Ve,QUADRATURE,properties,CONSTANT,GEOM,matyp,PLAST_element,...
        counter,KINEMATICS,MAT,DAMPING,dt);
        
        Jn_1_vec(ielement) = Jn_1;
        VolRate_vec(ielement) = VolRate;
    end
    formt = [repmat('% -1.4E ',1,3) '\n'];
%     fprintf("\nElementForce:\n");
%     fprintf("Element %u\n",ielement);
%     for i=1:3:24
%         fprintf(formt, T_internal(i:i+2));
%     end
%     fprintf('\n');
    %----------------------------------------------------------------------
    % Assemble element contribution into global internal force vector.   
    %----------------------------------------------------------------------
    Step_globalT_int = force_vectors_assembly(T_internal,global_nodes,...
                   Step_globalT_int,FEM(1).mesh.dof_nodes);
               
%    fprintf("\nStep_globalT_int:\n");
%     for i=1:3:48
%         fprintf(formt, Step_globalT_int(i:i+2));
%     end
%     fprintf('\n');
    %----------------------------------------------------------------------
    % Storage of updated value of the internal variables. 
    %----------------------------------------------------------------------    
    updated_PLAST = plasticity_storage(PLAST_element,updated_PLAST,matyp,...
                                       ielement);    

end
    
    GLOBAL.T_int = Step_globalT_int;
%     Step_globalT_int
%        fprintf("\nGlobalForce:\n");
%     for i=1:3:48
%         fprintf(formt, GLOBAL.T_int(i:i+2));
%     end
%     fprintf('\n');
%--------------------------------------------------------------------------
% Global tangent stiffness matrix sparse assembly except pressure contributions. 
%--------------------------------------------------------------------------
% GLOBAL.K = sparse(indexi,indexj,global_stiffness);               
%--------------------------------------------------------------------------
% Compute global residual force vector except pressure contributions.
%--------------------------------------------------------------------------
% GLOBAL.Residual = GLOBAL.T_int - GLOBAL.external_load;

% |-/

GLOBAL.external_load_effective(BC.fixdof) = GLOBAL.T_int(BC.fixdof);
% GLOBAL.external_load_effective(BC.tiedof) = GLOBAL.T_int(BC.tiedof);

% algorithm in box 6.1 gives f_n = f_ext - f_int
GLOBAL.Residual = GLOBAL.external_load_effective - GLOBAL.T_int;
GLOBAL.Reactions(BC.fixdof) =  GLOBAL.Residual(BC.fixdof) + GLOBAL.external_load_effective(BC.fixdof);
% GLOBAL.Reactions(BC.tiedof) =  GLOBAL.Residual(BC.tiedof) + GLOBAL.external_load_effective(BC.tiedof);
% 
% ffid = fopen('GlobalForce.txt','a+');
%     fprintf(ffid,"Global Internal Force:\n"); 
%     for i=1:3:60
%         fprintf(ffid,formt, GLOBAL.T_int(i:i+2));
%     end
%     fprintf(ffid,'\n');
%     fprintf(ffid,"Global Reaction Force:\n"); 
%     for i=1:3:60
%         fprintf(ffid,formt, GLOBAL.Reactions(i:i+2));
%     end
%     fprintf(ffid,'\n');
%         fprintf(ffid,"Global Residual Force:\n"); 
%     for i=1:3:60
%         fprintf(ffid,formt, GLOBAL.Residual(i:i+2));
%     end
%     fclose(ffid);

    
end


 


 
