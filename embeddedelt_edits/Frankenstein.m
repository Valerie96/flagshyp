%A bunch of pieces of flagshyp smashed together so I can actually figure
%out what all of the variables are. This may be a disaster
clear; clc; close all; 
% inputfile='explicit_embedded3D_new.dat';
% inputfile='explicit_embedded_4elt_new.dat';
inputfile='explicit_3D.dat';
% inputfile='seperate_embedded.dat';
% inputfile='seperate.dat';
inputfile='explicit_embedded_truss.dat';
basedir_fem='C:/Users/Valerie/Documents/GitHub/flagshyp/embeddedelt_edits/';
simtime = 0.01;
outputfreq=2;
DAMPING.b1 = 0.042; %Linear bulk viscosity damping
DAMPING.b2 = 0; %Quadratic bulk viscosity damping

ansmlv='y'; 
global explicit ;
explicit = 1;
global EmbedElt;
EmbedElt = 1;
tic
%% Input_data_and_initilaization.m

%--------------------------------------------------------------------------
% -Welcomes the user and determines whether the problem is being
%  restarted or a data file is to be read. 
% -Reads all necessary input data.
% -Initialises kinematic variables and internal variables.
% -Compute initial tangent matrix and equivalent force vector, excluding
%  pressure component.
%-------------------------------------------------------------------------- 
% function [PRO,FEM,GEOM,QUADRATURE,BC,MAT,LOAD,CON,CONSTANT,GLOBAL,...
%           PLAST,KINEMATICS] = input_data_and_initialisation(basedir_fem,ansmlv,inputfile)
%--------------------------------------------------------------------------
% Welcomes the user and determines whether the problem is being
% restarted or a data file is to be read.
%-------------------------------------------------------------------------- 
    PRO = welcome(basedir_fem,ansmlv,inputfile);
    fid = PRO.fid_input;
    %----------------------------------------------------------------------
    % Read input file, see textbook for user instructions.
    %----------------------------------------------------------------------
%     [FEM,GEOM,QUADRATURE,BC,MAT,LOAD,CON,PRO,GLOBAL] = ...
%      reading_input_file(PRO,fid);
%--------------------------------------------------------------------------
        % Read input data.
        %--------------------------------------------------------------------------
%         function [FEM,GEOM,QUADRATURE,BC,MAT,LOAD,CON,PRO,GLOBAL] = ...
%                   reading_input_file(PRO,fid)
        %--------------------------------------------------------------------------
        % Problem title.   
        %--------------------------------------------------------------------------
        PRO.title = strtrim(fgets(fid));
        %--------------------------------------------------------------------------
        % Element type.    
        %--------------------------------------------------------------------------
        [FEM,GEOM,QUADRATURE] = elinfo(fid);    
        %--------------------------------------------------------------------------
        % Obtain quadrature rules, isoparametric shape functions and their  
        % derivatives for the internal and boundary elements.
        %--------------------------------------------------------------------------
        for i = 1:FEM(1).n_elet_type
              QUADRATURE(i).element = element_quadrature_rules(FEM(i).mesh.element_type);
              QUADRATURE(i).boundary = edge_quadrature_rules(FEM(i).mesh.element_type);

              FEM(i).interpolation = [];
              FEM(i) = shape_functions_iso_derivs(QUADRATURE(i),FEM(i),GEOM.ndime);
        end
        %--------------------------------------------------------------------------
        % Read the number of mesh nodes, nodal coordinates and boundary conditions.  
        %--------------------------------------------------------------------------
        [GEOM,BC,FEM] = innodes(GEOM,fid,FEM);
        %--------------------------------------------------------------------------
        % Read the number of elements, element connectivity and material number.
        %--------------------------------------------------------------------------
        GEOM.total_n_elets = 0;
        for i = 1:FEM(1).n_elet_type
            [FEM(i),MATA(i)] = inelems(FEM(i) ,fid);
            GEOM.total_n_elets = GEOM.total_n_elets + FEM(i).mesh.nelem;
        end
        %--------------------------------------------------------------------------
        % Obtain fixed and free degree of freedom numbers (dofs).
        %--------------------------------------------------------------------------
        BC = find_fixed_free_dofs(GEOM,FEM(1),BC);
        %--------------------------------------------------------------------------
        % Read the number of materials and material properties.  
        %--------------------------------------------------------------------------
        for i = 1:FEM(1).n_elet_type
            MAT(i) = matprop(MATA(i),FEM(i),fid); %It didn't like returning such a 
                         %new looking MAT structure so MATA is a temp variable
        end
        clear MATA;
        %--------------------------------------------------------------------------
        % Read nodal point loads, prescribed displacements, surface pressure loads
        % and gravity (details in textbook).
        %--------------------------------------------------------------------------
        [LOAD,BC,FEM(1),GLOBAL] = inloads(GEOM,FEM(1),BC,fid);
        %--------------------------------------------------------------------------
        % Read control parameters.
        %--------------------------------------------------------------------------
        CON = incontr(BC,fid);
        fclose('all'); 
        

    %----------------------------------------------------------------------
    % Obtain entities which will be constant and only computed once.
    %----------------------------------------------------------------------
    CONSTANT = constant_entities(GEOM.ndime);
    %----------------------------------------------------------------------
    % Initialise load and increment parameters.
    %----------------------------------------------------------------------
    CON.xlamb = 0;
    CON.incrm = 0; 
    %----------------------------------------------------------------------
    % Initialises kinematic variables and compute initial tangent matrix 
    % and equivalent force vector, excluding pressure component.
    %----------------------------------------------------------------------
    
    [GEOM,LOAD,GLOBAL,PLAST,KINEMATICS] = ...
     initialisation(FEM,GEOM,QUADRATURE,MAT,LOAD,CONSTANT,CON,GLOBAL,BC);   
    
    %----------------------------------------------------------------------
    % |-/
    GLOBAL.external_load_effective = GLOBAL.external_load;
    %|-/
    %----------------------------------------------------------------------
    % Save into restart file.
    %----------------------------------------------------------------------
    cd(PRO.job_folder);
    save_restart_file(PRO,FEM,GEOM,QUADRATURE,BC,MAT,LOAD,CON,CONSTANT,...
                      GLOBAL,PLAST,KINEMATICS,'internal')    
    %output_vtk(PRO,CON,GEOM,FEM,BC,GLOBAL,MAT,PLAST,QUADRATURE.element,CONSTANT,KINEMATICS); 
    output_vtu(PRO,CON,GEOM,FEM,BC,GLOBAL,MAT,PLAST,QUADRATURE,CONSTANT,KINEMATICS);


%% ExplicitDynamics_algorithm
%--------------------------------------------------------------------------
% Explicit central diff. algorithm 
%--------------------------------------------------------------------------
% function ExplicitDynamics_algorithm(PRO,FEM,GEOM,QUADRATURE,BC,MAT,LOAD,...
%                                   CON,CONSTANT,GLOBAL,PLAST,KINEMATICS)
 d1=digits(64);
%step 1 - iniitalization
%       - this is done in the intialisation.m file, line 68
CON.xlamb = 0;
CON.incrm = 0; 

%step 2 - getForce
[GLOBAL,updated_PLAST,GEOM.Jn_1,GEOM.VolRate] = getForce_explicit(CON.xlamb,...
          GEOM,MAT,FEM,GLOBAL,CONSTANT,QUADRATURE,PLAST,KINEMATICS,BC,DAMPING,1);      
     
%step 3 - compute accelerations.
GLOBAL.accelerations = inv(GLOBAL.M)*(GLOBAL.external_load_effective - GLOBAL.T_int);

velocities_half = zeros(FEM(1).mesh.n_dofs,1);
disp_n = zeros(FEM(1).mesh.n_dofs,1);
disp_prev = zeros(FEM(1).mesh.n_dofs,1);

%step 4 - time update/iterations
Time = 0; 
tMax = simtime; % in seconds
prefactor = 0.8;
dt= prefactor * CalculateTimeStep(FEM(1),GEOM,MAT(1),DAMPING); % in seconds
plot_counter  = 0;
time_step_counter = 0;
plot_counter = 0;
nPlotSteps = 20 ;
nSteps = round(tMax/dt);
nsteps_plot = round(nSteps/nPlotSteps);

% start explicit loop
while(Time<=tMax)
    t_n=Time;
    t_np1 = Time + dt;
    Time = t_np1; % update the time by adding full time step
    dt_nphalf = dt; % equ 6.2.1
    t_nphalf = 0.5 *(t_np1 + t_n); %equ 6.2.1
    
% step 5 - update velocities
    velocities_half = GLOBAL.velocities + (t_nphalf - t_n) * GLOBAL.accelerations;
 
   % step 7 Update nodal displacments 
    % store old displacements for energy computation
    disp_prev = disp_n;
    % update nodal displacements 
    disp_n = disp_n + dt_nphalf *velocities_half;
    
%----------------------------------------------------------------
% Update stored coodinates.
  
  displ = disp_n-disp_prev; 
  GEOM.x = update_geometry_explicit(GEOM.x,GEOM.x0,1,disp_n(BC.freedof),BC.freedof);
  dx = GEOM.x - GEOM.x0;
%----------------------------------------------------------------   

% step 6 - enforce displacement BCs 
%   %--------------------------------------------------------------------
%   % Update nodal forces (excluding pressure) and gravity. 
%   %--------------------------------------------------------------------
%   [GLOBAL.Residual,GLOBAL.external_load] = ...
%    external_force_update(GLOBAL.nominal_external_load,...
%    GLOBAL.Residual,GLOBAL.external_load,CON.dlamb);
%   %--------------------------------------------------------------------
%   % Update nodal forces and stiffness matrix due to external pressure 
%   % boundary face (line) contributions. 
%   %--------------------------------------------------------------------      
%   if LOAD.n_pressure_loads      
%      GLOBAL = pressure_load_and_stiffness_assembly(GEOM,MAT,FEM,...
%               GLOBAL,LOAD,QUADRATURE.boundary,CON.dlamb);    
%   end
  %--------------------------------------------------------------------
  % For the case of prescribed geometry update coodinates.
  % -Recompute equivalent nodal forces and assembles residual force, 
  % excluding pressure contributions.
  % -Recompute and assembles tangent stiffness matrix components, 
  % excluding pressure contributions.
  %--------------------------------------------------------------------
  if  BC.n_prescribed_displacements > 0
      [GEOM.x ,velocities_half]  = update_prescribed_displacements_explicit(BC.dofprescribed,...
               GEOM.x0,GEOM.x,velocities_half,BC.presc_displacement,t_n,tMax); 
       disp_n(BC.fixdof) = GEOM.x(BC.fixdof) - GEOM.x0(BC.fixdof);  
%      |-/
%      Update coodinates of embedded nodes (if there are any)  
%           if ~isempty(FEM(1).mesh.embedded) || ~isempty(FEM(2).mesh.embedded)
          if EmbedElt == 1
              [GEOM.x,velocities_half, GLOBAL.accelerations ] = update_embedded_displacements_explicit(BC.tiedof, BC.tienodes,...
                    FEM,GEOM, velocities_half, GLOBAL.accelerations); 
              disp_n(BC.tiedof) = GEOM.x(BC.tiedof) - GEOM.x0(BC.tiedof);
          end

         dx = GEOM.x - GEOM.x0;
%       [GLOBAL,updated_PLAST] = residual_and_stiffness_assembly(CON.xlamb,...
%        GEOM,MAT,FEM,GLOBAL,CONSTANT,QUADRATURE.element,PLAST,KINEMATICS);
      %----------------------------------------------------------------
      % Update nodal forces due to pressure. 
      %----------------------------------------------------------------
%       if LOAD.n_pressure_loads
%           GLOBAL = pressure_load_and_stiffness_assembly(GEOM,MAT,FEM,...
%                    GLOBAL,LOAD,QUADRATURE.boundary,CON.xlamb);
%       end
  end
  
% %----------------------------------------------------------------   
  
  % save internal force, to be used in energy computation
  fi_prev = GLOBAL.T_int;
  fe_prev = GLOBAL.external_load_effective;
  
%step 8 - getForce
  [GLOBAL,updated_PLAST,GEOM.Jn_1,GEOM.VolRate] = getForce_explicit(CON.xlamb,...
          GEOM,MAT,FEM,GLOBAL,CONSTANT,QUADRATURE,PLAST,KINEMATICS,BC,DAMPING,dt);
      
  % updated stable time increment based on current deformation     
  dt_old=dt;
  dt = prefactor * CalculateTimeStep(FEM(1),GEOM,MAT(1),DAMPING);
  
%step 9 - compute accelerations.       
  AccOld = GLOBAL.accelerations;
  GLOBAL.accelerations = inv(GLOBAL.M)*(GLOBAL.external_load_effective - GLOBAL.T_int);
  
% step 10 second partial update of nodal velocities
  VelOld = GLOBAL.velocities;
  GLOBAL.velocities = velocities_half + (t_np1 - t_nphalf) * GLOBAL.accelerations;
  
%      |-/
%      Update v/a of embedded nodes (if there are any)  
%           if ~isempty(FEM(1).mesh.embedded) || ~isempty(FEM(2).mesh.embedded)
          if EmbedElt == 1
              [GEOM.x,GLOBAL.velocities, GLOBAL.accelerations ] = update_embedded_displacements_explicit(BC.tiedof, BC.tienodes,...
                    FEM,GEOM, GLOBAL.velocities, GLOBAL.accelerations); 
               disp_n(BC.tiedof) = GEOM.x(BC.tiedof) - GEOM.x0(BC.tiedof);
          end
  
%--------------------------------------------------------------
  
% step 11 check energy
  [energy_value, max_energy] = check_energy_explicit(PRO,FEM,CON,BC, ...
      GLOBAL,disp_n, disp_prev,GLOBAL.T_int,fi_prev,...
      GLOBAL.external_load_effective,fe_prev,t_n);
  
%Plot every # steps
  if( mod(time_step_counter,outputfreq) == 0 )
      plot_counter = plot_counter +1;
      
      output(PRO,CON,GEOM,FEM,BC,GLOBAL,MAT,PLAST,QUADRATURE,CONSTANT,KINEMATICS);
      output_vtu(PRO,CON,GEOM,FEM,BC,GLOBAL,MAT,PLAST,QUADRATURE,CONSTANT,KINEMATICS);

%       PLAST = save_output(updated_PLAST,PRO,FEM,GEOM,QUADRATURE,BC,...
%                 MAT,LOAD,CON,CONSTANT,GLOBAL,PLAST,KINEMATICS);  
  end
  
  time_step_counter = time_step_counter + 1;  
  
  % this is set just to get the removal of old vtu files in output.m
  % correct
  CON.incrm =  CON.incrm + 1;
  
  if( mod(time_step_counter,outputfreq) == 0 )
    disp(['step = ',sprintf('%d', time_step_counter-1),'     time = ',... 
  sprintf('%.2e', t_n), ' sec.     dt = ', sprintf('%.2e', dt_old) ,...
  ' sec.'])
  end

formt = [repmat('% -1.4E ',1,3) '\n'];

end % end on while loop
                   

fprintf(' Normal end of PROGRAM flagshyp. \n');


toc

delete *RESTART* 
digits(d1);

