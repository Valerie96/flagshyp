 %--------------------------------------------------------------------------
% Computes the element vector of global internal forces and the tangent
% stiffness matrix. 
%--------------------------------------------------------------------------
function [T_internal,counter,PLAST_element,geomJn_1,VolRate] = ...
          InternalForce_HostElt_explicit(ielement,FEM,xlocal,x0local,...
          element_connectivity,Ve,QUADRATURE,properties,CONSTANT,GEOM,...
          matyp,PLAST,counter,KINEMATICS,MAT,DAMPING,dt)
      
%define explicit as global variable in order to use it in function      
% it is assigned value in Flagshyp.m
global explicit
dim=GEOM.ndime;

%First, run Internal Force as usual for the host element

% step 2.II       
T_internal = zeros(FEM.mesh.n_dofs_elem,1);
%--------------------------------------------------------------------------
% Computes initial and current gradients of shape functions and various 
% strain measures at all the Gauss points of the element.
%--------------------------------------------------------------------------
KINEMATICS = gradients(xlocal,x0local,FEM.interpolation.element.DN_chi,...
             QUADRATURE,KINEMATICS);

%|-/
Jn_1=GEOM.Jn_1(ielement);
J=KINEMATICS.J(1);
eps_dot = (J-Jn_1)/dt;
b1=DAMPING.b1;
b2=DAMPING.b2;

%--------------------------------------------------------------------------
% Computes element mean dilatation kinematics, pressure and bulk modulus. 
%--------------------------------------------------------------------------
switch matyp
     case {5,7,17}
          [pressure,kappa_bar,DN_x_mean,ve] = ...
           mean_dilatation_pressure(FEM,dim,matyp,properties,Ve,...
                                    QUADRATURE,KINEMATICS);
     otherwise
          pressure = 0;
end
%--------------------------------------------------------------------------
% Gauss quadrature integration loop.
%--------------------------------------------------------------------------
for igauss=1:QUADRATURE.ngauss
    %----------------------------------------------------------------------
    % Extract kinematics at the particular Gauss point.
    %----------------------------------------------------------------------
    kinematics_gauss = kinematics_gauss_point(KINEMATICS,igauss);     
    %----------------------------------------------------------------------
    % Obtain stresses (for incompressible or nearly incompressible, 
    % only deviatoric component) and internal variables in plasticity.
    %----------------------------------------------------------------------    
    [Cauchy,PLAST,...
     plast_gauss] = Cauchy_type_selection(kinematics_gauss,properties,...
                                          CONSTANT,dim,matyp,PLAST,igauss);
    %----------------------------------------------------------------------
    % Obtain elasticity tensor (for incompressible or nearly incompressible, 
    % only deviatoric component).
    %----------------------------------------------------------------------   
    if(explicit==0)
        c = elasticity_modulus_selection(kinematics_gauss,properties,CONSTANT,...
                                      dim,matyp,PLAST,plast_gauss,igauss);
    else
        c = 0;
    end
    %----------------------------------------------------------------------
    % Add pressure contribution to stresses and elasticity tensor.
    %----------------------------------------------------------------------    
    [Cauchy,c] = mean_dilatation_pressure_addition(Cauchy,c,CONSTANT,pressure,matyp);    
    %----------------------------------------------------------------------
    %|-/ 
    % Calculate bulk viscosity damping
    le=calc_element_size(FEM,GEOM,ielement);
    rho=properties(1); mu=properties(2); lambda=properties(3);
    Cd=sqrt((lambda + 2*mu)/rho);
    
    p1 = rho*b1*le*Cd*eps_dot*CONSTANT.I;
    p2 = rho*(b2*le)^2*abs(eps_dot)*min(0,eps_dot)*CONSTANT.I;
    %
    %|-/
        
    %----------------------------------------------------------------------
    % Compute numerical integration multipliers.
    %----------------------------------------------------------------------
    JW = kinematics_gauss.Jx_chi*QUADRATURE.W(igauss)*...
         thickness_plane_stress(properties,kinematics_gauss.J,matyp);
    %----------------------------------------------------------------------
    % Compute equivalent (internal) force vector.
    %----------------------------------------------------------------------
%     T = Cauchy*kinematics_gauss.DN_x;
    T = (Cauchy+p1+p2)*kinematics_gauss.DN_x;
    T_internal = T_internal + T(:)*JW;
    
end

% |-/
    %Update previous Jacobian and element strain rate
    %Assuming that J and eps_dot are the same for all the element Gauss Pts
    %and I have now confirmed this
    geomJn_1=J;
    VolRate = eps_dot;
% |-/

%--------------------------------------------------------------------------
% Compute conttribution (and extract relevant information for subsequent
% assembly) of the mean dilatation term (Kk) of the stiffness matrix.
%--------------------------------------------------------------------------
switch matyp
    case {5,7,17}         
         [indexi,indexj,global_stiffness,...
          counter] = mean_dilatation_volumetric_matrix(FEM,dim,...
          element_connectivity,DN_x_mean,counter,indexi,indexj,...
          global_stiffness,kappa_bar,ve);
end 

%--------------------------------------------------------------------------
% Store internal variables.
%--------------------------------------------------------------------------
PLAST_element = PLAST;





end

