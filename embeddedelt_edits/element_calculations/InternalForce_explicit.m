%--------------------------------------------------------------------------
% Computes the element vector of global internal forces and the tangent
% stiffness matrix. 
%--------------------------------------------------------------------------
function [T_internal,counter,PLAST_element,geomJn_1,VolRate] = ...
          InternalForce_explicit(ielement,FEM,xlocal,x0local,...
          element_connectivity,Ve,QUADRATURE,properties,CONSTANT,GEOM,...
          matyp,PLAST,counter,KINEMATICS,MAT,GlobT_int,DAMPING,dt)
      
%define explicit as global variable in order to use it in function      
% it is assigned value in Flagshyp.m
global explicit
dim=GEOM.ndime;


% step 2.II       
T_internal = zeros(FEM(1).mesh.n_dofs_elem,1);
%--------------------------------------------------------------------------
% Computes initial and current gradients of shape functions and various 
% strain measures at all the Gauss points of the element.
%--------------------------------------------------------------------------
KINEMATICS(1) = gradients(xlocal,x0local,FEM(1).interpolation.element.DN_chi,...
             QUADRATURE(1).element,KINEMATICS(1));

%|-/
Jn_1=GEOM.Jn_1(ielement);
J=KINEMATICS(1).J(1);
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
for igauss=1:QUADRATURE(1).element.ngauss
    %----------------------------------------------------------------------
    % Extract kinematics at the particular Gauss point.
    %----------------------------------------------------------------------
    kinematics_gauss = kinematics_gauss_point(KINEMATICS(1),igauss);     
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
    le=calc_element_size(FEM(1),GEOM(1),ielement);
    rho=properties(1); mu=properties(2); lambda=properties(3);
    Cd=sqrt((lambda + 2*mu)/rho);
    
    p1 = rho*b1*le*Cd*eps_dot*CONSTANT.I;
    p2 = rho*(b2*le)^2*abs(eps_dot)*min(0,eps_dot)*CONSTANT.I;
    %
    %|-/
        
    %----------------------------------------------------------------------
    % Compute numerical integration multipliers.
    %----------------------------------------------------------------------
    JW = kinematics_gauss.Jx_chi*QUADRATURE(1).element.W(igauss)*...
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

%|-/
% Embedded Elt, Internal force modification, if this element has any
% embedded elements

if GEOM.embedded.HostTotals(ielement,2) > 0
    search = 0;
    for eelt = 1:GEOM.npoin
        if GEOM.embedded.ElementHost(eelt) == ielement
%             fprintf("Host: %u    Guest: %u\n", ielement, eelt);

testfid = fopen('Force.txt','a');
formt = [repmat('% -1.4E ',1,3)];
    fprintf(testfid,'\n');
    for i=1:3:24
        fprintf(testfid,formt, T_internal(i:3));
    end            
    
              T_internal = DistributeCorrectInternalForce_explicit(ielement,...
                           T_internal,FEM, xlocal,x0local,...
                           QUADRATURE,CONSTANT,GEOM,GlobT_int,...
                           PLAST,KINEMATICS,MAT,DAMPING,eelt);


% 
% 
%               T_internal = CorrectInternalForce_explicit(ielement,...
%                         T_internal,FEM,xlocal,x0local,QUADRATURE,CONSTANT,GEOM,...
%                             PLAST,KINEMATICS,MAT,DAMPING,eelt);
% 
    fprintf(testfid,'\n');
    for i=1:3:24
        fprintf(testfid,formt, T_internal(i:3));
    end  
    fprintf(testfid,'\n');
    fclose(testfid);
                        
                        
                        
              search = search + 1;
              if search >= GEOM.embedded.HostTotals(ielement,2)
                  break;
              end
        end
    end
end
%--------------------------------------------------------------------------
% Store internal variables.
%--------------------------------------------------------------------------
PLAST_element = PLAST;
end

