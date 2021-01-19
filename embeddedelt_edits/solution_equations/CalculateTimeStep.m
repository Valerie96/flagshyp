%--------------------------------------------------------------------------
% Check for equilibrium convergence.
%--------------------------------------------------------------------------
function [dt] = CalculateTimeStep(FEM,GEOM,MAT,DAMPING)

dtMax = 1e20;
% Main element loop.
%--------------------------------------------------------------------------
for ielement=1:FEM.mesh.nelem
    mu = MAT.props(2);
    lambda = MAT.props(3);
    rho = MAT.props(1);
    %|-/
    b1 = DAMPING.b1;
    b2 = DAMPING.b2;
    eps_dot = GEOM.VolRate(ielement);

    %Longitudinal Bulk wave speed
    ce = sqrt((lambda + 2*mu)/rho);
    
    [le, max_le] = calc_element_size(FEM,GEOM,ielement);
    dt_ielt = le/ce;
    
    %Add effect of damping 
    zeta = b1-(b2^2)*dt_ielt*min(0, eps_dot);
    dt_ielt = dt_ielt*(sqrt(1+zeta^2)-zeta);
    
    if(dt_ielt < dtMax)
        dt = dt_ielt;
        dtMax = dt;
    end
    
end

        lmax = le;
        zmax = b1-(b2^2)*dt_ielt*min(0, eps_dot);
        ee=min(0, eps_dot);

% Ffid = fopen('Damping.txt','a+');
% formt = [repmat('%1.4d ',1,3) '\n'];
% fprintf(Ffid, 'z=%d   le=%d  e=%d  dt=%d \n',zmax,lmax,ee,dt);
% fclose(Ffid);

end



