 %--------------------------------------------------------------------------
% Evaluates the Cauchy stress tensor (Mooney-Rivlin) for material type 99.
%--------------------------------------------------------------------------
function Cauchy = stress99(kinematics,properties,cons)
mu           = properties(2); %Assume mu1 = mu2 = mu
lambda       = properties(3); %Assume k = lambda
j            = kinematics.J;
F            = kinematics.F; %deformation gradient
B            = kinematics.b; 
F_bar        = j^(-1/3)*F;
B_bar        = F_bar*F_bar';
I_1          = trace(B_bar);
I_1_bar      = j^(-2/3)*I_1;
I_2          = 0.5*(I_1*I_1-trace(B*B));
I_2_bar      = j^(-4/3)*I_2;
sigma_vol    = (j-1)*lambda*cons.I;
sigma_iso    = j^(-1)*(-1/3*(mu*I_1_bar+2*mu*I_2_bar)*cons.I+(mu+mu*I_1_bar)*B_bar+mu*B_bar*B_bar);
Cauchy       = sigma_vol + sigma_iso;                      
end