%--------------------------------------------------------------------------
% Evaluates the constitutive tensor (Mooney-Rivlin) for material type 99.
%--------------------------------------------------------------------------
function c   = ctens99(kinematics,properties,cons)
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
Id           = cons.I;
c_vol        = lambda*(j*(2*j-1)*kron(Id,Id))-j*(j-1)*(kron(Id,Id)+kron(Id,Id));
c_iso        = 2*mu*(kron(B_bar,B_bar)-1/2*(kron(B_bar,B_bar)*2))-2/3*(mu+2*mu*I_1_bar)*(kron(B_bar,Id)+kron(B_bar,Id))+...
               4/3*mu*(kron(B_bar^2,Id)*2) + 2/9*(mu*I_1_bar+4*mu*I_2_bar)*kron(Id,Id)+...
               1/3*(mu*I_1_bar + 2*mu*I_2_bar)*kron(Id,Id)*2;
c            = c_vol + C_iso;
end