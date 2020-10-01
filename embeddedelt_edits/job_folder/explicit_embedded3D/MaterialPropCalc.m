%% Half Stiffness

mu=76.92E9; lam=115.4E9;
rho=7800;
K=lam+(2*mu/3);
C10=mu/2;
D1=2/K;


ce_m = sqrt((lam + 2*mu)/rho)
le = 0.2;
dt = 0.8*le/ce_m;

K=lam+(2*mu/3);
C10=mu/2;
D1=2/K;
% 

fprintf("Half Stiffness\n");
E=(mu*(3*lam+2*mu)/(lam+mu))
nu = lam/(2*(lam+mu))

E=E*0.5

lam2 = E*nu/((1+nu)*(1-2*nu))
mu2 = E/(2*(1+nu))
C10=mu2/2
D1=2/K

ce = sqrt((lam2 + 2*mu2)/rho)
le = 0.2;
dt = 0.8*le/ce

%% Half Density

mu=76.92E9; lam=115.4E9;
rho=7800;
K=lam+(2*mu/3);
C10=mu/2;
D1=2/K;

ce_m = sqrt((lam + 2*mu)/rho);
le = 0.2;
dt = 0.8*le/ce_m

K=lam+(2*mu/3);
C10=mu/2;
D1=2/K;
% 
fprintf("Half Density\n");
E=(mu*(3*lam+2*mu)/(lam+mu))
nu = lam/(2*(lam+mu))

rho=rho*0.5

lam2 = E*nu/((1+nu)*(1-2*nu))
mu2 = E/(2*(1+nu))
C10=mu2/2
D1=2/K

ce = sqrt((lam2 + 2*mu2)/rho)
le = 0.2;
dt = 0.8*le/ce

%%
E1=3E9; v1=0.4;

E2=900E9; v2=0.3;

lam1 = E1*v1/((1+v1)*(1-2*v1))
mu1 = E1/(2*(1+v1))
K1=lam1+(2*mu1/3);
C10=mu1/2;
D1=2/K1;

lam2 = E2*v2/((1+v2)*(1-2*v2))
mu2 = E2/(2*(1+v2))
K2=lam2+(2*mu2/3);
C10=mu2/2
D1=2/K2