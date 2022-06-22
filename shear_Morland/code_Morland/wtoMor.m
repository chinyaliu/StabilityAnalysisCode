function [ud_nd,delta_nd,lam_nd,Re] = wtoMor(u,lam)
T = 75;           % Surface tension dyne/cm
rho_w = 1;        % Water density g/(cm)^3
rho_a = 0.0012;   % Air density g/(cm)^3
nu_w = 0.01;      % kinematic viscosity cm^2/s
g = 980;          % Acceleration due to gravity cm/(s^2)

U_0 = 0.5*u;       
delta = 2*U_0*rho_w*nu_w/rho_a/(u^2);
lam_m = 2*pi*sqrt(T/rho_w/g);
c_m = (4*g*T/rho_w)^0.25;

ud_nd = 0.5*u/c_m;
delta_nd = delta/lam_m;
lam_nd = lam/lam_m;
Re = c_m*lam_m/nu_w;
end