function [oz] = MortoZ(lam_nd,om)
T = 75;           % Surface tension dyne/cm
rho_w = 1;        % Water density g/(cm)^3
rho_a = 0.0012;   % Air density g/(cm)^3
nu_w = 0.01;      % kinematic viscosity cm^2/s
g = 980;          % Acceleration due to gravity cm/(s^2)

lam_m = 2*pi*sqrt(T/rho_w/g);
c_m = (4*g*T/rho_w)^0.25;

lam = lam_nd*lam_m;
k = 2*pi/lam;
o0 = sqrt(g*k+(T*k^3)/rho_w);

o = om*c_m/lam_m;
oz = o/o0;
end