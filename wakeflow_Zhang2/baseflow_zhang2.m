function [U,z] = baseflow_zhang2(zeta,h,zL)
z = [0.5*zL*(zeta-1);0.5*(h-zL)*(zeta-1)-zL];
c1 = 0.9988; c2 = 0.8814;
U(:,1) = (1-c1*cosh(c2*z).^(-2));
U(:,2) = 2*c1*c2*tanh(c2*z).*(sech(c2*z)).^2;
U(:,3) = 2*c1*c2^2*(sech(c2*z)).^2.*((sech(c2*z)).^2-2*(tanh(c2*z)).^2);
end