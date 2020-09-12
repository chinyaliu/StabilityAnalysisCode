function U = baseflow_zhang(zeta,h,N)
z = 0.5*h*(zeta-1);
c1 = 0.9988; c2 = 0.8814;
U(:,:,1) = (1-c1*cosh(c2*z).^(-2))*ones(1,N+1);
U(:,:,2) = 2*c1*c2*tanh(c2*z).*(sech(c2*z)).^2*ones(1,N+1);
U(:,:,3) = 2*c1*c2^2*(sech(c2*z)).^2.*((sech(c2*z)).^2-2*(tanh(c2*z)).^2)*ones(1,N+1);
end