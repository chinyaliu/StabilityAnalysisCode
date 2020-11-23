function U = baseflow(obj)
% obj.z = 0.5*obj.h*(obj.zeta-1+1i*obj.del*(1-obj.zeta.^2));
obj.z = 0.5*obj.h*(obj.zeta-1+1i*obj.del*cospi(obj.zeta/2));
c1 = 0.9988; c2 = 0.8814;
U(:,1) = (1-c1*cosh(c2*obj.z).^(-2));
U(:,2) = 2*c1*c2*tanh(c2*obj.z).*(sech(c2*obj.z)).^2;
U(:,3) = 2*c1*c2^2*(sech(c2*obj.z)).^2.*((sech(c2*obj.z)).^2-2*(tanh(c2*obj.z)).^2);
end