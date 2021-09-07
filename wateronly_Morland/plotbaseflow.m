clear;
%%
ud_nd = 2;
delta_nd = 0.234;
lambda_nd = 0.817;
h = 4*lambda_nd;
z = -linspace(0,h,100);
U(:,1) = ud_nd*exp(2*z/delta_nd);
U(:,2) = (2/delta_nd)*U(:,1);
U(:,3) = (2/delta_nd)*U(:,2);
figure('position',[100,50,540,720]);
plot(U,z);

%%
hr = 100;
zr = -linspace(0,hr,100);
c1 = 0.9988; c2 = 0.8814;
Ur(:,1) = (1-c1*cosh(c2*zr).^(-2));
Ur(:,2) = 2*c1*c2*tanh(c2*zr).*(sech(c2*zr)).^2;
Ur(:,3) = 2*c1*c2^2*(sech(c2*zr)).^2.*((sech(c2*zr)).^2-2*(tanh(c2*zr)).^2);
figure('position',[100,50,540,720]);
plot(Ur,zr);

%%
ud_nd = 2;
delta_nd = 0.234;
lambda_nd = 0.817;
h = 4*lambda_nd;
z = -linspace(0,h,100).';
U(:,1) = ud_nd*erfc(-2*z/delta_nd/sqrt(pi));
U(:,2) = 4*ud_nd*exp(-4*z.^2/pi/delta_nd^2)/pi/delta_nd;
U(:,3) = U(:,2).*(-8*z/pi/delta_nd^2);
figure('position',[100,50,540,720]);
plot(U,z);
