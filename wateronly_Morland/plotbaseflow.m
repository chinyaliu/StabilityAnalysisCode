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
Ue(:,1) = ud_nd*erfc(-2*z/delta_nd/sqrt(pi));
Ue(:,2) = 4*ud_nd*exp(-4*z.^2/pi/delta_nd^2)/pi/delta_nd;
Ue(:,3) = Ue(:,2).*(-8*z/pi/delta_nd^2);
figure('position',[100,50,540,720]);
plot(Ue,z);

 %%
ud_nd = 2;
delta_nd = 0.234;
lambda_nd = 0.817;
h = 4*lambda_nd;
z = -linspace(0,h,100);
U = ud_nd*exp(2*z/delta_nd);
Ue = ud_nd*erfc(-2*z/delta_nd/sqrt(pi));
Uf = ud_nd*(exp(-0.25*pi*z.^2/delta_nd^2)+0.5*pi*z/delta_nd.*erfc(-0.5*sqrt(pi)*z/delta_nd));
figure('position',[100,50,540,720]);
hold on;
plot(U,z,'-k','Displayname','exp');
plot(Ue,z,'--b','Displayname','erf');
plot(Uf,z,':r','Displayname','I-erf');
hold off; box on; grid on;
xlabel('$z$','fontsize',36);
ylabel('$U(z)$','fontsize',36);
ylim([-0.8 0]);
legend('location','southeast');