close all; clear all;% clc
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [2,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
N = 1001;
k = 1.6;
Re = inf;
Fr2 = 2.25;
% h = @(k) 2*pi/k;
h = @(k) 6;
%% Run solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
[o, an, cA, errGEP, db] = wZhang_solver(N,k,h(k),Re,Fr2,method,alg,do_balancing);
%% Find critical height
if (strcmpi(method(1),'d4') && strcmpi(method(2),'schimd') && strcmpi(method(3),'d4'))
    M = N-2;
    [zeta,D] = Dcheb(N,1,'d4');
else
    M = N;
    [zeta,D] = Dcheb(N,1,'n');
end
z = 0.5*h(k)*(zeta-1);
c1 = 0.9988; c2 = 0.8814;
U = (1-c1*cosh(c2*z).^(-2));
% [minpt,minpos] = min(abs(U-real(o)/k));
%% Plot streamfunction phi
switch lower(method(1))
    case {'d2', 'd4'}
        phi = D(:,:,1)*an(1:N+1);
        up = D(:,:,2)*an(1:N+1);
        wp = -1i*k*phi;
    case 'uw'
        up = an(1:N+1);
        wp = an(N+2:end);
        phi = -wp/k/(1i);
end
% fprintf('z_c = %.4f\n', z(minpos));
%% Plot phi
figure('position',[0 0 1440 1280]);
subplot(1,2,1);
plot(real(phi),z,'-k','linewidth',1);
% hold on; plot(real(phi(minpos)),z(minpos),'ro'); hold off;
xlabel('real(\phi)','FontSize',30);
ylabel('z','FontSize',30);
ylim([-h(k) 0]);
title('real','FontSize',30);
grid on;
set(gca,'fontsize',24);
subplot(1,2,2);
plot(imag(phi),z,'-k','linewidth',1);
% hold on; plot(imag(phi(minpos)),z(minpos),'ro'); hold off;
xlabel('imag(\phi)','FontSize',30);
ylabel('z','FontSize',30);
ylim([-h(k) 0]);
title('imag','FontSize',30);
set(gca,'fontsize',24);
grid on;
sgtitle('\phi','FontSize',30);
%% Plot u
figure('position',[0 0 1440 1280]);
subplot(1,2,1);
plot(real(up),z,'.-k','linewidth',1);
% hold on; plot(real(up(minpos)),z(minpos),'ro'); hold off;
xlabel('real(u'')','FontSize',30);
ylabel('z','FontSize',30);
ylim([-h(k) 0]);
title('real','FontSize',30);
grid on;
set(gca,'fontsize',24);
subplot(1,2,2);
plot(imag(up),z,'.-k','linewidth',1);
% hold on; plot(imag(up(minpos)),z(minpos),'ro'); hold off;
xlabel('imag(u'')','FontSize',30);
ylabel('z','FontSize',30);
ylim([-h(k) 0]);
title('imag','FontSize',30);
set(gca,'fontsize',24);
grid on;
sgtitle('u''','FontSize',30);
%% Plot w
figure('position',[0 0 1440 1280]);
subplot(1,2,1);
plot(real(wp),z,'k','linewidth',1);
% hold on; plot(real(wp(minpos)),z(minpos),'ro'); hold off;
xlabel('real(w'')','FontSize',30);
ylabel('z','FontSize',30);
ylim([-h(k) 0]);
title('real','FontSize',30);
set(gca,'fontsize',24);
grid on;
subplot(1,2,2);
plot(imag(wp),z,'k','linewidth',1);
% hold on; plot(imag(wp(minpos)),z(minpos),'ro'); hold off;
xlabel('imag(w'')','FontSize',30);
ylabel('z','FontSize',30);
ylim([-h(k) 0]);
title('imag','FontSize',30);
set(gca,'fontsize',24);
grid on;
sgtitle('w''','FontSize',30);