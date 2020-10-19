close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray", "D4"]; %, "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];%, "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1; % solveGEPmethod
do_balancing = 'n';
N = 600;
Re = inf;
Fr2 = 2.25;
h = @(k) 2*pi/k;
% h = @(k) 6;
zL = 0.74708299;
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Base flow
c1 = 0.9988; c2 = 0.8814;
U = @(z) (1-c1*cosh(c2*z).^(-2));
Uzz = @(z) 2*c1*c2^2*(sech(c2*z)).^2.*((sech(c2*z)).^2-2*(tanh(c2*z)).^2);
%% Branch 1
k = 0.01:0.01:0.6;
for i = 1:length(k)
    [o, an, ~, ~, ~, z, phi, zc, ~] = wZhang_solver2(N,k(i),h(k(i)),Re,Fr2,method,alg,do_balancing,zL,'y');
    [~,pos] = min(abs(z+zL));
    zin(i) = z(pos);
    p2(i) = phi(pos,3);
end
fig1 = figure;
plot(k, real(p2), '-o');
fig2 = figure;
plot(k, imag(p2), '-o');
%% Plot
% c = o/k;
% [~,pos] = min(abs(z+zL));
% zin = z(pos);
% Uin = U(zin);
% Uzzin = Uzz(zin);
% pk2 = k*k*phi(pos,1);
% p2 = phi(pos,3);
% pdif = p2-pk2;
% infA = (Uin-c)*pdif;
% equ1 = (U(z)-c).*(phi(:,3)-k*k*phi(:,1));
% figure('position',[0 0 1920 1280]);
% subplot(1,4,1);
% plot(abs(phi(:,1)),z,'-k.','linewidth',1,'markersize',10);
% hold on;
% yline(zc, '-.r', 'linewidth', 1.5);
% yline(-zL, '--b', 'linewidth', 1.5);
% hold off;
% set(gca,'fontsize',20);
% if (h(k) > 6)
%     blim = -6;
% else
%     blim = fix(-h(k));
% end
% xlabel('$abs(\phi)$','FontSize',30, 'Interpreter', 'LaTeX');
% ylim([blim 0]);
% grid on;
% subplot(1,4,2);
% plot(abs(phi(:,3)),z,'-k.','linewidth',1,'markersize',10);
% hold on;
% yline(zc, '-.r', 'linewidth', 1.5);
% yline(-zL, '--b', 'linewidth', 1.5);
% hold off;
% xlabel('$\phi_{zz}$','FontSize',30, 'Interpreter', 'LaTeX');
% set(gca,'fontsize',20);
% grid on;
% ylim([blim 0]);
% subplot(1,4,3);
% plot(abs(U(z)-c),z,'-k.','linewidth',1,'markersize',10);
% hold on;
% yline(zc, '-.r', 'linewidth', 1.5);
% yline(-zL, '--b', 'linewidth', 1.5);
% hold off;
% xlabel('$abs(U-c)$','FontSize',30, 'Interpreter', 'LaTeX');
% set(gca,'fontsize',20);
% ylim([blim 0]);
% grid on;
% subplot(1,4,4);
% plot(abs(equ1),z,'-k.','linewidth',1,'markersize',10);
% hold on;
% yline(zc, '-.r', 'linewidth', 1.5);
% yline(-zL, '--b', 'linewidth', 1.5);
% hold off;
% xlabel('$abs(equ1)$','FontSize',30, 'Interpreter', 'LaTeX');
% set(gca,'fontsize',20);
% ylim([blim 0]);
% grid on;
% %% Branch 2
% k = 1.2;
% [o, an, cA, errGEP, db, z, phi, zc, ~] = wZhang_solver2(N,k,h(k),Re,Fr2,method,alg,do_balancing,zL,'y');
% c = o/k;
% [~,pos] = min(abs(z+zL));
% zin = z(pos);
% Uin = U(zin);
% Uzzin = Uzz(zin);
% pk2 = k*k*phi(pos,1);
% p2 = phi(pos,3);
% pdif = p2-pk2;
% infA = (Uin-c)*pdif;
% equ1 = (U(z)-c).*(phi(:,3)-k*k*phi(:,1));
% %% Plot
% figure('position',[0 0 1920 1280]);
% subplot(1,4,1);
% plot(real(phi(:,3)-k*k*phi(:,1)),z,'-k.','linewidth',1,'markersize',10);
% hold on;
% yline(zc, '-.r', 'linewidth', 1.5);
% yline(-zL, '--b', 'linewidth', 1.5);
% hold off;
% set(gca,'fontsize',20);
% if (h(k) > 6)
%     blim = -6;
% else
%     blim = fix(-h(k));
% end
% xlabel('$real(\phi_{zz}-k^2\phi)$','FontSize',30, 'Interpreter', 'LaTeX');
% ylim([blim 0]);
% grid on;
% subplot(1,4,2);
% plot(imag(phi(:,3)-k*k*phi(:,1)),z,'-k.','linewidth',1,'markersize',10);
% hold on;
% yline(zc, '-.r', 'linewidth', 1.5);
% yline(-zL, '--b', 'linewidth', 1.5);
% hold off;
% xlabel('$imag(\phi_{zz}-k^2\phi)$','FontSize',30, 'Interpreter', 'LaTeX');
% set(gca,'fontsize',20);
% grid on;
% ylim([blim 0]);
% subplot(1,4,3);
% plot(real(phi(:,3)),z,'-k.','linewidth',1,'markersize',10);
% hold on;
% yline(zc, '-.r', 'linewidth', 1.5);
% yline(-zL, '--b', 'linewidth', 1.5);
% hold off;
% xlabel('$\phi_{zzr}$','FontSize',30, 'Interpreter', 'LaTeX');
% set(gca,'fontsize',20);
% ylim([blim 0]);
% grid on;
% subplot(1,4,4);
% plot(imag(phi(:,3)),z,'-k.','linewidth',1,'markersize',10);
% hold on;
% yline(zc, '-.r', 'linewidth', 1.5);
% yline(-zL, '--b', 'linewidth', 1.5);
% hold off;
% xlabel('$\phi_{zzi}$','FontSize',30, 'Interpreter', 'LaTeX');
% set(gca,'fontsize',20);
% ylim([blim 0]);
% grid on;
% figure('position',[0 0 960 1280]);
% subplot(1,2,1);
% plot(real(equ1),z,'-k.','linewidth',1,'markersize',10);
% hold on;
% yline(zc, '-.r', 'linewidth', 1.5);
% yline(-zL, '--b', 'linewidth', 1.5);
% hold off;
% xlabel('$real[(U-c)*(\phi_{zz}-k^2\phi)]$','FontSize',30, 'Interpreter', 'LaTeX');
% set(gca,'fontsize',20);
% ylim([blim 0]);
% grid on;
% subplot(1,2,2);
% plot(imag(equ1),z,'-k.','linewidth',1,'markersize',10);
% hold on;
% yline(zc, '-.r', 'linewidth', 1.5);
% yline(-zL, '--b', 'linewidth', 1.5);
% hold off;
% xlabel('$imag[(U-c)*(\phi_{zz}-k^2\phi)]$','FontSize',30, 'Interpreter', 'LaTeX');
% set(gca,'fontsize',20);
% ylim([blim 0]);
% grid on;