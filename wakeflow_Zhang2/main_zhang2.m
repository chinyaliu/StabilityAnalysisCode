close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray", "D4", "Ray_match"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];
solveGEPmethod = ["qr", "qz", "eig"];
%% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1; % solveGEPmethod
do_balancing = 'n';
N = 400;
k = 3;
Re = inf;
Fr2 = 2.25;
if k < pi/3
    h = @(k) 2*pi/k;
else
    h = @(k) 6.5;
end
zL = 0.74708299;
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
t1 = tic;
[o, an, cA, errGEP, db, z, phi, zc, zLo] = wZhang_solver2(N,k,h(k),Re,Fr2,method,alg,do_balancing,zL,'y');
toc(t1);
%% Plot
for i = 1:3
figure('position',[0 0 1680 960]);
subplot(1,4,1);
plot(abs(phi(:,i)),z,'-k.','linewidth',1,'markersize',10);
hold on;
yline(-zLo, '-.r', 'linewidth', 1.5);
yline(-zL, '--b', 'linewidth', 1.5);
hold off;
set(gca,'fontsize',20);
if (h(k) > 6)
    blim = -6;
else
    blim = fix(-h(k));
end
xlabel('$abs$','FontSize',30, 'Interpreter', 'LaTeX');
ylabel('$z$','FontSize',30, 'Interpreter', 'LaTeX');
ylim([blim 0]);
% yax = sort([fix(-h(k)):0.5:0 -zL]); 
% yticks(yax);
% ax = gca;
% ax.YTickLabel{yax == -zL} = ['\color{red}' ax.YTickLabel{yax == -zL}];
grid on;
subplot(1,4,2);
ang_phi = unwrap(angle(phi(:,i)));
% rn = real(phi(:,i)) < 0;
% in = imag(phi(:,i)) < 0;
% ang_phi(rn & in) = ang_phi(rn & in) + 2*pi;
plot(ang_phi,z,'-k.','linewidth',1,'markersize',10);
hold on;
yline(-zLo, '-.r', 'linewidth', 1.5);
yline(-zL, '--b', 'linewidth', 1.5);
hold off;
set(gca,'fontsize',20);
ylim([blim 0]);
xlabel('$angle$','FontSize',30, 'Interpreter', 'LaTeX');
grid on;

subplot(1,4,3);
plot(real(phi(:,i)),z,'-k.','linewidth',1,'markersize',10);
hold on;
yline(-zLo, '-.r', 'linewidth', 1.5);
yline(-zL, '--b', 'linewidth', 1.5);
hold off;
set(gca,'fontsize',20);
ylim([blim 0]);
xlabel('$real$','FontSize',30, 'Interpreter', 'LaTeX');
% yticks(yax);
% ax = gca;
% ax.YTickLabel{yax == -zL} = ['\color{red}' ax.YTickLabel{yax == -zL}];
grid on;

subplot(1,4,4);
plot(imag(phi(:,i)),z,'-k.','linewidth',1,'markersize',10);
hold on;
yline(-zLo, '-.r', 'linewidth', 1.5);
yline(-zL, '--b', 'linewidth', 1.5);
hold off;
set(gca,'fontsize',20);
xlabel('$imag$','FontSize',30, 'Interpreter', 'LaTeX');
ylim([blim 0]);
grid on;
end