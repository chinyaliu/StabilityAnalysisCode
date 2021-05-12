close all; clear all;% clc
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver1 = [1,1,1]; % [order, diff_method, constructAB_method]
solver2 = [2,2,2];
algorithm = 1;
do_balancing = 'n';
N = 800;
k = linspace(0.01,4,400);
Re1 = 1000;
Re2 = inf;
Fr2 = 2.25;
h = @(k) 2*pi/k;
%% Run solver
method1 = [order(solver1(1)), diff_method(solver1(2)), constructAB_method(solver1(3))];
% method2 = [order(solver2(1)), diff_method(solver2(2)), constructAB_method(solver2(3))];
alg = solveGEPmethod(algorithm);
for i = 1:length(k)
    [o, ~, cA1(i), errGEP1(i)] = wZhang_solver(N,k(i),h(k(i)),Re2,Fr2,method1,alg);
    o1(i,:) = imag(o);
    o1r(i,:) = real(o);
%     fprintf('k = %.2f, growth rate = %.4f\n', k(i), o1(i));
    fprintf('k = %.2f, growth rate = %.4f\n', k(i), real(o)/k(i));
%     [o, ~, cA2(i), errGEP2(i), db2(i)] = wZhang_solver(N,k(i),h,Re2,Fr2,method2,alg,do_balancing);
%     o2(i,:) = imag(o);
end
o1(isnan(o1)) = nan;
o1r(imag(o1)<0) = nan;
%% Plot growth rate v.s. k
% fig1 = figure('position',[50,0,1000,720]);
fig1 = figure;
plot(k,o1,'-ro','linewidth',1,'markersize',3);
% plot(k,o1,'-k','linewidth',3,'DisplayName','original');
ylim([0 0.04]);
xlim([0 4]);
set(gca,'fontsize',24);
xlabel('$\tilde{k}$','fontsize',30);
ylabel('$\tilde{\omega_i}$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
ax.YAxis.Exponent = -2;
% legend('location','northeast');
grid on;
%% Plot omega_r v.s. k
fig2 = figure('position',[50,0,1000,720]);
plot(k,o1r,'linewidth',1,'DisplayName','Re = 1000');
% hold on;
% plot(k,o2,'linewidth',1,'DisplayName','Re = inf');
% hold off;
% ylim([0 0.04]);
% xlim([0 4]);
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{\omega_r}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% legend('location','northeast');
grid on;
%% Plot c_r v.s. k
fig3 = figure('position',[50,0,1000,720]);
plot(k,o1r./k.','linewidth',1,'DisplayName','Re = 1000');
% hold on;
% plot(k,o2,'linewidth',1,'DisplayName','Re = inf');
% hold off;
% ylim([0 0.04]);
% xlim([0 4]);
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{c_r}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% legend('location','northeast');
ylim([0 1]);
grid on;