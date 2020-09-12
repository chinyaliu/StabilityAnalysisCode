close all; clear all;% clc
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["Origin", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver1 = [2,1,2]; % [order, diff_method, constructAB_method]
solver2 = [2,2,2];
algorithm = 1;
do_balancing = 'n';
N = 800;
dk = 0.01;
k = dk:dk:1.2;
Re1 = 1000;
Re2 = 100000;
%% Run solver
method1 = [order(solver1(1)), diff_method(solver1(2)), constructAB_method(solver1(3))];
method2 = [order(solver2(1)), diff_method(solver2(2)), constructAB_method(solver2(3))];
alg = solveGEPmethod(algorithm);
for i = 1:length(k)
    [o, ~, cA1(i), errGEP1(i), db1(i)] = poiseuille_solver(N,k(i),Re2,method1,alg,do_balancing);
    o1(i,:) = imag(o);
    fprintf('k = %.2f, growth rate = %.4f\n', k(i), o1(i));
    [o, ~, cA2(i), errGEP2(i), db2(i)] = poiseuille_solver(N,k(i),Re2,method2,alg,do_balancing);
    o2(i,:) = imag(o);
end
% oi(isnan(oi)) = 0;
% oinf(isnan(oinf)) = 0;
%% Plot growth rate v.s. k
% fig1 = figure('position',[50,0,1000,720]);
% plot(k,oi,'linewidth',1,'DisplayName','Re = 1000');
% hold on;
% plot(k,oinf,'linewidth',1,'DisplayName','Re = inf');
% hold off;
% ylim([0 0.04]);
% xlim([0 4]);
% xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',14);
% ylabel('$\tilde{\omega_i}$', 'Interpreter', 'LaTeX','fontsize',14);
% ttext = sprintf('$Fr^2 = %.2f$', Fr2);
% title(ttext,'fontsize',16,'interpreter','latex');
% legend('location','northeast');
% grid on;