close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray", "D4"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];
solveGEPmethod = ["qr", "qz", "eig"];
%% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1; % solveGEPmethod
do_balancing = 'n';
N = 400;
Re = inf;
Fr2 = 2.25;
delt = 0.04;
h = @(k) 2*pi/k;
zL = 0.74708299;
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
m = linspace(0.5,10,100);
k = 0.3;
tic;
case1 = wZhang_solver(N,k,1,Re,Fr2,method);
% case1 = wZhang_complex(N,k,1,Re,Fr2,method,delt);
for j = 1:length(m)
    case1.h = m(j).*h(k);
    o = case1.solver(zL,'y',alg,do_balancing);
%     o = case1.solver(alg,do_balancing);
    oi(j,:) = imag(o);
    fprintf('m = %3d, growth rate = %.8f\n', m(j), oi(j));
end
doi = abs(diff(oi));
toc;
%% Plot figure
fig1 = figure('position',[50,0,1000,720]);
semilogy(m(1:end-1), doi(:,1), '-o', 'linewidth', 1);
set(gca,'fontsize',20);
xlabel('$n\lambda$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\ | \ \omega_i(h_m) - \omega_i(h_{m+1})\ |$', 'Interpreter', 'LaTeX','fontsize',30);
% legend('location','northeast');
grid on;
% exportgraphics(fig1, 'fig_convergence\branch2.png');
%% Plot figure
fig2 = figure('position',[50,0,1000,720]);
semilogy(m, abs(oi(:,1)-oi(end,1)), '-o', 'linewidth', 1);
set(gca,'fontsize',20);
xlabel('$n\lambda$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\ | \ \omega_i(h_m) - \omega_i(h_{end})\ |$', 'Interpreter', 'LaTeX','fontsize',30);
% legend('location','northeast');
grid on;