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
Ni = 600;
Re = inf;
Fr2 = 2.25;
h = @(k) 2*pi/k;
zL = 0.74708299;
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
% N = 20:10:600;
N = 20:20:1200;
% k = [0.01 0.2 0.4 0.6];
k = [0.87 1 2 3 4];
tic;
% case1 = wZhang_solver(N(1)/2,1,1,Re,Fr2,method);
case1 = wZhang_block(N,k(1),h(k(1)),Re,Fr2,method);
for j = 1:length(k)
    fprintf('k = %.2f\n',k(j));
%     case1.k = k(j); case1.h = h(k(j));
    case1.k = k(j); case1.h = h(k(j)); 
%     case1.N = Ni;
%     case1.solver(zL, 'y', alg, do_balancing);
%     case1.solver(zL, 'y', alg, do_balancing,0.05);
%     zLn = case1.zL;
    for i = 1:length(N)
        case1.N = N(i);
%         o = case1.solver(zL,'y',alg,do_balancing);
        o = case1.solver(zL, 'y', alg, do_balancing,0.05);
        oi(i) = imag(o(1));
        fprintf('N = %3d, growth rate = %.8f\n', N(i), oi(i));
    end
    oiall(j,:) = oi;
    doi{j} = abs(diff(oi));
end
toc;
% %% Plot figure
% fig1 = figure('position',[50,0,1000,720]);
% % doi = abs(diff(oi));
% nam = sprintf('k = %.2f',k(1));
% semilogy(N(1:end-1), doi{1}, '-o', 'linewidth', 1, 'Displayname', nam);
% hold on;
% for i = 2:length(k)
%     nam = sprintf('k = %.2f',k(i));
%     semilogy(N(1:end-1), doi{i}, '-o', 'linewidth', 1, 'Displayname', nam);
% end
% hold off;
% set(gca,'fontsize',20);
% xlabel('$N$', 'Interpreter', 'LaTeX','fontsize',30);
% ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{m+1})\ |$', 'Interpreter', 'LaTeX','fontsize',30);
% legend('location','northeast');
% grid on;
% % exportgraphics(fig1, 'fig_convergence\branch2.png');
%% Plot figure 2
dall = abs(oiall - oiall(:,end).*ones(size(oiall,1),size(oiall,2)));
fig1 = figure('position',[50,0,1000,720]);
nam = [sprintf('$k = %.2f, ',k(1)) '\' sprintf('omega = %0.9e$',oiall(1,end))];
semilogy(N, dall(1,:), '-o', 'linewidth', 1, 'Displayname', nam);
hold on;
for i = 2:length(k)
    nam = [sprintf('$k = %.2f, ',k(i)) '\' sprintf('omega = %0.9e$',oiall(i,end))];
    semilogy(N, dall(i,:), '-o', 'linewidth', 1, 'Displayname', nam);
end
hold off;
set(gca,'fontsize',20);
xlabel('$N$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{end})\ |$', 'Interpreter', 'LaTeX','fontsize',30);
legend('location','northeast','interpreter','latex');
grid on;
% exportgraphics(fig1, 'fig_convergence\branch2.png');