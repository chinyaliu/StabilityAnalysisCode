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
Ni = 600;
% k = 4;
Re = inf;
Fr2 = 2.25;
h = @(k) 2*pi/k;
zL = 0.74708299;
%% Find inflection point
% c1 = 0.9988; c2 = 0.8814;
% syms z
% uzz(z) = diff(diff(1-c1*cosh(c2*z).^(-2)));
% inflec_pt = double(solve(uzz,'MaxDegree',3));
% zL = -inflec_pt(imag(inflec_pt)==0 & real(inflec_pt) < 0);
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
N = 20:10:600;
% k = [0.01 0.2 0.4 0.6];
k = [0.87 1 2 3 4];
for j = 1:length(k)
    fprintf('k = %.2f\n',k(j));
    [~, ~, ~, ~, ~, ~, ~, ~, zL] = wZhang_solver2(Ni,k(j),h(k(j)),Re,Fr2,method,alg,do_balancing,zL,'y');
    for i = 1:length(N)
        [o, an, cA, errGEP, db, z, ~, ~, ~] = wZhang_solver2(N(i),k(j),h(k(j)),Re,Fr2,method,alg,do_balancing,zL,'n');
        oi(i,:) = imag(o);
        or(i,:) = real(o);
        cr = real(o)/k(j);
        fprintf('N = %3d, growth rate = %.8f\n', N(i), oi(i));
    end
    doi{j} = abs(diff(oi));
end
%% Plot figures
fig1 = figure('position',[50,0,1000,720]);
% doi = abs(diff(oi));
nam = sprintf('k = %.2f',k(1));
semilogy(N(1:end-1), doi{1}, '-o', 'linewidth', 1, 'Displayname', nam);
hold on;
for i = 2:length(k)
    nam = sprintf('k = %.2f',k(i));
    semilogy(N(1:end-1), doi{i}, '-o', 'linewidth', 1, 'Displayname', nam);
end
hold off;
set(gca,'fontsize',20);
xlabel('$N$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{m+1})\ |$', 'Interpreter', 'LaTeX','fontsize',30);
% titext = sprintf('z_L = %.8f', zL);
% text(350,1e-4,titext,'fontsize',30);
legend('location','northeast');
grid on;
% %% Run solver
% zL = 0.74708299;
% N = 20:10:600;
% for i = 1:length(N)
%     [o, an, cA, errGEP, db, z, ~, ~, z_c(i)] = wZhang_solver2(N(i),k,h(k),Re,Fr2,method,alg,do_balancing,zL,'y');
%     oi(i,:) = imag(o);
%     or(i,:) = real(o);
%     cr = real(o)/k;
%     fprintf('N = %3d, growth rate = %.8f\n', N(i), oi(i));
% end
% %% Plot figures
% fig1 = figure('position',[50,0,1000,720]);
% doi = abs(diff(z_c));
% semilogy(N(1:end-1), doi, '-ko', 'linewidth', 1);
% set(gca,'fontsize',20);
% xlabel('$N$', 'Interpreter', 'LaTeX','fontsize',30);
% ylabel('$\ | \ z_L(N_m) - z_L(N_{m+1})\ |$', 'Interpreter', 'LaTeX','fontsize',30);
% % titext = sprintf('z_L = %.8f', zL);
% % text(350,1e-4,titext,'fontsize',30);
% grid on;