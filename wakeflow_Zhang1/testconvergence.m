close all; clear all;% clc
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1; % solveGEPmethod
do_balancing = 'n';
Re = inf;
Fr2 = 2.25;
h = @(k) 2*pi/k;
%% Run solver
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
N = 40:20:1200;
% k = [0.01 0.2 0.4 0.6];
k = [0.87 1 2 3];
for j = 1:length(k)
    fprintf('k = %.2f\n',k(j));
    for i = 1:length(N)
        [o, ~, ~, ~] = wZhang_solver(N(i),k(j),h(k(j)),Re,Fr2,method,alg);
        oi(i,:) = imag(o);
        or(i,:) = real(o);
        cr = real(o)/k(j);
        fprintf('N = %3d, growth rate = %.8f\n', N(i), oi(i));
    end
    doi{j} = abs((oi-oi(end))./oi(end));
end
%% Plot figures
fig1 = figure('position',[50,0,1000,720]);
% doi = abs(diff(oi));
nam = sprintf('k = %.2f',k(1));
% semilogy(N(1:end-1), doi{1}, '-o', 'linewidth', 1, 'Displayname', nam);
semilogy(N, doi{1}, '-o', 'linewidth', 1, 'Displayname', nam);
hold on;
for i = 2:length(k)
    nam = sprintf('k = %.2f',k(i));
%     semilogy(N(1:end-1), doi{i}, '-o', 'linewidth', 1, 'Displayname', nam);
    semilogy(N, doi{i}, '-o', 'linewidth', 1, 'Displayname', nam);
end
hold off;
set(gca,'fontsize',20);
xlabel('$N$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{end})\ |/\omega_i(N_{end})$','fontsize',30);
% titext = sprintf('z_L = %.8f', zL);
% text(350,1e-4,titext,'fontsize',30);
legend('location','northeast');
grid on;