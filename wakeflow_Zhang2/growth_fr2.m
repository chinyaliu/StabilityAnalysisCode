close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray", "D4"]; %, "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];%, "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
N = 400;
dk = 0.05;
k = dk:dk:4;
Re = inf;
Fr2 = (0.5:0.5:2.5).^2;
inflec_pt = -0.74708299;
cutz(1) = -inflec_pt;
h1 = @(k) 2*pi/k;
h2 = @(k) 6.5;
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
for j = 1:length(Fr2)
    for i = 1:length(k)
        if (k(i) < pi/3)
            h = h1;
        else
            h = h2;
        end
        [o(i), ~, ~, ~, ~, z(:,i), phi{i}, z_c(i), cutz(i+1)] = wZhang_solver2(N,k(i),h(k(i)),Re,Fr2(j),method,alg,do_balancing,cutz(i),'y');
        fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o(i)));
    end
    base(j,:) = {o, z, phi, z_c, cutz(2:end)};
end
% save('modeshape','phi','z','N','k','cutz','o','z_c');
%% Plot growth rate v.s. omega_r
o = cell2mat(base(:,1));
fig1 = figure('position',[50,0,1000,720]);
plot(real(o(1,:)),imag(o(1,:)),'linewidth',1);
hold on;
for i = 2:length(Fr2)
    plot(real(o(i,:)),imag(o(i,:)),'linewidth',1);
end
% ylim([0 0.04]);
% xlim([0 4]);
set(gca,'fontsize',20);
xlabel('$\tilde{\omega_r}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{\omega_i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
ax.YAxis.Exponent = -2;
grid on;