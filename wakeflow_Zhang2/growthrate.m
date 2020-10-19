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
dk = 0.01;
k = dk:dk:4.5;
Re = inf;
Fr2 = 2.25;
inflec_pt = -0.74708299;
cutz(1) = -inflec_pt;
h1 = @(k) 2*pi/k;
h2 = @(k) 6.5;
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
for i = 1:length(k)
    if (k(i) < pi/3)
        h = h1;
    else
        h = h2;
    end
    [o(i), ~, ~, ~, ~, z(:,i), phi{i}, z_c(i), cutz(i+1)] = wZhang_solver2(N,k(i),h(k(i)),Re,Fr2,method,alg,do_balancing,cutz(i),'y');
    fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o(i)));
end
cutz = cutz(2:end);
% save('modeshape','phi','z','N','k','cutz','o','z_c');
%% Plot growth rate v.s. k
o1 = imag(o);
o1r = real(o);
o1(isnan(o1)) = nan;
o1r(imag(o1)<0) = nan;
fig1 = figure('position',[50,0,1000,720]);
plot(k,o1,'linewidth',1);
ylim([0 0.04]);
% xlim([0 4]);
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{\omega_i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
ax.YAxis.Exponent = -2;
grid on;
% %% Plot omega_r v.s. k
% fig2 = figure('position',[50,0,1000,720]);
% plot(k,o1r,'linewidth',1);
% set(gca,'fontsize',20);
% xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
% ylabel('$\tilde{\omega_r}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% grid on;
% %% Plot c_r v.s. k
% fig3 = figure('position',[50,0,1000,720]);
% plot(k,o1r./k,'linewidth',1);
% set(gca,'fontsize',20);
% xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
% ylabel('$\tilde{c_r}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% ylim([0 1]);
% grid on;
%% Plot z_c v.s. k
fig4 = figure('position',[50,0,1000,720]);
plot(k,z_c,'linewidth',1);
hold on; yline(inflec_pt, '-.r', 'linewidth', 1.5); hold off;
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{z_c}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
grid on;
yt = sort([-3.5:0.5:0 inflec_pt]);
yticks(yt);
ind = find(yt==inflec_pt);
ax = gca;
ax.YTickLabel{ind} = ['\color{red}' ax.YTickLabel{ind}];
xlim([0 4.5]);
ylim([-3.6 0]);
% %% Plot domain height v.s. k
% fig5 = figure('position',[50,0,1000,720]);
% plot(k,-cutz,'linewidth',1);
% set(gca,'fontsize',20);
% xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
% ylabel('$(\frac{d^2 \phi }{dz^2})_{max}$', 'Interpreter', 'LaTeX','fontsize',30);
% grid on;