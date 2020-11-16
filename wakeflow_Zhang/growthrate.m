close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray", "D4"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];%, "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
N = 400;
k = linspace(0.01,4,400);
Re = inf;
Fr2 = 2.25;
inflec_pt = -0.74708299;
cutz = NaN(1,length(k)+1);
cutz(1) = -inflec_pt;
h = 6.5*ones(1,length(k));
h(k<pi/3) = 2*pi./k(k<pi/3);
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
tic;
p1 = wZhang_solver(N,1,1,Re,Fr2,method);
o = NaN(1,length(k)); z_c = NaN(1,length(k));
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i);
%     [o(i), an] = p1.solver(cutz(i), 'y', alg, do_balancing);
%     z(:,i) = p1.z; phi{i} = p1.phi; 
    o(i) = p1.solver(cutz(i), 'y', alg, do_balancing);
    z_c(i) = p1.zc; cutz(i+1)=p1.zL;
    fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o(i)));
end
toc;
cutz = cutz(2:end);
%% Plot growth rate v.s. k
fig1 = figure('position',[50,0,1000,720]);
plot(k,imag(o),'linewidth',2);
ylim([0 0.04]);
% xlim([0 4]);
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{\omega_i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
ax.YAxis.Exponent = -2;
grid on;
%% Plot z_c v.s. k
fig2 = figure('position',[50,0,1000,720]);
plot(k,z_c,'linewidth',2);
hold on; yline(inflec_pt, '-.r', 'linewidth', 2); hold off;
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{z_c}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0);
grid on;
yt = sort([-3.5:0.5:0 inflec_pt]);
yticks(yt);
ind = find(yt==inflec_pt);
ax = gca;
ax.YTickLabel{ind} = ['\color{red}' ax.YTickLabel{ind}];
xlim([0 4.5]);
ylim([-3.6 0]);
%% Save data & figures
% save('diffk','phi','z','N','k','cutz','o','z_c');
exportgraphics(fig1, 'fig_growthrate\omega_i.png');
exportgraphics(fig2, 'fig_growthrate\criticalheight.png');