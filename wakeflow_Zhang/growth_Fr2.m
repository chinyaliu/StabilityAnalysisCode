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
k = linspace(0.01,4,100);
Re = inf;
Fr = [0 0.5:2:6.5 inf];
Fr2 = Fr.^2;
inflec_pt = -0.74708299;
cutz = NaN(1,length(k)+1);
cutz(1) = -inflec_pt;
h = 6*ones(1,length(k));
h(k<pi/3) = 2*pi./k(k<pi/3);
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
tic;
p1 = wZhang_solver(N,1,1,Re,Fr2,method);
o = NaN(1,length(k)); z_c = NaN(1,length(k));
for j = 1:length(Fr2)
    fprintf('Fr = %.1f\n', Fr(j));
    p1.Fr2 = Fr2(j);
    for i = 1:length(k)
        p1.k = k(i); p1.h = h(i);
        o(i) = p1.solver(cutz(i), 'y', alg, do_balancing);
        z_c(i) = p1.zc; cutz(i+1)=p1.zL;
        fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o(i)));
    end
    o(isnan(imag(o))) = 0;
    frcase(j) = struct('cutz',cutz(2:end),'o',o,'z_c',z_c);
end
toc;
%% Plot growth rate v.s. k
% cl = {'#0E054D','#0E0D63','#182279','#273E8F','#395EA6','#4E80BC','#65A5D2','#81CAE9'};
cl = {'#5ed5f2','#27a0ef','#126ef1','#072ec5','#020294','k'};
ls = {'-','-.','-',':','-.','--'};
fig1 = figure('position',[50,0,1000,720]);
plot(k,imag(frcase(1).o),'linewidth',2,'Displayname',num2str(Fr(1)),'color',cl{1},'linestyle',ls{1});
hold on;
for i = 2:length(Fr)
    p = plot(k,imag(frcase(i).o),'linewidth',2,'Displayname',num2str(Fr(i)),'color',cl{i},'linestyle',ls{i});
end
hold off;
xlim([0 4]);
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{\omega_i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
ax.YAxis.Exponent = -2;
grid on;
leg = legend('location','northeast');
title(leg,'Fr');
%% Save data
save('growth_fr','N','k','Fr','frcase');