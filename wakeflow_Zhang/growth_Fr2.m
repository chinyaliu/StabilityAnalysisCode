clear all;% clc
%% Solver & Algorithm list
order = ["Ray","D4"];
diff_meth = ["Schimd", "Trefethen"];
makeAB_meth = ["D4", "Schimd"];
solveGEPmeth = ["qr", "qz", "eig"];
% Set solver
method = [order(1), diff_meth(1), makeAB_meth(1)];
alg = solveGEPmeth(1);
do_balancing = 'n';
N = 400;
k = linspace(0.01,4,100);
Re = inf;
Fr = [0.5 1.5 2.5 5.5 inf];
Fr2 = Fr.^2;
inflec_pt = -0.74708299;
cutz = NaN(1,length(k)+1);
zL = 0.74708299*ones(length(k),1);
cutz(1) = -inflec_pt;
h = 2*pi./k;
eps = 0.15;
inflec_pt = -0.74708299;
addvar = struct('zL1',zL(1),'eps',eps);
numberofDDM = 4;
f = wZhang_ddm.ddmtype(numberofDDM);
in_init = {N,k,h,Re,Fr2,method};
%% Run solver
tic;
p1 = wZhang_ddm(in_init{:});
o = NaN(1,length(k));
for j = 1:length(Fr2)
    fprintf('Fr = %.1f\n', Fr(j));
    p1.Fr2 = Fr2(j);
    for i = 1:length(k)
        p1.k = k(i); p1.h = h(i);
        addvar.zL1 = zL(i);
        o(i) = p1.solver(alg, do_balancing, f, addvar);
        fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o(i)));
    end
    o(isnan(imag(o))) = 0;
    frcase(j) = struct('o',o);
end
toc;
%% Plot growth rate v.s. k
% cl = {'#0E054D','#0E0D63','#182279','#273E8F','#395EA6','#4E80BC','#65A5D2','#81CAE9'};
cl = {'#5ed5f2','#27a0ef','#072ec5','#020294','k'};
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
% save('growth_fr','N','k','Fr','frcase');