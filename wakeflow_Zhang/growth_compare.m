close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray","D4"];
diff_meth = ["Schimd", "Trefethen"];
makeAB_meth = ["D4", "Schimd"];
solveGEPmeth = ["qr", "qz", "eig"];
% Set solver
method = [order(1), diff_meth(1), makeAB_meth(1)];
alg = solveGEPmeth(1);
% Inputs
do_balancing = 'n';
Re = inf;
Fr2 = 2.25;
N = 600;
k = linspace(0.01,4,400);
h = 2*pi./real(k);
eps = 0.15;
inflec_pt = -0.74708299;
c0 = 1./sqrt(k*Fr2);
% zL = real(wZhang_ddm.g(c0)); 
zL = 0.74708299*ones(length(k),1);
cutz = NaN(1,length(k)+1);
cutz(1) = -inflec_pt;
in_init = {N,k,h,Re,Fr2,method};
%% DDM 1
addvar = struct('zL1',zL(1));
numberofDDM = 1;
f = wZhang_ddm.ddmtype(numberofDDM);
tic;
p1 = wZhang_ddm(in_init{:});
o1 = NaN(1,length(k));
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i);
    oall{i} = p1.solver(alg, do_balancing, f, addvar);
    o1(i) = oall{i}(1);
    fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o1(i)));
end
toc;
%% DDM 4
addvar = struct('zL1',zL(1),'eps',eps);
numberofDDM = 2;
f = wZhang_ddm.ddmtype(numberofDDM);
tic;
p1.N = 400;
o4 = NaN(1,length(k)); z_c = NaN(1,length(k));
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i);
    addvar.zL1 = cutz(i);
    oall{i} = p1.solver(alg, do_balancing, f, addvar);
    o4(i) = oall{i}(1);
    %
    z_c(i) = p1.zc; 
    if isnan(z_c(i))
        cutz(i+1)=cutz(1);
    else
        cutz(i+1)=-p1.zc;
    end
    fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o4(i)));
end
toc;
%% Plot growth rate v.s. k
fig1 = figure('position',[50,0,1000,720]);
hold on;
plot(k,imag(o1),'-ro','linewidth',2,'markersize',2,'DisplayName','Original');
plot(k,imag(o4),'-bx','linewidth',2,'markersize',2,'DisplayName','DDM 4');
hold off;
ylim([0 0.04]);
xlim([0 4]);
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{\omega_i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
ax.YAxis.Exponent = -2;
grid on; box on;
legend('location','northeast');