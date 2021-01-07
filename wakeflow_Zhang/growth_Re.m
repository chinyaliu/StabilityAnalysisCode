close all; clear all;% clc
%% Solver & Algorithm list
wZhangMethod = ["DDM", "Complex"];
order = ["Ray","D4"];
diff_meth = ["Schimd", "Trefethen"];
makeAB_meth = ["D4", "Schimd"];
solveGEPmeth = ["qr", "qz", "eig"];
%% Set solver
meth = wZhangMethod(1);
method = [order(1), diff_meth(1), makeAB_meth(1)];
alg = solveGEPmeth(1);
%% Inputs
do_balancing = 'n';
N = 400;
k = linspace(0.01,4,400);
Re = inf;
Fr2 = 2.25;
h = 2*pi./k;
eps = 0.15;
c0 = 1./sqrt(k*Fr2);
zL = real(wZhang_ddm.g(c0)); % zL = 0.74708299;
addvar = struct('zL1',zL(1),'eps',eps);
switch lower(meth)
    case 'ddm' % Additional parameters for DDM
        wZhang = @wZhang_ddm;
        numberofDDM = 4;
        f = wZhang_ddm.ddmtype(numberofDDM);
        in_init = {N,k,h,Re,Fr2,method};
        in_solver = {alg, do_balancing, f, addvar};
    case 'complex' % Additional parameters for complex
        wZhang = @wZhang_complex;
        delt = 0.02;
        in_init = {N,k,h,Re,Fr2,method,delt};
        in_solver = {alg, do_balancing};
    otherwise
        error('Method not defined');
end
%% Run solver
tic;
p1 = wZhang(in_init{:});
o = NaN(1,length(k)); z_c = NaN(1,length(k));
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i); 
    if i ~= 1
        addvar.zL1 = -z_c(i-1);%zL(i);
    end
%     addvar.zL1 = zL(i);
    in_solver = {alg, do_balancing, f, addvar};
    oall = p1.solver(in_solver{:});
    o(i) = oall(1);
    z_c(i) = p1.zc;
    fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o(i)));
end
toc;
%% Plot growth rate v.s. k
fig1 = figure('position',[50,0,1000,720]);
plot(k,imag(o),'linewidth',2);
% ylim([0 0.04]);
xlim([0 4]);
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{\omega_i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
ax.YAxis.Exponent = -2;
grid on;
%% Plot z_c v.s. k
z_c(isnan(o)) = NaN;
fig2 = figure('position',[50,0,1000,720]);
plot(k,z_c,'linewidth',2);
inflec_pt = -0.74708299;
hold on; 
yline(inflec_pt, '-.r', 'linewidth', 2); 
plot(k,-real(zL),'--k','linewidth',2);
hold off;
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{z_c}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0);
grid on;
yt = sort([-3.5:0.5:0 inflec_pt]);
yticks(yt);
ind = find(yt==inflec_pt);
ax = gca;
ax.YTickLabel{ind} = ['\color{red}' ax.YTickLabel{ind}];
xlim([0 4]);