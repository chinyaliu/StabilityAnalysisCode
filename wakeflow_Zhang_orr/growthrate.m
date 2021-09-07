% close all; clear all;% clc
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
do_balancing = 'y';
eig_spectrum = 'max';
Re = 1e7;
Fr2 = 2.25;
N = 600;
k = linspace(0.01,4,100);
h = 2*2*pi./k;
% h = 6*ones(1,length(k));
% h(h<10) = 10;
inflec_pt = -0.74708299;
% c0 = 1./sqrt(k*Fr2);
% zL = real(wZhang_ddm.g(c0)); 
zL = 0.74708299*ones(length(k),1);
cutz = NaN(1,length(k));
cutz(1) = -inflec_pt;
numberofDDM = 4;
eps = 0.01;
f = wZhang_ddm.ddmtype(numberofDDM);
in_init = {N,k(1),h(1),Re,Fr2};
addvar = struct('zL1',zL(1),'eps',eps);
%% Run solver
tic;
p1 = wZhang_ddm(in_init{:});
p1.numMeth(method);
o = NaN(1,length(k));
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i);
%     addvar.zL1 = zL(i);
    o(i) = p1.solver(alg, do_balancing, eig_spectrum, f, addvar);
    if isnan(p1.zc)
        cutz(i)=cutz(1);
    else
        cutz(i)=-p1.zc;
    end
    addvar.zL1 = cutz(i);
    fprintf('k = %.2f\n', k(i));
end
toc;

%% Plot oi vs or
fig1 = figure('position',[50,0,1000,720]);
plot(real(o),imag(o),'k.','Markersize',6);
hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
xlabel('$\tilde{\omega _r}$','fontsize',30);
ylabel('$\tilde{\omega _i}$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot oi vs k
fig2 = figure('position',[50,0,1000,720]);
plot(k,imag(o),'k.','Markersize',6);
hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
xlim([0 4]);
xlabel('$\tilde{k}$','fontsize',30);
ylabel('$\tilde{\omega _i}$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot critical height vs k
fig3 = figure('position',[50,0,1000,720]);
plot(k,-cutz,'k.','Markersize',6);
hold on; 
yline(inflec_pt,'linewidth',1.5,'color','b');
xlim([0 4]);
xlabel('$\tilde{k}$','fontsize',30);
ylabel('$\tilde{z_c}$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
%% Plot ci vs cr
c = o./k;
fig1 = figure('position',[50,0,1000,720]);
plot(real(c),imag(c),'k.','Markersize',6);
hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
xlabel('$\tilde{c _r}$','fontsize',30);
ylabel('$\tilde{c _i}$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');