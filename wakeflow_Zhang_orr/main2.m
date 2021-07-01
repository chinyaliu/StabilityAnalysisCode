clear all;% clc
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
% Inputs
do_balancing = 'n';
eig_spectrum = 'all';
N = 600;
k = 0.01;
Re = inf;
Fr2 = 1.5^2;
h = 2*pi/real(k);
% h = 6;
c0 = 1./sqrt(k*Fr2);
% zL = wZhang_ddm.criticalH(c0); 
zL = 0.74708299;
% DDM numbers
numberofDDM = 4;
eps = 0.01;
f = wZhang_ddm.ddmtype(numberofDDM);

%% Run solver
t1 = tic;
case1 = wZhang_ddm(N,k,h,Re,Fr2);
case1.numMeth(method);
[o, an] = case1.solver(alg, do_balancing, eig_spectrum, f, struct('zL1',zL,'eps',eps));
c = o/k;
% Run second time with different N
case1.N = N+50;
o2 = case1.solver(alg, do_balancing, eig_spectrum, f, struct('zL1',zL,'eps',eps));
toc(t1);

%% Compare between eigenvalues of N and N2
O = repmat(o,1,length(o2));
O2 = repmat(o2.',length(o),1);
o_min = min(abs((O-O2)./O2),[],2);
o_chosen = o(o_min<1e-5);

%% Plot two eigenvalue spectrums
figure;
scatter(real(o),imag(o),30,'+k','DisplayName',sprintf('N = %d',N));
hold on;
scatter(real(o2),imag(o2),30,'xb','DisplayName',sprintf('N = %d',N+50));
scatter(real(o_chosen),imag(o_chosen),30,'or','DisplayName','Selected');
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off; box on;
xlabel('$\omega _r$');
ylabel('$\omega _i$','rotation',0, 'HorizontalAlignment','right');
legend('location','southeast');
title(sprintf('k = %.2f',k));

%% Plot oi_or
figure;
scatter(real(o),imag(o),'ok');
hold on;
yline(0,'linewidth',1.5,'color','#898989');
scatter(real(o_chosen),imag(o_chosen),'b','filled');
line([0,real(k)],[0,imag(k)],'color','r','linewidth',1);
hold off;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
titext = sprintf('$k=%.2f%+.2fi$',real(k),imag(k));
title(titext);

%% Plot ci_cr
c = o/k;
c2 = o2/k;
c_chosen = o_chosen/k;
figure;
g(1) = scatter(real(c),imag(c),100,'+k','DisplayName',sprintf('N = %d',N));
hold on;
g(2) = scatter(real(c2),imag(c2),100,'xb','DisplayName',sprintf('N = %d',N+50));
g(3) = scatter(real(c_chosen),imag(c_chosen),100,'or','DisplayName','Selected');
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off;box on;
xlabel('$c_r$');
ylabel('$c_i$','rotation',0, 'HorizontalAlignment','right');
legend('location','southwest');
titext = sprintf('$k=%.2f%+.2fi$',real(k),imag(k));
title(titext);