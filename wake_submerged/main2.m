clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,~,N,H,k,Fr2,Re,eps,c0,h,f] = pars_wake;
eig_spectrum = 'all';
Nad = 200;

%% Run solver
t1 = tic;
case1 = wSubmerged(N,H,k,h,Re,Fr2,bflow);
case1.numMeth(method);
zL = 0.74708299+H;
addvar = struct('zL1',case1.invbf(c0),'eps',eps);
o = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
% Run second time with different N
case1.setprop('N',N+Nad);
o2 = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
o = o(real(o)>-50); o2 = o2(real(o2)>-50);
toc(t1);

%% Compare between eigenvalues of N and N2
O = repmat(o,1,length(o2));
O2 = repmat(o2.',length(o),1);
o_min = min(abs((O-O2)./O2),[],2);
o_chosen = o(o_min<1e-5);

%% Plot two eigenvalue spectrums
% figure;
fig2 = figure('position',[50,0,810,540]);
scatter(real(o),imag(o),100,'+k','DisplayName',sprintf('N = %d',N));
hold on;
scatter(real(o2),imag(o2),100,'xb','DisplayName',sprintf('N = %d',N+Nad));
scatter(real(o_chosen),imag(o_chosen),100,'or','DisplayName','Selected');
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off; box on;
xlabel('$\omega _r$');
ylabel('$\omega _i$','rotation',0, 'HorizontalAlignment','right');
legend('location','southeast');
title(sprintf('$k$ = %.2f',k));

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
g(1) = scatter(real(c),imag(c),150,'+k','DisplayName',sprintf('N = %d',N));
hold on;
g(2) = scatter(real(c2),imag(c2),150,'xb','DisplayName',sprintf('N = %d',N+Nad));
g(3) = scatter(real(c_chosen),imag(c_chosen),150,'or','DisplayName','Selected');
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off;box on;
xlabel('$c_r$');
ylabel('$c_i$','rotation',0, 'HorizontalAlignment','right');
legend('location','southwest');
titext = sprintf('$k=%.2f$',real(k));
% titext = sprintf('$k=%.2f%+.2fi$',real(k),imag(k));
title(titext);