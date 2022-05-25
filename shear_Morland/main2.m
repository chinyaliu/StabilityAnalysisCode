clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,~,N,ud_nd,delta_nd,lambda_nd,c0,h,f,epss,Re] = pars_Morland;
eig_spectrum = 'all';
N2 = N+100;

%% Run solver
t1 = tic;
case1 = wMorland(N,h,ud_nd,delta_nd,lambda_nd,method,bflow,Re);
case1.setprop('k',case1.k-0.1i);
addvar = struct('zL1',-0.26067,'eps',epss);
% addvar = struct('zL1',case1.invbf(c0),'eps',eps);
o = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
% Run second time with different N
case1.setprop('N',N2);
o2 = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
o = o(real(o)>-50); o2 = o2(real(o2)>-50);
k = case1.k;
c = oi./k;
c2 = o2./k;
toc(t1);

%% Compare between eigenvalues of N and N2
C = repmat(c,1,length(c2));
C2 = repmat(c2.',length(c),1);
c_min = min(abs((C-C2)./C2),[],2);
c_chosen = c(c_min<1e-5);
o_chosen = c_chosen*k;

%% Plot two eigenvalue spectrums
figure;
scatter(real(o),imag(o),30,'+k','DisplayName',sprintf('N = %d',N));
hold on;
scatter(real(o2),imag(o2),30,'xb','DisplayName',sprintf('N = %d',N2));
scatter(real(o_chosen),imag(o_chosen),30,'or','DisplayName','Selected');
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off; box on;
xlabel('$\omega _r$');
ylabel('$\omega _i$','rotation',0, 'HorizontalAlignment','right');
legend('location','southeast');
title(sprintf('$k=%.2f%+.2fi$',real(case1.k),imag(case1.k)));

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
titext = sprintf('$k=%.2f%+.2fi$',real(case1.k),imag(case1.k));
title(titext);

%% Plot ci_cr
figure;
g(1) = scatter(real(c),imag(c),100,'+k','DisplayName',sprintf('N = %d',N));
hold on;
g(2) = scatter(real(c2),imag(c2),100,'xb','DisplayName',sprintf('N = %d',N2));
g(3) = scatter(real(c_chosen),imag(c_chosen),100,'or','DisplayName','Selected');
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off;box on;
xlabel('$c_r$');
ylabel('$c_i$','rotation',0, 'HorizontalAlignment','right');
legend('location','southwest');
titext = sprintf('$k=%.2f%+.2fi$',real(case1.k),imag(case1.k));
title(titext);