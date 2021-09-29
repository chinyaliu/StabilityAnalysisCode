clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,delta_nd,lambda_nd,c0,h,f] = pars_Morland(1);
fprintf('u_d = %1.2f, delta = %1.3f, lambda = %1.3f\n',ud_nd,delta_nd,lambda_nd);

%% Run solver
t1 = tic;
case1 = wMorland(N,h,ud_nd,delta_nd,lambda_nd,method,bflow);
addvar = struct('zL1',-case1.criticalH(c0),'eps',0.2);
c = case1.solvers(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
case1.N = N+20;
c2 = case1.solvers(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
toc(t1);
if strcmpi(de_singularize,'y')
    c = c(real(c)>-100);
    c2 = c2(real(c2)>-100);
end

%% Compare between eigenvalues of N and N2
C = repmat(c,1,length(c2));
C2 = repmat(c2.',length(c),1);
c_min = min(abs((C-C2)./C2),[],2);
c_chosen = c(c_min<1e-6);

%% Plot eigenvalue spectrum ci_cr
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
titext = sprintf('$\\hat{u_d}=%1.2f,\\ \\hat{\\delta}=%1.3f, \\ \\hat{\\lambda}=%1.3f$',ud_nd,delta_nd,lambda_nd);
title(titext);