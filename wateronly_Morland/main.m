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
% case1.k = case1.k-0.1i;
% [c,an,cA,errz,dob] = case1.solver(alg, do_balancing, eig_spectrum, f, addvar);
addvar = struct('zL1',-case1.criticalH(c0),'eps',0.1);
[c1,an,cA,errz,dob] = case1.solvers(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
toc(t1);
% Discard eigenvalues set by de-singularizing process
if strcmpi(de_singularize,'y')
    c = c1(real(c1)>-100);
else
    c = c1;
end

%% Choose eigenvalue
a = 1:length(c); 
crange = (real(c)>0) & (real(c)<ud_nd);
aa = a(crange);
abch = isoutlier(imag(c(aa)),'movmedian',5);
aa = [a(real(c)<0) aa(abch) a(real(c)>ud_nd)];
% aa = aa(abch);
c_chosen = c(aa);
an_c = an(:,aa);
o = 2*pi*c/lambda_nd;
o_chosen = 2*pi*c_chosen/lambda_nd;
fprintf('c_r = %1.3f, growth rate = %1.3f\n',real(c_chosen),imag(o_chosen));

%% Plot eigenvalue spectrum ci_cr
figure;
plot(real(c),imag(c),'ok');
hold on; 
yline(0,'linewidth',1.5,'color','#898989'); 
scatter(real(c_chosen),imag(c_chosen),'b','filled');
hold off;
xlabel('$c_r$','fontsize',30);
ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ymax = max(abs(imag(c)));
ylim([-1.5*ymax 1.5*ymax]);
titext = sprintf('$\\hat{u_d}=%1.2f,\\ \\hat{\\delta}=%1.3f, \\ \\hat{\\lambda}=%1.3f$',ud_nd,delta_nd,lambda_nd);
title(titext);

%% Plot eigenvalue spectrum oi_or
figure;
plot(real(o),imag(o),'ok');
hold on;       
yline(0,'linewidth',1.5,'color','#898989'); 
scatter(real(o_chosen),imag(o_chosen),'b','filled');
hold off;
xlabel('$\omega_r$','fontsize',30);
ylabel('$\omega_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ymax = max(abs(imag(o)));
ylim([-1.5*ymax 1.5*ymax]);
titext = sprintf('$\\hat{u_d}=%1.2f,\\ \\hat{\\delta}=%1.3f, \\ \\hat{\\lambda}=%1.3f$',ud_nd,delta_nd,lambda_nd);
title(titext);

%% Plot modeshape
[~,ind] = max(imag(c_chosen));
[z, phi] = case1.findmodeshape(an_c(:,ind));
figtitle = ["$\phi$", "$\phi_ z$", "$\phi_ {zz}$"];
xlab = {'$magnitude$','$angle$','$real$','$imag$'};
if (h > 6)
    blim = -6;
else
    blim = -h;
end
for i = 1:3
    fig(i) = figure('position',[0 0 1680 960]);
    plotvar = {abs(phi(:,i)),unwrap(angle(phi(:,i))),real(phi(:,i)),imag(phi(:,i))};
    for j = 1:4
        subplot(1,4,j);
        xline(0,'--','linewidth',1.5,'color','#606060');
        hold on;
        plot(plotvar{j},z,'-k.','linewidth',1,'markersize',10);
        if ~isnan(case1.zc)
            yline(case1.zc, '-r', 'linewidth', 1.5);
        end
        arr = -case1.getcut;
        for k = 2:length(arr)-1
            yline(arr(k), '--r', 'linewidth', 1);
        end
        hold off;
        xlabel(xlab{j});
        ylabel('$z$','rotation',0, 'HorizontalAlignment','right');
        ylim([blim 0]);
        grid on; box on;
    end
    sgtitle(figtitle(i), 'Interpreter', 'LaTeX');
end