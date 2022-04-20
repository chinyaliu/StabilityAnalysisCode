clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,delta_nd,lambda_nd,c0,h,f,eps] = pars_Morland(1);
fprintf('u_d = %1.2f, delta = %1.3f, lambda = %1.3f\n',ud_nd,delta_nd,lambda_nd);

%% Run solver
t1 = tic;
case1 = wMorland(N,h,ud_nd,delta_nd,lambda_nd,method,bflow);
case1.setprop('k',case1.k-0.1i);
addvar = struct('zL1',case1.invbf(c0),'eps',eps);
[c,an,call,anall] = case1.solver_RE(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
toc(t1);
k = 2*pi./lambda_nd;
o = c.*k;
oall = call.*k;

%% Plot eigenvalue spectrum ci_cr
figure;
plot(real(call),imag(call),'ok');
hold on; 
yline(0,'linewidth',1.5,'color','#898989'); 
scatter(real(c),imag(c),'b','filled');
hold off;
xlabel('$c_r$','fontsize',30);
ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ymax = max(abs(imag(c)));
if ymax>0
    ylim([-1.5*ymax 1.5*ymax]);
end
titext = sprintf('$k=%.2f%+.2fi$',real(k),imag(k));
title(titext);

%% Plot oi_or
figure;
scatter(real(oall),imag(oall),'ok');
hold on;
yline(0,'linewidth',1.5,'color','#898989');
scatter(real(o),imag(o),'b','filled');
line([0,real(k)],[0,imag(k)],'color','r','linewidth',1);
hold off;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
titext = sprintf('$k=%.2f%+.2fi$',real(k),imag(k));
title(titext);

%% Plot modeshape
[z, phi] = case1.getprop('modeshape',an(:,1));
% phi = -phi;
figtitle = ["$\phi$", "$\phi_ z$", "$\phi_ {zz}$"];
xlab = {'$magnitude$','$angle$','$real$','$imag$'};
if (h > 6)
    blim = -6;
else
    blim = fix(-h);
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
            yline(-case1.zc, '-r', 'linewidth', 1.5);
        end
        arr = -case1.getprop('cut');
        for kk = 2:length(arr)-1
            yline(arr(kk), '--r', 'linewidth', 1);
        end
        hold off;
        xlabel(xlab{j});
        ylabel('$z$','rotation',0, 'HorizontalAlignment','right');
        ylim([blim 0]);
        grid on; box on;
    end
    sgtitle(figtitle(i), 'Interpreter', 'LaTeX');
end