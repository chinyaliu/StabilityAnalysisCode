clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,H,k,Fr2,Re,eps,c0,h,f] = pars_wake;

%% Run solver
t1 = tic;
case1 = wSubmerged(N,H,k,h,Re,Fr2);
case1.numMeth(method);
zL = 0.74708299+H;
[o, an] = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',zL,'eps',eps));
o = o(real(o)>-50);
c = o/k;
toc(t1);

%% Choose eigenvalue
a = 1:length(c); 
crange = (real(c)>0) & (real(c)<0.5) & (imag(c)<0);
aa = a(crange);
abch = isoutlier(imag(c(aa)),'movmedian',5);
aa = [a(imag(c)>0) a(real(c)<0) aa(abch)];
% aa = aa(abch);
o_chosen = o(aa); c_chosen = c(aa);
an_c = an(:,aa);

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
if ymax>0
    ylim([-1.5*ymax 1.5*ymax]);
end
titext = sprintf('$k=%.2f%+.2fi$',real(k),imag(k));
title(titext);

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

%% Plot modeshape
[z, phi] = case1.findmodeshape(an(:,1));
figtitle = ["$\phi$", "$\phi_ z$", "$\phi_ {zz}$"];
xlab = {'$magnitude$','$angle$','$real$','$imag$'};
if (h > 6)
    blim = -6-H;
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
            yline(case1.zc-2*H, '-r', 'linewidth', 1.5);
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