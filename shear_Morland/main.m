clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,delta_nd,lambda_nd,c0,h,f,epss,Re] = pars_Morland(5,'vis');
fprintf('u_d = %1.2f, delta = %1.3f, lambda = %1.3f\n',ud_nd,delta_nd,lambda_nd);

%% Run solver
t1 = tic;
case1 = wMorland(N,h,ud_nd,delta_nd,lambda_nd,method,bflow,Re);
addvar = struct('zL1',case1.invbf(c0),'eps',epss);
[o,an] = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
toc(t1);
k = case1.k;
c = o./k;

%% Choose eigenvalue
a = 1:length(c); 
crange = ((real(c)>-1e-5) & (real(c)-ud_nd<=1e-5));
cimag = imag(c)>-1;
aa = a(crange&cimag);
abci = isoutlier(imag(c(aa)),'movmedian',5);
abcr = isoutlier(real(c(aa)),'movmedian',20);
aa = [aa(abci|abcr) a(~crange)];
[~, ind] = sort(imag(c(aa)),'descend');
aa = aa(ind);
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
ymax = imag(c_chosen(1));
ymin = imag(c_chosen(end));
if ymax>0
    ylim([1.5*ymin 2*ymax]);
else
    ylim([-1 0.05]);
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

%% Plot simple mode shape
[z, phi] = case1.getprop('modeshape',an_c(:,1));
figtitle = ["$\phi$", "$\phi_ z$", "$\phi_ {zz}$"];
for i = 1:3
    f = figure('position',[0 0 480 960]);
    plotvar = {abs(phi(:,i)),unwrap(angle(phi(:,i))),real(phi(:,i)),imag(phi(:,i))};
    xline(0,'--','linewidth',1.5,'color','#606060','HandleVisibility','off');
    hold on;
    if ~isnan(case1.zc)
        yline(case1.zc, '-k','linewidth',1.5,'HandleVisibility','off');
    end
    arr = -case1.getprop('cut');
    for k = 2:length(arr)-1
        yline(arr(k), '--k','linewidth',1.5,'HandleVisibility','off');
    end
    plot(real(phi(:,i)),z,'b','DisplayName','real', 'linewidth', 2);
    plot(imag(phi(:,i)),z,'r','DisplayName','imag', 'linewidth', 2);
    hold off;
    ylabel('$z$','rotation',0, 'HorizontalAlignment','right');
    xticks(0);
    hold off;
    if (h > 6)
        blim = -6;
    else
        blim = fix(-h);
    end
    ylim([blim 0]);
    grid on; box on;
    legend('location','southeast');
    title(figtitle(i));
end

% %% Plot modeshape
% [z, phi] = case1.getprop('modeshape',an_c(:,1));
% % phi = -phi;
% figtitle = ["$\phi$", "$\phi_ z$", "$\phi_ {zz}$"];
% xlab = {'$magnitude$','$angle$','$real$','$imag$'};
% if (h > 6)
%     blim = -6;
% else
%     blim = fix(-h);
% end
% for i = 1:3
%     fig(i) = figure('position',[0 0 1680 960]);
%     plotvar = {abs(phi(:,i)),unwrap(angle(phi(:,i))),real(phi(:,i)),imag(phi(:,i))};
%     for j = 1:4
%         subplot(1,4,j);
%         xline(0,'--','linewidth',1.5,'color','#606060');
%         hold on;
%         plot(plotvar{j},z,'-k.','linewidth',1,'markersize',10);
%         if ~isnan(case1.zc)
%             yline(-case1.zc, '-r', 'linewidth', 1.5);
%         end
%         arr = -case1.getprop('cut');
%         for kk = 2:length(arr)-1
%             yline(arr(kk), '--r', 'linewidth', 1);
%         end
%         hold off;
%         xlabel(xlab{j});
%         ylabel('$z$','rotation',0, 'HorizontalAlignment','right');
%         ylim([blim 0]);
%         grid on; box on;
%     end
%     sgtitle(figtitle(i), 'Interpreter', 'LaTeX');
% end