% close all; clear all;% clc
%% Input list
order = ["Ray","D4"];
diff_meth = ["Schimd", "Trefethen"];
makeAB_meth = ["D4", "Schimd"];
solveGEPmeth = ["qr", "qz", "eig"];
% Set solver
method = [order(1), diff_meth(1), makeAB_meth(1)];
alg = solveGEPmeth(1);
% Inputs
do_balancing = 'all';
N = 1000;
k = 0.2;
Re = inf;
Fr2 = 2.25;
% Fr2 = 30.25;
h = 2*pi/real(k);
numberofDDM = 4;
eps = 0.15;
c0 = 1./sqrt(k*Fr2);
% zL = wZhang_ddm.g(c0); 
zL = 0.74708299;
f = wZhang_ddm.ddmtype(numberofDDM);
%% Run solver
t1 = tic;
case1 = wZhang_ddm(N,k,h,Re,Fr2,method);
[o, an] = case1.solver(alg, do_balancing, f, struct('zL1',zL,'eps',eps));
toc(t1);
c = o/k;
%% Choose eigenvalue
% dis = (real(c)-0.5).^2 +imag(c).^2;
% aa = (abs(imag(c))>5e-3) & (dis <= 0.25) & (real(c)>0);
% oa = o(aa); ca = c(aa);
% if sum(aa)==0
%     o_chosen = NaN; c_chosen = NaN;
% else
%     [~,bb] = max(imag(oa));
%     o_chosen = oa(bb); c_chosen = ca(bb);
%     an_c = an(:,bb);
% end
a = 1:length(c); 
crange = ((real(c)-0.0012>-1e-5) & (real(c)-1<=1e-5));
dis = ((real(c)-1).^2 +imag(c).^2)>1e-5;
aa = a(crange&dis);
abch = isoutlier(imag(c(aa)),'movmedian',20);
aa = [a(dis&~crange) aa(abch)];
o_chosen = o(aa); c_chosen = c(aa);
an_c = an(:,aa);
%% Plot eigenvalue spectrum ci_cr
figure;
h1 = viscircles([0.5,0], 0.5,'color','w','LineStyle','none');
xd = h1.Children(1).XData(1:end-1);
yd = h1.Children(1).YData(1:end-1);
% set(gca,'Color','#e7e7e7');
hold on;
% fill(xd, yd, 'w','LineStyle','--','edgecolor', '#898989');
% fill([0 0 1 1], [0 -1 -1 0],[0.9 0.9 0.9],'LineStyle','-','edgecolor','#e7e7e7');
plot(real(c),imag(c),'ok');
line([0,1],[0,0],'color','r','linewidth',2.5);
scatter(real(c_chosen),imag(c_chosen),'b','filled');
hold off;
xlabel('$c_r$','fontsize',30);
ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ymax = max(abs(imag(c)));
ylim([-1.5*ymax 1.5*ymax]);
titext = sprintf('$k=%.2f%+.2fi$',real(k),imag(k));
% titext = sprintf('$k=%.2f%+.2fi,\\ \\omega =%.2f%+.2fi$',real(k),imag(k),real(o_chosen),imag(o_chosen));
title(titext);
box on;
%% Plot oi_or
figure;
scatter(real(o),imag(o),'ok');
% plot(real(o),imag(o),'o');
hold on;
line([0,real(k)],[0,imag(k)],'color','r','linewidth',2);
scatter(real(o_chosen),imag(o_chosen),'b','filled');
hold off;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
titext = sprintf('$k=%.2f%+.2fi$',real(k),imag(k));
% titext = sprintf('$k=%.2f%+.2fi,\\ \\omega =%.2f%+.2fi$',real(k),imag(k),real(o_chosen),imag(o_chosen));
title(titext);
%% Plot modeshape
xlab = {'$magnitude$','$angle$','$real$','$imag$'};
if (h > 6)
    blim = -6;
else
    blim = fix(-h);
end
for i = 1:size(an_c,2)
    [z, phi] = case1.findmodeshape(an_c(:,i));
    fig(i) = figure('position',[0 0 1680 960]);
    plotvar = {abs(phi(:,1)),unwrap(angle(phi(:,1))),real(phi(:,1)),imag(phi(:,1))};
    for j = 1:4
        subplot(1,4,j);
        xline(0,'--','linewidth',1.5,'color','#606060');
        hold on;
        plot(plotvar{j},z,'-k.','linewidth',1,'markersize',10);
        yline(case1.zc, '-r', 'linewidth', 1.5);
        arr = -case1.getcut;
        for k = 2:length(arr)-1
            yline(arr(k), '--r', 'linewidth', 1);
        end
        hold off;
        set(gca,'fontsize',20);
        xlabel(xlab{j},'FontSize',30, 'Interpreter', 'LaTeX');
        ylabel('$z$','FontSize',30, 'Interpreter', 'LaTeX');
        ylim([blim 0]);
        grid on; box on;
    end
    titext = ['$\phi ,\ \omega = ' sprintf('%.3f%+.3fi$', real(o_chosen(i)),imag(o_chosen(i)))];
    sgtitle(titext,'FontSize',32, 'Interpreter', 'LaTeX');
end