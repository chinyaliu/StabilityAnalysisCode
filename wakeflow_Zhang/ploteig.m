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
do_balancing = 'n';
N = 400;
k = 1.5-1.5i;
Re = inf;
Fr2 = 2.25;
% Fr2 = 30.25;
h = 2*pi/real(k);
numberofDDM = 4;
eps = 0.01;
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
dis = (real(c)-0.5).^2 +imag(c).^2;
aa = (abs(imag(c))>5e-3) & (dis <= 0.25);
if sum(aa)==0
    o_chosen = NaN; c_chosen = NaN;
else
    [~,bb] = max(abs(imag(o.*aa)));
    o_chosen = o(bb); c_chosen = c(bb);
    an_c = an(:,bb);
end
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
% scatter(real(c_chosen),imag(c_chosen),'k','filled');
hold off;
xlabel('$c_r$','fontsize',30);
ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ymax = max(abs(imag(c)));
ylim([-1.5*ymax 1.5*ymax]);
titext = sprintf('$k=%.2f%+.2fi,\\ \\omega =%.2f%+.2fi$',real(k),imag(k),real(o_chosen),imag(o_chosen));
title(titext);
%% Plot oi_or
figure;
scatter(real(o),imag(o),'ok');
% plot(real(o),imag(o),'o');
hold on;
line([0,real(k)],[0,imag(k)],'color','r','linewidth',2);
% scatter(real(o_chosen),imag(o_chosen),'k','filled');
hold off;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
titext = sprintf('$k=%.2f%+.2fi,\\ \\omega =%.2f%+.2fi$',real(k),imag(k),real(o_chosen),imag(o_chosen));
title(titext);