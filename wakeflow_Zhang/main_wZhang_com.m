close all; clear all;% clc
%% Input list
order = ["Ray", "D4"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];
solveGEPmethod = ["qr", "qz", "eig"];
% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1; % solveGEPmethod
do_balancing = 'n';
N = 800;
k = 2.5;
Re = inf;
Fr2 = 2.25;
delt = 0.01;
h = 2*pi/k;
% if k > pi/3
%     h = 6;
% end
zL = 0.74708299;
% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
t1 = tic;
case1 = wZhang_complex(N,k,h,Re,Fr2,method,delt);
% [o, an, cA, errGEP, dob] = case1.solver(zL, 'y', alg, do_balancing);
[o, an] = case1.solver(alg, do_balancing);
[z, phi] = case1.findmodeshape(an(:,1));
toc(t1);
c = o/k;
%% Plot eigenvalue spectrum ci_cr
figure;
h1 = viscircles([0.5,0], 0.5,'color','w','LineStyle','none');
xd = h1.Children(1).XData(1:end-1);
yd = h1.Children(1).YData(1:end-1);
set(gca,'Color','#e7e7e7');
hold on;
fill(xd, yd, 'w','LineStyle','--','edgecolor', '#898989');
fill([0 0 1 1], [0 -1 -1 0],[0.9 0.9 0.9],'LineStyle','-','edgecolor','#e7e7e7');
plot(real(c),imag(c),'ok');
Ub = case1.baseflow();
cb = Ub(:,1);
plot(real(cb),imag(cb),'r','linewidth',2);
% scatter(real(c_chosen),imag(c_chosen),'k','filled');
hold off;
xlabel('$c_r$','fontsize',30);
ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ymax = max(abs(imag(c)));
ylim([-1.5*ymax 1.5*ymax]);
% titext = sprintf('$k=%.1f%+.1fi,\\ \\omega =%.2f%+.2fi$',real(k),imag(k),real(o_chosen),imag(o_chosen));
% title(titext);
%% Plot oi_or
figure;
scatter(real(o),imag(o),'ok');
% plot(real(o),imag(o),'o');
hold on;
ob = cb*k;
plot(real(ob),imag(ob),'r','linewidth',2);
% scatter(real(o_chosen),imag(o_chosen),'k','filled');
hold off;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% titext = sprintf('$k=%.1f%+.1fi,\\ \\omega =%.2f%+.2fi$',real(k),imag(k),real(o_chosen),imag(o_chosen));
% title(titext);