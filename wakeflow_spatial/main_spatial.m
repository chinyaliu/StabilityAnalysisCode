% close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray", "Ray2"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];
solveGEPmethod = ["qr", "qz", "eig"];
% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 2; % solveGEPmethod
do_balancing = 'n';
N = 400;
o = 0.5;
Re = inf;
Fr2 = 2.25;
h = 60;
% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
t1 = tic;
case1 = wZhang_spatial(N,o,h,Re,Fr2,method);
[k, an] = case1.solver(alg, do_balancing);
toc(t1);
%% Choose eigenvalue
c = o./k;
aa = (abs(imag(c))>1e-3) & (real(c)>1e-3) & (imag(k)<0);
if sum(aa)==0
    k_chosen = NaN; c_chosen = NaN;
else
    [~,dd] = max(real(c.*aa));
    k_chosen = k(dd); c_chosen = c(dd);
    an_c = an(:,dd);
end
%% Plot ci vs cr
c = o./k;
figure;
plot(real(c),imag(c),'ok');
hold on;
xline(0,'Color','#898989','linewidth',2);
yline(0,'Color','#898989','linewidth',2);
scatter(real(c_chosen),imag(c_chosen),'r','filled');
hold off;
xlabel('$c_r$','fontsize',30);
ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
%% Plot ki vs kr
figure;plot(real(k),imag(k),'ok');
hold on;
xline(0,'Color','#898989','linewidth',2);
yline(0,'Color','#898989','linewidth',2);
scatter(real(k_chosen),imag(k_chosen),'r','filled');
hold off;
xlabel('$k_r$','fontsize',30);
ylabel('$k_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
%% Plot modeshape
phi = case1.modeshape(an_c);
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
        plot(plotvar{j},case1.z,'-k.','linewidth',1,'markersize',10);
        hold on;
%         yline(case1.zc, '-.r', 'linewidth', 1.5);
%         yline(-zL, '--b', 'linewidth', 1.5);
        xline(0,'--b','linewidth',1.5);
        hold off;
        set(gca,'fontsize',20);
        xlabel(xlab{j},'FontSize',30, 'Interpreter', 'LaTeX');
        ylabel('$z$','FontSize',30, 'Interpreter', 'LaTeX');
%         ylim([blim 0]);
        grid on;
    end
    sgtitle(figtitle(i),'FontSize',32);
end