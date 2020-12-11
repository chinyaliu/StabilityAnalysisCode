close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray"];
diff_method = ["Schimd", "Trefethen"];
solveGEPmethod = ["qr", "qz", "eig"];
%% Inputs
solver = [1,1]; % [order, diff_method]
algorithm = 1; % solveGEPmethod
do_balancing = 'n';
N = 800;
k = 0.01;
Re = inf;
Fr2 = 2.25;
h = 2*pi/k;
% h = 6;
zL = 0.74708299;
eps = 0.15;
%% Set solver
method = [order(solver(1)), diff_method(solver(2))];
alg = solveGEPmethod(algorithm);
%% Run solver
t1 = tic;
case1 = wZhang_block(N,k,h,Re,Fr2,method);
[o, an] = case1.solver(zL, 'y', alg, do_balancing,eps);
toc(t1);
%% Plot
figtitle = ["$\phi$", "$\phi_ z$", "$\phi_ {zz}$"];
xlab = {'$magnitude$','$angle$','$real$','$imag$'};
if (h > 6)
    blim = -6;
else
    blim = fix(-h);
end
for i = 1:3
    fig(i) = figure('position',[0 0 1680 960]);
    plotvar = {abs(case1.phi(:,i)),unwrap(angle(case1.phi(:,i))),real(case1.phi(:,i)),imag(case1.phi(:,i))};
    for j = 1:4
        subplot(1,4,j);
        plot(plotvar{j},case1.z,'-k.','linewidth',1,'markersize',10);
        hold on;
        yline(case1.zc, '-.r', 'linewidth', 1.5);
        yline(case1.zc-case1.cL, '--r', 'linewidth', 1);
        yline(case1.zc+case1.cL, '--r', 'linewidth', 1);
        xline(0,'--b','linewidth',1.5);
        hold off;
        set(gca,'fontsize',20);
        xlabel(xlab{j},'FontSize',30, 'Interpreter', 'LaTeX');
        ylabel('$z$','FontSize',30, 'Interpreter', 'LaTeX');
        ylim([blim 0]);
        grid on;
    end
    sgtitle(figtitle(i),'FontSize',32, 'Interpreter', 'LaTeX');
end