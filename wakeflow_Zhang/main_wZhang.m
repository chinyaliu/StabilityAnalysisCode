close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray", "D4"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];
solveGEPmethod = ["qr", "qz", "eig"];
%% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1; % solveGEPmethod
do_balancing = 'n';
N = 400;
k = 0.3;
Re = inf;
Fr2 = 2.25;
h = 2*pi/k;
zL = 0.74708299;
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
t1 = tic;
case1 = wZhang_solver(N,k,h,Re,Fr2,method);
[o, an, cA, errGEP, dob] = case1.solver(zL, 'y', alg, do_balancing);
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
%         yline(case1.zc, '-.r', 'linewidth', 1.5);
        yline(-zL, '--b', 'linewidth', 1.5);
        hold off;
        set(gca,'fontsize',20);
        xlabel(xlab{j},'FontSize',30, 'Interpreter', 'LaTeX');
        ylabel('$z$','FontSize',30, 'Interpreter', 'LaTeX');
        ylim([blim 0]);
        grid on;
    end
    sgtitle(figtitle(i),'FontSize',32, 'Interpreter', 'LaTeX');
end
%% Save plots
folderpath = "fig_modeshape";
for i = 1:3
    figname = sprintf("k%02d_%d",10*k,i);
    namepath = folderpath + "\" + figname + ".png";
    exportgraphics(fig(i), namepath);
end