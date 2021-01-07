close all; clear all;% clc
%% Solver & Algorithm list
wZhangMethod = ["DDM", "Complex"];
order = ["Ray","D4"];
diff_meth = ["Schimd", "Trefethen"];
makeAB_meth = ["D4", "Schimd"];
solveGEPmeth = ["qr", "qz", "eig"];
%% Set solver
meth = wZhangMethod(1);
method = [order(1), diff_meth(1), makeAB_meth(1)];
alg = solveGEPmeth(1);
%% Inputs
do_balancing = 'n';
N = 600;
k = 0.6;
Re = inf;
Fr2 = 2.25;
h = 2*pi/k;
switch lower(meth)
    case 'ddm' % Additional parameters for DDM
        wZhang = @wZhang_ddm;
        numberofDDM = 4;
        eps = 0.01;
        c0 = 1./sqrt(k*Fr2);
        zL = wZhang_ddm.g(c0); 
%         zL = 0.74708299;
        f = wZhang_ddm.ddmtype(numberofDDM);
        in_init = {N,k,h,Re,Fr2,method};
        in_solver = {alg, do_balancing, f, struct('zL1',zL,'eps',eps)};
    case 'complex' % Additional parameters for complex
        wZhang = @wZhang_complex;
        delt = 0.02;
        in_init = {N,k,h,Re,Fr2,method,delt};
        in_solver = {alg, do_balancing};
    otherwise
        error('Method not defined');
end
%% Run solver
t1 = tic;
case1 = wZhang(in_init{:});
[o, an] = case1.solver(in_solver{:});
[z, phi] = case1.findmodeshape(an);
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
    plotvar = {abs(phi(:,i)),unwrap(angle(phi(:,i))),real(phi(:,i)),imag(phi(:,i))};
    for j = 1:4
        subplot(1,4,j);
        xline(0,'--','linewidth',1.5,'color','#606060');
        hold on;
        plot(plotvar{j},z,'-k.','linewidth',1,'markersize',10);
        if strcmpi(meth,'ddm')
            yline(case1.zc, '-r', 'linewidth', 1.5);
            arr = -case1.getarr;
            for k = 2:length(arr)-1
                yline(arr(k), '--r', 'linewidth', 1);
            end
        end
        hold off;
        set(gca,'fontsize',20);
        xlabel(xlab{j},'FontSize',30, 'Interpreter', 'LaTeX');
        ylabel('$z$','FontSize',30, 'Interpreter', 'LaTeX');
        ylim([blim 0]);
        grid on; box on;
    end
    sgtitle(figtitle(i),'FontSize',32, 'Interpreter', 'LaTeX');
end