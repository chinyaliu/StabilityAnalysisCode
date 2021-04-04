close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray", "D4"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];%, "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
N = 400;
k = linspace(0.01,4,100);
Re = inf;
Fr2 = 2.25;
delt = 0;
inflec_pt = -0.74708299;
% h = 6*ones(1,length(k));
% h(k<pi/3) = 2*pi./k(k<pi/3);
h = 2*pi./k;
% h = 1*ones(length(k),1);
% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
tic;
p1 = wZhang_complex(N,1,1,Re,Fr2,method,delt);
o = NaN(1,length(k)); z_c = NaN(1,length(k));
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i);
%     [o(i), an] = p1.solver(cutz(i), 'y', alg, do_balancing);
%     z(:,i) = p1.z; phi{i} = p1.phi; 
    op = p1.solver(alg, do_balancing);
    call{i} = op/k(i);
    o(i) = op(1);
    z_c(i) = -p1.zc;
    fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o(i)));
end
toc;
%% Plot growth rate v.s. k
fig1 = figure('position',[50,0,1000,720]);
plot(k,imag(o),'linewidth',2);
ylim([0 0.04]);
xlim([0 4]);
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{\omega_i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
ax.YAxis.Exponent = -2;
grid on;
%% Plot c_i vs c_r
fig2 = figure('position',[50,0,1000,720]);
for i = 1:length(k)
    plot(real(call{i}),imag(call{i}),'.','Markersize',8);
    hold on
    set(gca,'fontsize',20);
    xlabel('$\tilde{c_r}$', 'Interpreter', 'LaTeX','fontsize',30);
    ylabel('$\tilde{c_i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    ylim([-0.12 0.12]);
    xlim([-2 4]);
    grid on;
    titext = sprintf('k = %.2f',k(i));
    title(titext,'FontSize',30);
    F(i) = getframe(gcf);
end
hold off;
% nam = sprintf('comp1_h%03d.mp4',100*delt);
% vname = ['mov_eigenvalue\' nam];
% writerObj = VideoWriter(vname,'MPEG-4');
% writerObj.FrameRate = 10;
% writerObj.Quality = 100;
% open(writerObj);
% for i = 1:length(F)
%     frame = F(i);
%     writeVideo(writerObj, frame);
% end
% close(writerObj);
%% Plot z_c v.s. k
fig2 = figure('position',[50,0,1000,720]);
plot(k,z_c,'linewidth',2);
hold on; yline(inflec_pt, '-.r', 'linewidth', 2); hold off;
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{z_c}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0);
grid on;
yt = sort([-3.5:0.5:0 inflec_pt]);
yticks(yt);
% ind = find(yt==inflec_pt);
% ax = gca;
% ax.YTickLabel{ind} = ['\color{red}' ax.YTickLabel{ind}];
xlim([0 4]);
% ylim([-3.6 0]);
%% Save data & figures
% save('diffk','phi','z','N','k','cutz','o','z_c');
% exportgraphics(fig1, 'fig_growthrate\omega_i.png');
% exportgraphics(fig2, 'fig_growthrate\criticalheight.png');