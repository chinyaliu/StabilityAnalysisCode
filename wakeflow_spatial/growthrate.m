% close all; clear all; clc
%% Solver & Algorithm list
order = ["Ray", "Ray2"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];
solveGEPmethod = ["qr", "qz", "eig"];
% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 2; % solveGEPmethod
do_balancing = 'n';
N = 200;
o = linspace(0.001,1.3,30);
Re = inf;
Fr2 = 2.25;
h = 20;
% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
tic;
p1 = wZhang_spatial(N,o(1),h,Re,Fr2,method);
k = NaN(1,length(o));
c = NaN(1,length(o));
kall = cell(1,length(o));
for i = 1:length(o)
    p1.o = o(i);
    kall{i} = p1.solver(alg, do_balancing);
    % Choose eigenvalue
    c = o(i)./kall{i};
    aa = (abs(imag(c))>1e-3) & (real(c)>1e-3) & imag(kall{i})<0;
    if sum(aa)==0
        k(i) = NaN; c(i) = NaN;
    else
        [~,dd] = max(real(c.*aa));
        k(i) = kall{i}(dd); c(i) = c(dd);
    end
    fprintf('o = %.2f, growth rate = %.4f\n', o(i), -imag(k(i)));
end
toc;
%% Plot k_i v.s. k_r
fig1 = figure('position',[50,0,1000,720]);
plot(real(k),imag(k),'-ko');
% ylim([0 0.04]);
% xlim([0 4]);
xlabel('$k_r$','fontsize',30);
ylabel('$k_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% ax = gca;
% ax.YAxis.Exponent = -2;
grid on;
%% Plot growth rate v.s. omega
fig2 = figure('position',[50,0,1000,720]);
plot(o,imag(k),'-ko');
% ylim([0 0.04]);
% xlim([0 4]);
xlabel('$\omega $','fontsize',30);
ylabel('$k_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
grid on;
%% Plot whole eigenvalue spectrum
% fig4 = figure('position',[50,0,1000,720]);
figure;
for i = 1:length(o)
    ca = o(i)./kall{i};
    plot(real(ca),imag(ca),'.','Markersize',8);
    hold on
    set(gca,'fontsize',20);
    xlabel('$c_r$','fontsize',30);
    ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    ylim([-1.5 1.5]);
    xlim([-1 1]);
    grid on;
    titext = sprintf('o = %.3f',o(i));
    title(titext,'FontSize',30);
    F(i) = getframe(gcf);
end
hold off;
writerObj = VideoWriter('spec.mp4','MPEG-4');
writerObj.FrameRate = 10;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);