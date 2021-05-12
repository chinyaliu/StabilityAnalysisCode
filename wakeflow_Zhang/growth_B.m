close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray","D4"];
diff_meth = ["Schimd", "Trefethen"];
makeAB_meth = ["D4", "Schimd"];
solveGEPmeth = ["qr", "qz", "eig"];
% Set solver
method = [order(1), diff_meth(1), makeAB_meth(1)];
alg = solveGEPmeth(1);
% Inputs
do_balancing = 'n';
Re = inf;
Fr2 = 2.25;
N = 400;
k = linspace(0.01,4,100);
h = 2*pi./real(k);
eps = 0.15;
inflec_pt = -0.74708299;
c0 = 1./sqrt(k*Fr2);
% zL = real(wZhang_ddm.g(c0)); 
zL = 0.74708299*ones(length(k),1);
cutz = NaN(1,length(k)+1);
cutz(1) = -inflec_pt;
addvar = struct('zL1',zL(1),'eps',eps);
numberofDDM = 4;
f = wZhang_ddm.ddmtype(numberofDDM);
in_init = {N,k,h,Re,Fr2,method};
%% Run solver
tic;
p1 = wZhang_ddm(in_init{:});
o = NaN(1,length(k)); z_c = NaN(1,length(k));
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i);
%     addvar.zL1 = zL(i);
    addvar.zL1 = cutz(i);
%     o(i) = p1.solver(alg, do_balancing, f, addvar);
    % test
    oall{i} = p1.solver(alg, do_balancing, f, addvar);
    o(i) = oall{i}(1);
    
    z_c(i) = p1.zc; 
    if isnan(z_c(i))
        cutz(i+1)=cutz(1);
    else
        cutz(i+1)=-p1.zc;
    end
    fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o(i)));
end
toc;
%%
figure; hold on;
for i = 1:length(k)
    ca = oall{i}/k(i);
    plot(real(oall{i}),imag(oall{i}),'.','Markersize',8);
%     dis = (real(ca)-0.5).^2 +imag(ca).^2;
%     aa = (imag(ca)>1e-2) & (dis <= 0.25);
    aa = imag(ca)>1e-2;
    if sum(aa) ~= 0
        [~,bb] = max(abs(imag(oall{i}.*aa)));
        o_c = oall{i}(bb); c_c = ca(bb);
%         plot(real(o_c),imag(o_c),'.','Markersize',8);
    end
%     plot(real(ch),imag(ch),'.','Markersize',8);
%     oh = ch*k(i);
end
hold off; box on;
%%
fig2 = figure('position',[50,0,1000,720]);
for i = 1:length(k)
    plot(real(oall{i}),imag(oall{i}),'.','Markersize',8);
    hold on
    set(gca,'fontsize',20);
    xlabel('$\tilde{\omega _r}$', 'Interpreter', 'LaTeX','fontsize',30);
    ylabel('$\tilde{\omega _i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    ylim([-2 0.5]);
    xlim([-0.5 2]);
    grid on;
    titext = sprintf('k = %.2f',k(i));
    title(titext,'FontSize',30);
    F(i) = getframe(gcf);
end
hold off;
writerObj = VideoWriter('test3.mp4','MPEG-4');
writerObj.FrameRate = 10;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);
%% Plot omega_i v.s. omega_r
fig3 = figure('position',[50,0,1000,720]);
plot(real(o),imag(o),'linewidth',2);
% ylim([0 0.04]);
% xlim([0 4]);
set(gca,'fontsize',20);
xlabel('$\tilde{\omega_r}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{\omega_i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% ax = gca;
% ax.YAxis.Exponent = -2;
grid on;
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
%% Plot z_c v.s. k
z_c(isnan(imag(o))) = NaN;
fig2 = figure('position',[50,0,1000,720]);
plot(k,z_c,'-ko','linewidth',2,'markersize',2);
hold on; yline(inflec_pt, '-.b', 'linewidth', 2); hold off;
set(gca,'fontsize',20);
xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
ylabel('$\tilde{z_c}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0);
grid on;
yt = sort([-3.5:0.5:0 inflec_pt]);
yticks(yt);
ind = find(yt==inflec_pt);
ax = gca;
ax.YTickLabel{ind} = ['\color{red}' ax.YTickLabel{ind}];
xlim([0 4]);