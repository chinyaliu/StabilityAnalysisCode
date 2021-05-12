% close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray","D4"];
diff_meth = ["Schimd", "Trefethen"];
makeAB_meth = ["D4", "Schimd"];
solveGEPmeth = ["qr", "qz", "eig"];
% Set solver
method = [order(1), diff_meth(1), makeAB_meth(1)];
alg = solveGEPmeth(1);
% Inputs
Re = inf;
Fr2 = 2.25;
N = 600;
k = 1.25-2.7i;
h = 2*pi./real(k);
eps = 0.15;
inflec_pt = -0.74708299;
% c0 = 1./sqrt(kr*Fr2);
% zL = real(wZhang_ddm.g(c0)); 
zL = 0.74708299;
addvar = struct('zL1',zL(1),'eps',eps);
numberofDDM = 4;
f = wZhang_ddm.ddmtype(numberofDDM);
in_init = {N,k,h,Re,Fr2,method};
%% Run solver
N = 100:50:2000;
tic;
p1 = wZhang_ddm(in_init{:});
fig2 = figure('position',[50,50,1000,720]);
line([0,real(k)],[0,imag(k)],'color','r','linewidth',1.5);
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
xlim([0.1 0.26]);
ylim([-0.6 -0.2]);
hold on;
for i = 1:length(N)
    fprintf('N = %4d\n', N(i));
    p1.N = N(i);
    % test
    o = p1.solver(alg, 'all', f, addvar);
    ca = o/k;

    a = 1:length(ca); 
    crange = ((real(ca)-0.0012>-1e-5) & (real(ca)-1<=1e-5));
    dis = ((real(ca)-1).^2 +imag(ca).^2)>1e-5;
    aa = a(crange&dis);
    abch = isoutlier(imag(ca(aa)),'movmedian',20);
    aa = [a(dis&~crange) aa(abch)];
    o_c = o(aa);
    
    g(1) = scatter(real(o),imag(o),'ok');
%     g(2) = scatter(real(o_c),imag(o_c),'b','filled');
    titext = sprintf('$k=%.2f%+.2fi,\\ N=%4d$',real(k),imag(k), N(i));
    title(titext);
    
    F(i) = getframe(gcf);
    delete(g);
end
hold off;
toc;
%%
writerObj = VideoWriter('testk.mp4','MPEG-4');
writerObj.FrameRate = 10;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);