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
% kr = 0.3;
% ki = -linspace(0,2,100);
kr = linspace(0.02,4,100);
ki = -2.5;
k = kr + 1i*ki;
h = ones(1,length(k)).*2*pi./kr;
eps = 0.01;
inflec_pt = -0.74708299;
% c0 = 1./sqrt(kr*Fr2);
% zL = real(wZhang_ddm.g(c0)); 
zL = 0.74708299*ones(length(k),1);
addvar = struct('zL1',zL(1),'eps',eps);
numberofDDM = 4;
f = wZhang_ddm.ddmtype(numberofDDM);
in_init = {N,k(1),h(1),Re,Fr2,method};
%% Run solver
tic;
p1 = wZhang_ddm(in_init{:});
oall{length(k),1} = [];
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i);
    addvar.zL1 = zL(i);
    % test
    o = p1.solver(alg, 'all', f, addvar);
    ca = o/k(i);

    a = 1:length(ca); 
    crange = ((real(ca)-0.0012>-1e-5) & (real(ca)-1<=1e-5));
    dis = ((real(ca)-1).^2 +imag(ca).^2)>1e-5;
    aa = a(crange&dis);
    abch = isoutlier(imag(ca(aa)),'movmedian',20);
    aa = [a(dis&~crange) aa(abch)];
    oall{i} = o(aa);

    fprintf('k = %.2f%+.2fi\n', real(k(i)), imag(k(i)));
end
toc;
%% Find k_i = 0
kk = linspace(0.01,4,100);
h = 2*pi./kk;
cutz = NaN(1,length(kk)+1);
cutz(1) = -inflec_pt;
o = NaN(1,length(kk));
p1.N = 400;
for i = 1:length(kk)
    p1.k = kk(i); p1.h = h(i);
%     addvar.zL1 = zL(i);
    addvar.zL1 = cutz(i);
    o(i) = p1.solver(alg, 'n', f, addvar);
    if isnan(p1.zc)
        cutz(i+1)=cutz(1);
    else
        cutz(i+1)=-p1.zc;
    end
    fprintf('k = %.2f\n', kk(i));
end
%%
fig2 = figure('position',[50,0,1000,720]);
% Read Dimas's results
imagedata = imread('dimas_fr15.bmp');
im2 = imbinarize(imagedata).*255;
im2(im2==0) = 200;
im2 = cast(im2,'uint8');
% fig1 = figure('position',[50,50,1280,720]);
% imagesc([0,1],[0.15,-0.1],im2); % Fr = 4.5
% imagesc([0,1],[0.15,-0.1],im2); % Fr = 3.5
% imagesc([0,1],[0.1,-0.4],im2); % Fr = 2.5
imagesc([0,1.4],[0.1,-0.6],im2); % Fr = 1.5
% imagesc([0,3.5],[0.2,-1.2],im2); % Fr = 0.5
set(gca,'YDir','normal');
hold on;
% plot
plot(real(o),imag(o),'k.','Markersize',6);
% xlim([0 1.5]);
% ylim([-1.2 .1]);
xlabel('$\tilde{\omega _r}$','fontsize',30);
ylabel('$\tilde{\omega _i}$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
for i = 1:length(k)
    plot(real(oall{i}),imag(oall{i}),'.','Markersize',8);
    h = line([0,real(k(i))],[0,imag(k(i))],'color','r','linewidth',2);
    grid on;
    titext = sprintf('$k=%.2f%+.2fi$',real(k(i)),imag(k(i)));
    title(titext,'FontSize',30);
    F(i) = getframe(gcf);
    delete(h);
end
hold off;
writerObj = VideoWriter('test25.mp4','MPEG-4');
writerObj.FrameRate = 10;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);