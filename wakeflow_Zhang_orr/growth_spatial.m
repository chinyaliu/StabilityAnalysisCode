% close all; clear all;% clc
%% Solver & Algorithm list
diff_meth = ["Schimd", "Trefethen"];
solveGEPmeth = ["qr", "qz", "eig"];
% Set solver
method = diff_meth(1);
alg = solveGEPmeth(1);
% Inputs
Re = 1e3;
Fr2 = 1.5^2;
N = 600;
% Cusp method from k_i
kr = linspace(0.01,4,100);
ki = [-0.5 -1 -1.5 -1.7];
% ki = -2.5;
% Cusp method from k_i
% kr = [0.5 0.55 0.6 0.65 0.7];
% ki = linspace(-2,-0.01,100);
h = 2*pi./kr;
eps = 0.01;
inflec_pt = -0.74708299;
% c0 = 1./sqrt(kr*Fr2);
% zL = real(wZhang_ddm.g(c0)); 
zL = 0.74708299*ones(length(kr),1);
cutz = NaN(1,length(kr)+1);
cutz(1) = -inflec_pt;
addvar = struct('zL1',zL(1),'eps',eps);
numberofDDM = 4;
f = wZhang_ddm.ddmtype(numberofDDM);
in_init = {N,kr,h,Re,Fr2};
p1 = wZhang_ddm(in_init{:});
p1.numMeth(method);
%% Read Dimas's results
imagedata = imread('dimas_fr15.bmp');
im2 = imbinarize(imagedata).*255;
im2(im2==0) = 200;
im2 = cast(im2,'uint8');
fig1 = figure('position',[50,50,1280,720]);
% imagesc([0,1],[0.15,-0.1],im2); % Fr = 4.5
% imagesc([0,1],[0.15,-0.1],im2); % Fr = 3.5
% imagesc([0,1],[0.1,-0.4],im2); % Fr = 2.5
imagesc([0,1.4],[0.1,-0.6],im2); % Fr = 1.5
% imagesc([0,3.5],[0.2,-1.2],im2); % Fr = 0.5
set(gca,'YDir','normal');
hold on;
%% Plot real k
kk = linspace(0.01,4,100);
hk = 2*pi./kk;
cutz = NaN(1,length(kk)+1);
cutz(1) = -inflec_pt;
o = NaN(1,length(kk));
for i = 1:length(kk)
    p1.k = kk(i); p1.h = hk(i);
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
% textk = sprintf('$k_i=%+.3fi$',-imag(k(1)));
textk = sprintf('$k_r=%.3fi$',imag(kk(1)));
plot(real(o),imag(o),'k.','Markersize',8,'DisplayName',textk);
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
%% Run solver
tic;
% figure; hold on;
o = NaN(1,length(ki));
% o = NaN(1,length(kr)); z_c = NaN(1,length(kr));
for j = 1:length(ki)
    k = kr+ki(j)*1i;
% for j = 1:length(kr)
%     k = kr(j)+ki*1i;  
%     p1.h = h(j); addvar.zL1 = zL(j);
%     od = NaN(1,length(k));
    od = [];
    for i = 1:length(k)
%         p1.k = k(i);
        p1.k = k(i); p1.h = h(i);
        addvar.zL1 = zL(i);
        % test
        oall = p1.solver(alg, 'all', f, addvar);
        o(i) = oall(1);
        od = [od; oall];
%         ca = oall/k(i);
        
%         % Criteria 1
%         aa = abs(imag(ca))>5e-3;% & real(oall)>0;
%         suma = sum(aa);
%         if suma ~= 0
% %             [~,bb] = max(abs(imag(ca.*aa)));
% %             od(i) = oall(bb);
%             oo = oall(aa);
%             od = [od; oo];
%         end
        
%             % Criteria 2
%             a = 1:length(ca); 
%             crange = ((real(ca)-0.0012>-1e-5) & (real(ca)-1<=1e-5));
%             dis = ((real(ca)-1).^2 +imag(ca).^2)>1e-5;
%             aa = a(crange&dis);
%             abch = isoutlier(imag(ca(aa)),'movmedian',20);
%             aa = [a(dis&~crange) aa(abch)];
%             od = [od; oall(aa)];

            fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o(i)));
        end
    textk = sprintf('$k_i=%+.3fi$',imag(k(i)));
%     textk = sprintf('$k_r=%.3fi$',real(k(i)));
    plot(real(od),imag(od),'.','Markersize',8,'DisplayName',textk);
end
[~, objh] = legend('interpreter','latex','location','southeastoutside');
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 20);
toc;
hold off; box on;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');