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
N = 400;
% Cusp method from k_i
kr = [0.3 1.5 2];
ki = -linspace(0.01,1.5,100);
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
figure;
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
cutz = NaN(1,length(kk));
ok = NaN(1,length(kk));
for i = 1:length(kk)
    p1.k = kk(i); p1.h = hk(i);
    ok(i) = p1.solver(alg, 'n', 'max', f, addvar);
    if isnan(p1.zc)
        cutz(i)=-inflec_pt;
    else
        cutz(i)=-p1.zc;
    end
    addvar.zL1 = cutz(i);
    fprintf('k = %.2f\n', kk(i));
end
% textk = sprintf('$k_i=%+.3fi$',-imag(k(1)));
textk = sprintf('$k_r=%.2fi$',imag(kk(1)));
plot(real(ok),imag(ok),'k.','Markersize',6,'DisplayName',textk);
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
for i = 1:length(kr)
    [~,ind(i)] = min(abs(kk-kr(i)));
end
kr = kk(ind);
cutz = cutz(ind);
h = 2*pi./kr;
or = ok(ind);
%% Run solver
col = get(gca,'colororder');
tic;
o_sel = cell(length(kr),1);
for j = 1:length(kr)
    k = kr(j)+ki*1i;  
    p1.h = h(j); addvar.zL1 = cutz(j);
    od = [];
    oref = or(j);
    o = NaN(1,length(k));
    for i = 1:length(k)
        p1.k = k(i);
        oall = p1.solver(alg, 'n', 'all', f, addvar);
        od = [od; oall];
        
        [~,a] = min(abs(oall-oref));
        o(i) = oall(a);
        oref = o(i);
        
        fprintf('k_i = %.2f, growth rate = %.4f\n', imag(k(i)), imag(o(i)));
    end
    textk = sprintf('$k_r=%.2f$',round(real(k(i)),2));
    plot(real(od),imag(od),'.','Markersize',6,'HandleVisibility','off','color',col(j,:));
    plot(real(o),imag(o),'o','Markersize',4,'DisplayName',textk,'color',col(j,:));
    o_sel{j} = o;
end
[~, objh] = legend('interpreter','latex','location','southeastoutside');
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 6);
toc;
hold off; box on;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
%% Plot selected eigenvalue
figure;
imagesc([0,1.4],[0.1,-0.6],im2); % Fr = 1.5
set(gca,'YDir','normal');
hold on;
plot(real(ok),imag(ok),'k.','Markersize',6,'HandleVisibility','off');
for i = 1:length(kr)
    textk = sprintf('$k_r=%.2f$',round(real(k(i)),2));
    plot(real(o_sel{i}),imag(o_sel{i}),'.','Markersize',8,'DisplayName',textk,'color',col(i,:));
end
[~, objh] = legend('interpreter','latex','location','southeastoutside');
objhl = findobj(objh, 'type', 'line'); 
set(objhl, 'Markersize', 6);
toc;
hold off; box on;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');