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
do_balancing = 'n';
Re = inf;
% Fr2 = 2.25;
Fr2 = inf;
N = 400;
kr = linspace(0.01,4,400);
h = 2*pi./kr;
eps = 0.15;
inflec_pt = -0.74708299;
c0 = 1./sqrt(kr*Fr2);
% zL = real(wZhang_ddm.g(c0)); 
zL = 0.74708299*ones(length(kr),1);
cutz = NaN(1,length(kr)+1);
cutz(1) = -inflec_pt;
addvar = struct('zL1',zL(1),'eps',eps);
numberofDDM = 4;
f = wZhang_ddm.ddmtype(numberofDDM);
in_init = {N,kr,h,Re,Fr2,method};
ki = [-0.2 -0.5 -0.8 -1];
%% Run solver
tic;
figure; hold on;
p1 = wZhang_ddm(in_init{:});
o = NaN(1,length(kr)); z_c = NaN(1,length(kr));
for j = 1:length(ki)
    k = kr+ki(j)*1i;
    od = NaN(1,length(k));
    for i = 1:length(k)
        p1.k = k(i); p1.h = h(i);
        addvar.zL1 = zL(i);
    %     addvar.zL1 = cutz(i);
    %     o(i) = p1.solver(alg, do_balancing, f, addvar);
        % test
        oall = p1.solver(alg, do_balancing, f, addvar);
        o(i) = oall(1);

        ca = oall/k(i);
        aa = abs(imag(ca))>1e-2 & real(oall)>0;
        suma = sum(aa);
        if suma ~= 0
            [~,bb] = max(abs(imag(ca.*aa)));
            od(i) = oall(bb);
%             oo = oall(aa);
%             od = [od; oo];
        end
        fprintf('k = %.2f, growth rate = %.4f\n', k(i), imag(o(i)));
    end
    textk = sprintf('$k_i=%+.3fi$',imag(k(i)));
    plot(real(od),imag(od),'.','Markersize',6,'DisplayName',textk);
end
legend('interpreter','latex','location','southwest');
toc;
hold off; box on;
%% Plot real k
k = kr;
od = []; o = NaN(1,length(k));
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i);
    addvar.zL1 = zL(i);
    % test
    oall = p1.solver(alg, do_balancing, f, addvar);

    ca = oall/k(i);
    dis = (real(ca)-0.5).^2 +imag(ca).^2;
    aa = (imag(ca)>1e-5) & (dis <= 0.25) & (real(oall)>0);
    if sum(aa)~=0
        [~,bb] = max(abs(imag(oall.*aa)));
        o(i) = oall(bb);
    end
    fprintf('k = %.2f\n', k(i));
end
textk = sprintf('$k_i=%+.3fi$',-imag(k(1)));
hold on;
plot(real(o),imag(o),'k.','Markersize',6,'DisplayName',textk);
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');