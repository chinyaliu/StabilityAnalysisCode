close all; clear all;% clc
tic;
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [2,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
Re = logspace(3,9,20);
N = 300;
k = 0.1:0.005:1.1;
%% Run solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
snum = 1;
fnum = 0;
for ki = 1:length(k)
    kn = k(ki);
    [o,~,~,~,~] = poiseuille_solver(N,kn,Re(1),method,alg,do_balancing);
    oi(1) = imag(o);
    for i = 2:length(Re)
        [o,~,~,~,~] = poiseuille_solver(N,kn,Re(i),method,alg,do_balancing);
        oi(i) = imag(o);
    end
    cs = (oi(2:end).*oi(1:end-1))<0;
    pos = find(cs);
    for i = 1:length(pos)
        j = pos(i);
        x1 = Re(j+1);
        x0 = Re(j);
        s1 = oi(j+1);
        s0 = oi(j);
        for count = 1:1000
            xn = x1-s1*(x1-x0)/(s1-s0);
            [sn,~,~,~,~] = poiseuille_solver(N,kn,xn,method,alg,do_balancing);
            sn = imag(sn);
            if (abs(sn) <= 1e-7)
                xs0(snum) = xn;
                ys0(snum) = kn;
                ss0(snum) = sn;
                snum = snum+1;
                flag = 1;
                break;
            end
            if sn*s0 > 0
                x0 = xn;
                s0 = sn;
            else
                x1 = xn;
                s1 = sn;
            end
            flag = 0;
        end
        if (flag == 0)
            fnum = fnum+1;
            xf(fnum) = xn;
            yf(fnum) = kn;
            sf(fnum) = sn;
        end
    end
end
toc;
fig1 = figure('position',[50,50,720,540]);
semilogx(xs0,ys0,'ko');
xlim([1e+3 1e+9]);
ylim([0.1 1.1]);
xlabel('Re','fontsize',14');
ylabel('k','fontsize',14');
grid on;
[a,b] = min(xs0);
tx = sprintf('Re_{cr} = %.3f',a);
text(2*1e+3,0.4,tx,'fontsize',14);
% print(fig1,'msc_poiseuille','-r800','-dpng');
print(fig1,'msc_long','-r800','-dpng');
[Re0,I] = sort(xs0);
k0 = ys0(I);
save('neutral_pts.mat','Re0','k0');