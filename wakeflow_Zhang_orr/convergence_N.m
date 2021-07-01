clear all;% close all; clc
%% Solver & Algorithm list
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
do_balancing = 'n';
eig_spectrum = 'max';
N = 20:20:1200;
% k = [0.01 0.2 0.4 0.6];
% k = [0.87 1 2 3 4];
k = 0.2;
Re = 1000;
Fr2 = 2.25;
h = @(k) 2*pi/k*4.79;
c0 = 1./sqrt(k*Fr2);
% zL = wZhang_ddm.g(c0); 
zL = 0.74708299*ones(length(k),1);
% DDM numbers
numberofDDM = 4;
eps = 0.2;
f = wZhang_ddm.ddmtype(numberofDDM);
in_init = {N(1),k(1),h(k(1)),Re,Fr2};
addvar = struct('zL1',zL(1),'eps',eps);
%% Run solver
tic;
case1 = wZhang_ddm(in_init{:});
case1.numMeth(method);
for j = 1:length(k)
    fprintf('k = %.2f\n',k(j));
    case1.k = k(j); case1.h = h(k(j)); 
    addvar.zL1 = zL(j);
    for i = 1:length(N)
        case1.N = N(i);
        o = case1.solver(alg, do_balancing, eig_spectrum, f, addvar);
        oi(i) = imag(o(1));
        fprintf('N = %3d, growth rate = %.8f\n', N(i), oi(i));
    end
    oiall(j,:) = oi;
    doi{j} = abs(diff(oi));
end
toc;
%% Plot figure
dall = abs(oiall - oiall(:,end).*ones(size(oiall,1),size(oiall,2)))./oiall(:,end);
% fig1 = figure('position',[50,0,1080,840]);
figure;
nam = [sprintf('$k = %.2f, ',k(1)) '\' sprintf('omega = %0.9e$',oiall(1,end))];
semilogy(N, dall(1,:), '-o', 'Displayname', nam);
hold on;
for i = 2:length(k)
    nam = [sprintf('$k = %.2f, ',k(i)) '\' sprintf('omega = %0.9e$',oiall(i,end))];
    semilogy(N, dall(i,:), '-o', 'Displayname', nam);
end
hold off;
xlabel('$N$');
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{end})\ |/\omega_i(N_{end})$');
legend('location','northeast');
grid on;
% exportgraphics(fig1, 'fig_convergence\branch2.png');