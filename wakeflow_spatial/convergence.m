% close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray", "Ray2"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];
solveGEPmethod = ["qr", "qz", "eig"];
% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 2; % solveGEPmethod
do_balancing = 'n';
N = 400;
o = 0.15;
Re = inf;
Fr2 = 2.25;
h = 20;
% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Convergence on grid number
Nt = 50:50:600;
t1 = tic;
case1 = wZhang_spatial(Nt(1),o,h,Re,Fr2,method);
kn = NaN(1,length(Nt));
for i = 1:length(Nt)
    case1.N = Nt(i);
    knspec = case1.solver(alg, do_balancing);
    kn(i) = eigchoose(o,knspec);
    fprintf('N = %3d, growth rate = %.8f\n', Nt(i), -imag(kn(i)));
end
toc(t1);
dkn = abs(imag(kn)-imag(kn(end)));
figure;
semilogy(Nt,dkn,'-ko');
xlabel('$N$');
ylabel('$\ | \ k_i(N_m) - k_i(N_{end})\ |$');
titext = sprintf('$o = %.3f,\\ k=%.2f%+.2fi$',o,real(kn(end)),imag(kn(end)));
title(titext,'FontSize',30);
grid on;
%% Convergence on truncation height
ht = 5:5:80;
t2 = tic;
case2 = wZhang_spatial(N,o,ht(1),Re,Fr2,method);
kh = NaN(1,length(ht));
for i = 1:length(ht)
    case2.h = ht(i);
    knspec = case2.solver(alg, do_balancing);
    kh(i) = eigchoose(o,knspec);
    fprintf('h = %3d, growth rate = %.8f\n', ht(i), -imag(kh(i)));
end
toc(t2);
dkh = abs(imag(kh)-imag(kh(end)));
figure;
semilogy(ht,dkh,'-o');
xlabel('$h$');
ylabel('$\ | \ k_i(h_m) - k_i(h_{end})\ |$');
titext = sprintf('$o = %.3f,\\ k=%.2f%+.2fi$',o,real(kh(end)),imag(kh(end)));
title(titext,'FontSize',30);
grid on;
%% Choose eigenvalue
function k_chosen = eigchoose(o,kspec)
c = o./kspec;
aa = (abs(imag(c))>1e-3) & (real(c)>1e-3) & (imag(kspec)<0);
if sum(aa)==0
    k_chosen = NaN;
else
    [~,dd] = max(real(c.*aa));
    k_chosen = kspec(dd);
end
end