% close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray", "D4"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd"];
solveGEPmethod = ["qr", "qz", "eig"];
% Inputs
solver = [1,1,1]; % [order, diff_method, constructAB_method]
algorithm = 2; % solveGEPmethod
do_balancing = 'n';
N = 400;
o = 0.08+0.03i;
Re = inf;
Fr2 = 2.25;
h = 10;
% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Run solver
t1 = tic;
case1 = wZhang_spatial(N,o,h,Re,Fr2,method);
[k, an] = case1.solver(alg, do_balancing);
toc(t1);
%% Choose eigenvalue
aa = (abs(imag(k))>1e-2) & (abs(real(k))>1e-2);
bb = (abs(imag(k))<10) & (abs(real(k))<10);
 k1 = k(aa&bb);
 figure;plot(real(k1),imag(k1),'o');