close all; clear all;% clc
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [2,1,2]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
N = 500;
k = 0.3;
Re = 1000;
Fr2 = 2.25;
h = 6;
%% Run solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
[o, ~, cA, errGEP, db] = wZhang_solver(N,k,h,Re,Fr2,method,alg,do_balancing);