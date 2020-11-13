close all; clear all;% clc
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [2,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1; % solveGEPmethod
do_balancing = 'n';
N = 1000;
k = 3;
Re = inf;
Fr2 = 2.25;
h = @(k,m) 2*m*pi/k;
m = 0.5:0.1:4;
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
% %% change h
% for i = 1:length(m)
%     [o, ~, cA1(i), errGEP1(i), db1(i)] = wZhang_solver(N,k,h(k,m(i)),Re,Fr2,method,alg,do_balancing);
%     o1(i,:) = imag(o);
%     o1r(i,:) = real(o);
% %     fprintf('k = %.2f, growth rate = %.4f\n', k(i), o1(i));
%     fprintf('h = %.2f lambda, growth rate = %.8f\n', m(i), imag(o));
% %     [o, ~, cA2(i), errGEP2(i), db2(i)] = wZhang_solver(N,k(i),h,Re2,Fr2,method2,alg,do_balancing);
% %     o2(i,:) = imag(o);
% end
%% change N
m = 2;
N = 100:10:800;
for i = 1:length(N)
    [o, ~, cA1(i), errGEP1(i), db1(i)] = wZhang_solver(m*N(i),k,h(k,m),Re,Fr2,method,alg,do_balancing);
    o1(i,:) = imag(o);
    o1r(i,:) = real(o);
%     fprintf('k = %.2f, growth rate = %.4f\n', k(i), o1(i));
    fprintf('N = %4d, growth rate = %.8f\n', N(i), imag(o));
%     [o, ~, cA2(i), errGEP2(i), db2(i)] = wZhang_solver(N,k(i),h,Re2,Fr2,method2,alg,do_balancing);
%     o2(i,:) = imag(o);
end
plot(N,o1,'-ko');