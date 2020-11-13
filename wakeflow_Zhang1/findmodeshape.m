close all; clear all;% clc
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [2,2,2]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
N = 1001;
k = 0.1:0.01:3.5;
Re = inf;
Fr2 = 2.25;
% h = @(k) 2*pi/k;
h = @(k) 6;
%% Set solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Find critical height
if (strcmpi(method(1),'d4') && strcmpi(method(2),'schimd') && strcmpi(method(3),'d4'))
    M = N-2;
    [zeta,D] = Dcheb(M,1,'d4');
else
    M = N;
    [zeta,D] = Dcheb(N,1,'n');
end
if (strcmpi(method(2),'trefethen'))
    D = repmat(eye(N+1),1,1,2);
end
z = 0.5*h(k)*(zeta-1);
c1 = 0.9988; c2 = 0.8814;
U = (1-c1*cosh(c2*z).^(-2));
%% Run Solver
for i = 1:length(k)
    [o, an, cA, errGEP, db] = wZhang_solver(N,k(i),h(k(i)),Re,Fr2,method,alg,do_balancing); 
    fprintf('k = %.2f, growth rate = %.4f\n', k(i), real(o)/k(i));
    switch lower(method(1))
        case {'d2', 'd4'}
            mk(i).phi = D(:,:,1)*an(1:N+1);
            mk(i).up = D(:,:,2)*an(1:N+1);
            mk(i).wp = -1i*k(i)*mk(i).phi;
        case 'uw'
            mk(i).up = an(1:N+1);
            mk(i).wp = an(N+2:end);
            mk(i).phi = -mk(i).wp/k(i)/(1i);
    end
end
save('modeshape','mk','z','U','N','k','h');