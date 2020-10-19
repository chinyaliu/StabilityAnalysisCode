close all; clear all;% clc
tic;
load('neutral_pts.mat');
%% Solver & Algorithm list
order = ["D2", "D4", "uw"];
diff_method = ["Schimd", "Trefethen"];
constructAB_method = ["D4", "Schimd", "Herbert"];
solveGEPmethod = ["qr", "qz", "eig", "eigs", "polyeig", "singgep", "jdqz"];
%% Inputs
solver = [2,1,1]; % [order, diff_method, constructAB_method]
algorithm = 1;
do_balancing = 'n';
N = 401;
%% Run solver
method = [order(solver(1)), diff_method(solver(2)), constructAB_method(solver(3))];
alg = solveGEPmethod(algorithm);
%% Find k with largest growth rate at each Re
R_tar = logspace(4,9,46);
for i = 1:length(R_tar)
    fprintf('Re = %.3e\n',R_tar(i));
    [~,pos] = min(abs(Re0-R_tar(i)));
    if ( (pos+5)>length(Re0) || (pos-5)<1)
        k_tar = sort(k0(pos-10:pos));
    else
        k_tar = sort(k0(pos-5:pos+5));
    end
    kn = k_tar(1)-0.01:0.01:k_tar(end)+0.01;
    [o,~,~,~,~] = poiseuille_solver(N,kn(1),R_tar(i),method,alg,do_balancing);
    oi(1) = imag(o);
    for j = 2:length(kn)
        [o,~,~,~,~] = poiseuille_solver(N,kn(j),R_tar(i),method,alg,do_balancing);
        oi(j) = imag(o);
        if (oi(j)<oi(j-1))
            kn_tar(i) = kn(j-1);
            break;
        end
    end
end
%% Find minimum point zj(U == c_r)
if (strcmpi(method(1),'d4') && strcmpi(method(2),'schimd') && strcmpi(method(3),'d4'))
    M = N-2;
else
    M = N;
end
z = cospi((0:1:M)/M)'; 
U = 1-z.^2;
for i = 1:length(kn_tar)
    [o, ~, cA(i), errGEP(i), db(i)] = poiseuille_solver(N,kn_tar(i),R_tar(i),method,alg,do_balancing);
    c_r(i) = real(o)/kn_tar(i);
    omax(i) = o;
    [minpt(i),minpos(i)] = min(abs(U(1:fix(N+1)/2)-c_r(i)));
    fprintf('k = %.3f, c_r = %.8f, c_i = %.8f\n', kn_tar(i), c_r(i),imag(o)/kn_tar(i));
end
Rmax = R_tar;
kmax = kn_tar;
save('neutral_pts.mat','Re0','k0','Rmax','kmax','omax','U','z','N');
toc;