clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
baseflowlist = ["exponential", "error function"];
bflow = baseflowlist(1);
% Inputs
do_balancing = 'y';
eig_spectrum = 'max';
ud_nd = 2;
delta_nd = 0.291;
lambda_nd = 0.817;
h = 2.5*lambda_nd;
ddm_number = 4;
addvar = struct('zL1',delta_nd,'eps',0.2);
f = wMorland.ddmtype(ddm_number);
init_var = {1000,h,ud_nd,delta_nd,lambda_nd,method,bflow};
sol_var = {alg, do_balancing, eig_spectrum, f, addvar};
%% Run solver
N_list = 100:100:2100;
c_list = NaN(1,length(N_list));
tic;
case1 = wMorland(init_var{:});
for i = 1:length(N_list)
    fprintf('N = %3d\n', N_list(i));
    case1.N = N_list(i);
    c_list(i) = case1.solver(sol_var{:});
end
ci = imag(c_list);
dci = abs(ci - ci(end).*ones(1,length(ci)))./ci(end);
toc;
%% Plot figure
figure;
semilogy(N_list, dci, '-o');
xlabel('$N$');
ylabel('$\ | \ c_i(N_m) - c_i(N_{end})\ |/c_i(N_{end})$');
grid on;

%% Plot c_i vs h
figure;
plot(N_list,ci,'-o');
xlabel('$N$');
ylabel('$c_i$');
grid on;