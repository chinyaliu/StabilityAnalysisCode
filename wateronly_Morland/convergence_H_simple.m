clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig", "invB"];
alg = solveGEPmeth(4);
baseflowlist = ["exponential", "error function"];
bflow = baseflowlist(1);
% Inputs
de_singularize = 'y';
do_balancing = 'y';
eig_spectrum = 'max';
N = 1000;
ud_nd = 2;
delta_nd = 0.291;
lambda_nd = 0.817;
ddm_number = 44;
f = wMorland.ddmtype(ddm_number);
init_var = {N,0,ud_nd,delta_nd,lambda_nd,method,bflow};
sol_var = {alg, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.2)};
sol_vars = {alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.2)};

%% Run solver
h_list = linspace(0.5,5,30);
c_list = NaN(1,length(h_list));
tic;
case1 = wMorland(init_var{:});
for i = 1:length(h_list)
    fprintf('h = %.2f times wave length.\n',h_list(i));
    case1.h = h_list(i)*lambda_nd;
%     c = case1.solver(sol_var{:});
    c = case1.solvers(sol_vars{:});
    c_list(i) = c;
end
toc;

%% Plot c_i vs h
figure;
plot(h_list,imag(c_list),'-bo');
xlabel('$n\lambda$');
ylabel('$c_i$');
grid on;

%% Plot diff(c_i) vs h
dc = abs(diff(imag(c_list)));
figure;
semilogy(h_list(1:end-1),dc,'-bo');
xlabel('$n\lambda$');
ylabel('$\| c_i(m)-c_i(m+1) \|$');
grid on;

%% Plot figure
ci = imag(c_list);
dci = abs(ci - ci(end).*ones(1,length(ci)))./ci(end);
figure;
semilogy(h_list, dci, '-bo');
xlabel('$n\lambda$');
ylabel('$\ | \ c_i(N_m) - c_i(N_{end})\ |/c_i(N_{end})$');
grid on;