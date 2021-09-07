clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Test solver
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
baseflowlist = ["exponential", "error function"];
bflow = baseflowlist(1);
% Inputs
do_balancing = 'y';
eig_spectrum = 'max';
N = 600;
ud_nd = 2;
delta_nd = 0.6;
ddm_number = 4;
f = wMorland.ddmtype(ddm_number);
init_var = {N,6,ud_nd,delta_nd,1,method,bflow};
sol_var = {alg, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',2)};

t1 = tic;
[lambda,c] = findmaxgrowth(init_var, sol_var, 10);
toc(t1);