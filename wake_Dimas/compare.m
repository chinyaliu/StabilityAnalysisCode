%% Link needed functions
if ~contains(path,'code_Dimas;')
    addpath('code_Dimas');
end
pat = fullfile('../wakeflow_Zhang_orr');
if ~contains(path,'wakeflow_Zhang_orr')
    addpath(pat);
end    
%% Input arguments
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
do_balancing = 'y';
eig_spectrum = 'max';
Re = inf;
N_fd = 1500;
h = @(k) 2*2*pi./k;

%% Input arguments specifically for Chebyshev
N_cheby = 600;
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
inflec_pt = -0.74708299;
numberofDDM = 4;
eps = 0.01;
f = wZhang_ddm.ddmtype(numberofDDM);
addvar = struct('zL1',-inflec_pt,'eps',eps);

%% Run solver
klist = linspace(0.01,4,400);
o_cheby = nan(1,length(klist));
o_fd = nan(1,length(klist));
p1 = wZhang_ddm(N_cheby,1,1,Re,Fr^2);
p1.numMeth(method);
tic;
for i = 1:length(klist)
    k = klist(i);
    fprintf('k = %.2f\n', k);
    p1.k = k; p1.h = h(k);
    o_cheby(i) = p1.solver(alg, do_balancing, eig_spectrum, f, addvar);
    if isnan(p1.zc)
        addvar.zL1 = -inflec_pt;
    else
        addvar.zL1 = -p1.zc;
    end
    o_fd(i) = runFD(N_fd,k,h(k),Fr,alg,do_balancing,eig_spectrum);
end
toc;

%% Plot growth rate
figure;
plot(klist,imag(o_cheby),'-ko','markersize',4);
hold on;
plot(klist,imag(o_fd),'--r.','markersize',8);
hold off;
xlabel('$k$');
ylabel('$\omega_i$');
grid on;
legend("Chebyshev","Finite Difference",'location','northeast');

function o = runFD(N,k,h,Fr,alg,bal,eigspec)
    [A, B] = makeAB_fs(k, N, h, Fr);
    if strcmpi(bal,'y')
        o = balanceAB(A, B, eigspec, alg);
    else
        o = solveGEP(A, B, eigspec, alg);
    end
end