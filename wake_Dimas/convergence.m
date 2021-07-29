clear all;
if ~contains(path,'code_Dimas;')
    addpath('code_Dimas');
end    
%% Inout arguments
meth = ["qr", "qz", "eig"];
alg = meth(1);
bal = 'y';
eigspec = 'max';
k = 3.2;
Fr = 1.5;
h = @(k) 2*pi/real(k);

%% Run solver
nlist = 250:250:5250;
olist = nan(1,length(nlist));
tic;
for i = 1:length(nlist)
    N = nlist(i);
    fprintf('N = %.2f\n', N);
    [A, B] = makeAB(k, N, h(k), Fr);
    if strcmpi(bal,'y')
        [o, an] = balanceAB(A, B, eigspec, alg);
    else
        [o, an] = solveGEP(A, B, eigspec, alg);
    end
    if  strcmpi(eigspec,'all')
        c = o./k;
        o1 = k.*filt(c);
    else
        o1 = o;
    end
    olist(i) = o1(1);
end
toc;

%% Plot convergence
erro = abs((olist-olist(end))/olist(end));
figure;
semilogy(nlist, erro, '-o');
xlabel('$N$');
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{end})\ |/\omega_i(N_{end})$');
grid on;
title(sprintf('$k=%.2f$',k));
%% Plot convergence 2
erro2 = abs(diff(olist));
figure;
semilogy(nlist(1:end-1), erro2, '-o');
xlabel('$N$');
ylabel('$\ | \ \omega_i(N_{m+1}) - \omega_i(N_m)\ |$');
title(sprintf('$k=%.2f$',k));
grid on;