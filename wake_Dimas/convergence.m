clear all;
%% Inout arguments
meth = ["qr", "qz", "eig"];
alg = meth(1);
bal = 'n';
eigspec = 'max';
k = 3;
Fr = 1.5;
h = @(k) 2*pi/real(k);

%% Run solver
nlist = 500:500:5000;
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
erro = abs((olist-olist(end))/olist(end));
toc;

%% Plot convergence
figure;
semilogy(nlist, erro, '-o');
xlabel('$N$');
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{end})\ |/\omega_i(N_{end})$');
grid on;