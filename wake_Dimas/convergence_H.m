clear all;
if ~contains(path,'code_Dimas;')
    addpath('code_Dimas');
end    
%% Inout arguments
meth = ["qr", "qz", "eig"];
alg = meth(1);
bal = 'y';
eigspec = 'max';
k = 0.25;
Fr = 1.5;

%% Run solver
hlist = 0.5:0.25:4;
Nlist = 500:250:4000;
hfunc = hlist.*2*pi/real(k);
olist = nan(1,length(hlist));
tic;
for i = 1:length(hlist)
    h = hfunc(i);
    fprintf('h = %.2f times wave length\n', hlist(i));
    oi = [];
    for j = 1:length(Nlist)
        N = Nlist(j);
        [A, B] = makeAB(k, N, h, Fr);
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
        oi(j) = o1(1);
        if j~=1
            if (abs((oi(j)-oi(j-1))/oi(j))<1e-6 || j == length(Nlist))
                Nc(i) = Nlist(j);
                olist(i) = oi(j);
                break;
            end
        end
    end
end
toc;

%% Plot convergence
erro = abs((olist-olist(end))/olist(end));
figure;
semilogy(hlist, erro, '-o');
xlabel('$h/\lambda$');
ylabel('$\ | \ \omega_i(h_m) - \omega_i(h_{end})\ |/\omega_i(h_{end})$');
grid on;
title(sprintf('$k=%.2f$',k));
%% Plot convergence 2
erro2 = abs(diff(olist));
figure;
semilogy(hlist(1:end-1), erro2, '-o');
xlabel('$h/\lambda$');
ylabel('$\ | \ \omega_i(h_{m+1}) - \omega_i(h_m)\ |$');
title(sprintf('$k=%.2f$',k));
grid on;
%% Plot grid
figure;
semilogy(hlist, Nc, '-o');
xlabel('$h/\lambda$');
ylabel('grids required to converge');
grid on;