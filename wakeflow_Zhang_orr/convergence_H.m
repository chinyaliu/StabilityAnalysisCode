close all; clear all;% clc
%% Solver & Algorithm list
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
do_balancing = 'y';
eig_spectrum = 'max';
N = 300:100:1500;
k = 3;
Re = inf;
Fr2 = 2.25;
c0 = 1./sqrt(k*Fr2);
zL = -wZhang_ddm.criticalH(c0); 
% zL = 0.74708299;
% DDM numbers
numberofDDM = 4;
eps = 0.15;
f = wZhang_ddm.ddmtype(numberofDDM);
% truncation height
nh = linspace(0.5,5,20);
% nh = linspace(0.1,8,30);
h = 2*pi/k*nh;
in_init = {N(1),k(1),h(1),Re,Fr2};
addvar = struct('zL1',zL,'eps',eps);
%% Run solver
tic;
case1 = wZhang_ddm(in_init{:});
case1.numMeth(method);
for j = 1:length(nh)
    fprintf('h = %.2f times wave length.\n',nh(j));
    case1.h = h(j); 
    for i = 1:length(N)
        case1.N = N(i);
        o = case1.solver(alg, do_balancing, eig_spectrum, f, addvar);
        oi(i) = imag(o);
        fprintf('N = %3d, growth rate = %.8f\n', N(i), oi(i));
        if i~=1
            if (abs((oi(i)-oi(i-1))/oi(i))<1e-10 || i == length(N))
                Nc(j) = N(i);
                break;
            end
        end
        if ~isnan(case1.zc)
            addvar.zL1 = -case1.zc;
        end
    end
    oih(j) = oi(i);
end
toc;
%% Plot diff(c_i) vs h
dc = abs(diff(imag(oih)));
figure;
semilogy(nh(1:end-1),dc,'-bo');
xlabel('$n\lambda$');
ylabel('$\| \omega_i(m)-\omega_i(m+1) \|$');
grid on;
%% Plot figure
% dall = abs(oih-oih(end))./oih(end);
dall = abs(oih-oih(end));
figure;
nam = [sprintf('$k = %.2f, ',k(1)) '\' sprintf('omega = %0.9e$',oih(end))];
semilogy(nh, dall, '-o');
xlabel('$h/\lambda$');
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{end})\ |/\omega_i(N_{end})$');
nam = [sprintf('$k = %.2f, ',k) '\' sprintf('omega = %0.9e$',oih(end))];
legend(nam,'location','northeast');
xlim([1 8]);
grid on;
%% Plot figure
figure;
nam = [sprintf('$k = %.2f, ',k(1)) '\' sprintf('omega = %0.9e$',oih(end))];
semilogy(nh, oih, '-o');
xlabel('$h/\lambda$');
ylabel('$\omega_i$');
nam = [sprintf('$k = %.2f, ',k) '\' sprintf('omega = %0.9e$',oih(end))];
legend(nam,'location','northeast');
grid on;
%% Plot figure
figure;
nam = [sprintf('$k = %.2f, ',k(1)) '\' sprintf('omega = %0.9e$',oih(end))];
semilogy(nh, Nc, '-o');
xlabel('$h/\lambda$');
ylabel('grids required to converge');
nam = [sprintf('$k = %.2f, ',k) '\' sprintf('omega = %0.9e$',oih(end))];
legend(nam,'location','southeast');
grid on;