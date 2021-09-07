clear;
%% Solver & Algorithm list
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
do_balancing = 'y';
eig_spectrum = 'max';
N = 300:100:1500;
k = 2;
Re = 1e3;
Fr2 = 2.25;
c0 = 1./sqrt(k*Fr2);
% zL = -wZhang_ddm.criticalH(c0); 
zL = 0.74708299;
% DDM numbers
numberofDDM = 4;
eps = 0.2;
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
    o_temp = 0;
    for n = N
        fprintf('N = %3d\n', n);
        case1.N = n;
        o = case1.solver(alg, do_balancing, eig_spectrum, f, addvar);
        if abs(o-o_temp)<1e-8
            break;
        elseif ~isnan(o)
            o_temp = o;
        end
        if ~isnan(case1.zc)
            addvar.zL1 = -case1.zc;
        end
    end
    oih(j) = o;
    h_Nlist(j) = n;
end
toc;
%% Plot diff(c_i) vs h
dc = abs(diff(imag(oih)));
hold on;
% figure;
semilogy(nh(1:end-1),dc,'-go', 'Displayname', sprintf('$k = %.1f$',k));
xlabel('$n\lambda$');
ylabel('$\| \omega_i(m)-\omega_i(m+1) \|$');
grid on;

%% Plot N vs h
figure;
plot(nh,h_Nlist,'-bo');
xlabel('$n\lambda$');
ylabel('$N$');
grid on;
