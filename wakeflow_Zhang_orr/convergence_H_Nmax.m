clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,~,k,Fr2,Re,eps,c0,~,f] = pars_wake(2);
N = 300:100:1500;
nh = linspace(0.5,5,20);
h = 2*pi/k*nh;
in_init = {N(1),k(1),h(1),Re,Fr2};
%% Run solver
tic;
case1 = wZhang_ddm(in_init{:});
case1.numMeth(method);
addvar = struct('zL1',-case1.criticalH(c0),'eps',eps);
for j = 1:length(nh)
    fprintf('h = %.2f times wave length.\n',nh(j));
    case1.h = h(j); 
    o_temp = 0;
    for n = N
        fprintf('N = %3d\n', n);
        case1.N = n;
        o = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
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
figure;
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
