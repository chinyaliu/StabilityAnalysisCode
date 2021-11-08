clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,~,k,Fr2,Re,eps,~,~,f] = pars_wake(2);
N = 300:100:1500;
% k = linspace(0.1,4,10);
h = 3*2*pi./real(k);
c0 = min(0.99,1./sqrt(k*Fr2));
in_init = {N(1),k(1),h(1),Re,Fr2};
%% Run solver
tic;
case1 = wZhang_ddm(in_init{:});
case1.numMeth(method);
zL = -case1.criticalH(c0);
addvar = struct('zL1',zL(1),'eps',eps);
for j = 1:length(k)
    fprintf('k = %.2f\n',k(j));
    case1.k = k(j); case1.h = h(j); 
    addvar.zL1 = zL(j); 
    for i = 1:length(N)
        case1.N = N(i);
        [o,~,cA(i)] = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oi(i) = imag(o(1));
        fprintf('N = %3d, growth rate = %.8f\n', N(i), oi(i));
    end
    oiall(j,:) = oi;
    cAall(j,:) = cA;
    doi{j} = abs(diff(oi));
end
toc;
% %% Plot figure
% dall = abs(oiall - oiall(:,end));
% figure;
% % hold on;
% % nam = [sprintf('$k = %.2f, ',k(1)) '\' sprintf('omega = %0.9e$',oiall(1,end))];
% nam = ['$\' sprintf('omega = %0.9e,',oiall(1,end)) ' \ ' sprintf('Re = %1.1e$',Re)];
% semilogy(N, dall(1,:), '-go', 'Displayname', nam);
% hold on;
% for i = 2:length(k)
%     nam = [sprintf('$k = %.2f, ',k(i)) '\' sprintf('omega = %0.9e$',oiall(i,end))];
%     semilogy(N, dall(i,:), '-o', 'Displayname', nam);
% end
% hold off;
% xlabel('$N$');
% ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{end})\ |$');
% legend('location','northeast');
% grid on;
%% Plot figure
figure;
% hold on;
% nam = sprintf('$k = %.2f$',k(1));
nam = sprintf('$Re = %1.1e$',Re);
semilogy(N, cAall(1,:), '-ro', 'Displayname', nam);
hold on;
for i = 2:length(k)
    nam = sprintf('$k = %.2f, ',k(i));
    semilogy(N, cAall(1,:), '-o', 'Displayname', nam);
end
hold off;
xlabel('$N$');
ylabel('$cond(A)$');
legend('location','northeast');
grid on;
%% Plot figure
figure;
% hold on;
% nam = sprintf('$k = %.2f, ',k(1));
nam = ['$\' sprintf('omega = %0.9e,',oiall(1,end)) ' \ ' sprintf('Re = %1.1e$',Re)];
semilogy(N(1:end-1), doi{1}, '-ro', 'Displayname', nam);
hold on;
for i = 2:length(k)
    nam = sprintf('$k = %.2f$',k(i));
    semilogy(N(1:end-1), doi{i}, '-o', 'Displayname', nam);
end
hold off;
xlabel('$N$');
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{m+1})\ |$');
legend('location','northeast');
grid on;