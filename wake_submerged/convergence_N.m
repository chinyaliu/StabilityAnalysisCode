clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,~,H,k,Fr2,Re,eps,~,h,f] = pars_wake(2);
N = 300:100:1500;

%% Run solver
tic;
oi = nan(length(k),length(N));
cAall = nan(length(k),length(N));
ln = length(N);
parfor j = 1:length(k)
    fprintf('k = %.2f\n',k(j));
    addvar = struct('zL1',0.74708299+H,'eps',eps);
    case1 = wSubmerged(100,H,k(j),h(j),Re,Fr2);
    case1.numMeth(method);
    nn = 300:100:1500;
    for i = 1:ln
        case1.N = nn(i);
        [o,~,cA] = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        cAall(j,i) = cA;
        oi(j,i) = imag(o(1));
%         fprintf('N = %3d, growth rate = %.8f\n', N(i), oi(i));
    end
%     cAall(j,:) = cA;
%     doi{j} = abs(diff(oi));
%     doi{j} = abs(oi-oi(end));
end
doi = abs(oi-oi(:,end));
toc;
%% Plot figure
col = {[0 0 0],[0 0 1],[1 0 0],"#77AC30"};
figure;
semilogy(N, cAall(1,:), '-ko', 'Displayname', sprintf('$k = %.2f$',k(1)));
hold on;
for i = 2:length(k)
    nam = sprintf('$k = %.2f$',k(i));
    semilogy(N, cAall(i,:), '-o', 'Displayname', nam, 'color',col{i});
end
hold off;
xlabel('$N$');
ylabel('$cond(A)$');
legend('location','southeast');
grid on;
%% Plot figure
figure;
semilogy(N, doi(1,:), '-ko', 'Displayname', sprintf('$k = %.2f$',k(1)));
hold on;
for i = 2:length(k)
    nam = sprintf('$k = %.2f$',k(i));
    semilogy(N, doi(i,:), '-o', 'Displayname', nam, 'color',col{i});
end
hold off;
xlabel('$N$');
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{end})\ |$');
legend('location','northeast');
grid on;