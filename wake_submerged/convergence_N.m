clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,~,H,k,Fr2,Re,eps,c0,h,f] = pars_wake(2);
N = 100:100:1500;

%% Run solver
tic;
oi = nan(length(k),length(N));
cAall = nan(length(k),length(N));
ln = length(N);
parfor j = 1:length(k)
    nn = N;
    for i = 1:ln
        case1 = wSubmerged(nn(i),H,k(j),h(j),Re,Fr2,bflow);
        case1.numMeth(method);
        addvar = struct('zL1',case1.invbf(c0(j)),'eps',eps);
        [o,~,cA] = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        cAall(j,i) = cA;
        oi(j,i) = imag(o(1));
    end
end
doi = abs(oi-oi(:,end));
toc;
%% Plot figure
col = {[0 0 0],[0 0 1],[1 0 0],'#4ab712'};
lsy = {'-o','--o',':o','-.o'};
figure;
yf(1) = semilogy(N, cAall(1,:), lsy{1}, 'color',col{1}, 'Displayname', sprintf('$k = %.2f$',k(1)));
hold on;
for i = 2:length(k)
    nam = sprintf('$k = %.2f$',k(i));
    yf(i) = semilogy(N, cAall(i,:), lsy{i}, 'color',col{i}, 'Displayname', nam);
end
hold off;
set(yf(:),'linewidth',3,'markersize',8);
set(gca,'fontsize',30);
xticks(500:500:1500);
xlabel('$N$','fontsize',36);
ylabel('$cond(A)$','fontsize',36);
legend('location','southeast');
grid on;
%% Plot figure
figure;
yf(1) = semilogy(N, doi(1,:), lsy{1}, 'color',col{1}, 'Displayname', sprintf('$k = %.2f$',k(1)));
hold on;
for i = 2:length(k)
    nam = sprintf('$k = %.2f$',k(i));
    yf(i) = semilogy(N, doi(i,:), lsy{i}, 'color',col{i}, 'Displayname', nam);
end
hold off;
set(yf(:),'linewidth',3,'markersize',8);
set(gca,'fontsize',30);
ylim([1e-17 1e-2]);
xlim([N(1) N(end)]);
xticks(500:500:1500);
xlabel('$N$','fontsize',36);
ylabel('$\ | \ \omega_i - \omega_{0,i}\ |$','fontsize',36);
legend('location','northeast');
grid on;