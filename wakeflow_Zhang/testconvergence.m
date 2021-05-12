close all; clear all;% clc
%% Solver & Algorithm list
order = ["Ray","D4"];
diff_meth = ["Schimd", "Trefethen"];
makeAB_meth = ["D4", "Schimd"];
solveGEPmeth = ["qr", "qz", "eig"];
%% Set solver
method = [order(1), diff_meth(1), makeAB_meth(1)];
alg = solveGEPmeth(1);
%% Inputs
do_balancing = 'n';
Re = inf;
Fr2 = 2.25;
N = 20:20:1200;
% k = [0.01 0.2 0.4 0.6];
% k = [0.87 1 2 3];
k = 3;
h = 2*pi./k;
eps = 0.15;
c0 = 1./sqrt(k*Fr2);
% zL = real(wZhang_ddm.g(c0)); 
zL = 0.74708299;
numberofDDM = 4;
addvar = struct('zL1',zL(1),'eps',eps); % ddm=4
% addvar = struct('zL1',zL(1)); % ddm = 2
f = wZhang_ddm.ddmtype(numberofDDM);
in_init = {N,k,h,Re,Fr2,method};
%% Run solver
tic;
case1 = wZhang_ddm(in_init{:});
for j = 1:length(k)
    fprintf('k = %.2f\n',k(j));
    case1.k = k(j); case1.h = h(j); 
%     addvar.zL1 = zL(j);
%     in_solver = {alg, do_balancing, f, addvar};
    for i = 1:length(N)
        case1.N = N(i);
        o = case1.solver(alg, do_balancing, f, addvar);
        oi(i) = imag(o(1));
        fprintf('N = %4d, growth rate = %.8f\n', N(i), oi(i));
    end
    oiall(j,:) = oi;
    doi{j} = abs(diff(oi));
end
toc;
% %% Plot figure
% fig1 = figure('position',[50,0,1000,720]);
% % doi = abs(diff(oi));
% nam = sprintf('k = %.2f',k(1));
% semilogy(N(1:end-1), doi{1}, '-o', 'linewidth', 1, 'Displayname', nam);
% hold on;
% for i = 2:length(k)
%     nam = sprintf('k = %.2f',k(i));
%     semilogy(N(1:end-1), doi{i}, '-o', 'linewidth', 1, 'Displayname', nam);
% end
% hold off;
% set(gca,'fontsize',20);
% xlabel('$N$', 'Interpreter', 'LaTeX','fontsize',30);
% ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{m+1})\ |$', 'Interpreter', 'LaTeX','fontsize',30);
% legend('location','northeast');
% grid on;
% % exportgraphics(fig1, 'fig_convergence\branch2.png');
%% Plot figure 2
dall = abs(oiall - oiall(:,end).*ones(size(oiall,1),size(oiall,2)))./oiall(:,end);
% fig1 = figure('position',[50,0,1080,840]);
figure;
nam = [sprintf('$k = %.2f, ',k(1)) '\' sprintf('omega = %0.4e$',oiall(1,end))];
% nam = sprintf('$k = %.2f$',k(1));
semilogy(N, dall(1,:), '-o', 'Displayname', nam);
hold on;
for i = 2:length(k)
%     nam = sprintf('$k = %.2f$',k(i));
    nam = [sprintf('$k = %.2f, ',k(i)) '\' sprintf('omega = %0.4e$',oiall(i,end))];
    semilogy(N, dall(i,:), '-o', 'Displayname', nam);
end
hold off;
xlabel('$N$');
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{end})\ |/\omega_i(N_{end})$');
ylim([1e-15 1])
legend('location','northeast');
grid on;