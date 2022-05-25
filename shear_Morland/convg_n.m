clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,~,ud_nd,delta_nd,lambda_nd,c0,h,f,epss,Re] = pars_Morland(3);
N = 100:100:1500;
lnk = length(lambda_nd);
ln = length(N);
if length(ud_nd)==1
    ud_nd = ud_nd*ones(1,lnk);
    delta_nd = delta_nd*ones(1,lnk);
end

%% Run solver
tic;
oi = nan(lnk,ln);
parfor j = 1:lnk
    nn = N;
    for i = 1:ln
        case1 = wMorland(nn(i),h(j),ud_nd(j),delta_nd(j),lambda_nd(j),method,bflow,Re);
        addvar = struct('zL1',case1.invbf(c0(j)),'eps',epss);
        o = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oi(j,i) = imag(o(1));
    end
end
doi = abs(oi-oi(:,end));
toc;

%% Plot figure
lsy = {'-ko',':ro','--bo','-.go'};
figure;
yf(1) = semilogy(N, doi(1,:), lsy{1}, 'Displayname', "$\"+sprintf("lambda = %.3f$",lambda_nd(1)));
hold on;
for i = 2:lnk
    nam = "$\"+sprintf('lambda = %.3f$',lambda_nd(i));
    yf(i) = semilogy(N, doi(i,:), lsy{i}, 'Displayname', nam);
end
hold off;
set(yf(:),'linewidth',3,'markersize',8);
set(gca,'fontsize',30);
ylim([1e-15 1e-1]);
xlim([N(1) N(end)]);
xticks(500:500:1500);
xlabel('$N$','fontsize',36);
ylabel('$\ | \ \omega_i - \omega_{0,i}\ |$','fontsize',36);
legend('location','northeast');
grid on;