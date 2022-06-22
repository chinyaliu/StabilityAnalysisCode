clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,~,~,ud_nd,delta_nd,lambda_nd,c0,h,f,epss,Re] = pars_Morland(6,'vis');
N = [100 160 200:100:1500];
% N = 100:100:1500;
eig_spectrum = 'max';
lnk = length(lambda_nd);
ln = length(N);
if length(ud_nd)==1
    ud_nd = ud_nd*ones(1,lnk);
    delta_nd = delta_nd*ones(1,lnk);
    Re = Re*ones(1,lnk);
end

%% Run solver
tic;
oi = nan(lnk,ln);
oi2 = nan(lnk,ln);
parfor j = 1:lnk
    nn = N;
    oii = nan(1,ln);
    for i = 1:ln
        case1 = wMorland(nn(i),h(j),ud_nd(j),delta_nd(j),lambda_nd(j),method,bflow,Re(j));
        addvar = struct('zL1',case1.invbf(c0(j)),'eps',epss);
        o = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        oii(i) = imag(o(1));
    end
    oi2(j,:) = oii;
    oi(j,:) = Mortow(ud_nd(j),delta_nd(j),lambda_nd(j),Re(j),oii);
end
doi = abs(oi-oi(:,end));
doi2 = abs(oi2-oi2(:,end));
toc;

%% Plot figure
lsy = {'-ko',':ro','--bo','-.go'};
figure;
% nam = sprintf("$U_0 = %.2f, ",ud_nd(1))+"\"+sprintf("Delta = %.3f,",delta_nd(1))+"\"+sprintf("lambda = %.3f$",lambda_nd(1));
u = [30 30 80 80];
lam = [4 20 4 20];
nam = sprintf("$u_{*}=%3d\\ \\textrm{cm/s}, \\lambda=%3.0f\\ \\textrm{cm}$",u(1),lam(1));
yf(1) = semilogy(N, doi(1,:), lsy{1}, 'Displayname',nam);
hold on;
for i = 2:lnk
%     nam = sprintf("$U_0 = %.2f, ",ud_nd(i))+"\"+sprintf("Delta = %.3f,",delta_nd(i))+"\"+sprintf("lambda = %.3f$",lambda_nd(i));
    nam = sprintf("$u_{*}=%3d\\ \\textrm{cm/s}, \\lambda=%3.0f\\ \\textrm{cm}$",u(i),lam(i));
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