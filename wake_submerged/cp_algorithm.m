clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Solver & Algorithm list
[method,~,bflow,de_singularize,do_balancing,eig_spectrum,~,H,k,Fr2,Re,eps,c0,h,f] = pars_wake;
N = 100:100:1500;
in_init = {N(1),H,k,h,Re,Fr2,bflow};
% alg = ["eig", "qr"];
alg = ["eig", "qr", "invB"];
%% Run
tic;
case1 = wSubmerged(in_init{:});
case1.numMeth(method);
zL = case1.invbf(c0);
addvar = struct('zL1',zL,'eps',eps);
for i = 1:length(N)
    case1.N = N(i);
    for j = 1:length(alg)
        t1 = tic;
        o = case1.solver(alg(j), de_singularize, do_balancing, eig_spectrum, f, addvar);
        tt(i,j) = toc(t1);
        oi(i,j) = imag(o(1));
    end
end
doi = abs(oi-oi(end,:));
toc;
%% Plot convergence
lsy = {'-ko',':ro','--bo'};
disname = {'qz', 'inv(A)', 'inv(B)'};
figure;
yf(1) = semilogy(N, doi(:,1), lsy{1}, 'Displayname', disname{1});
hold on;
for i = 2:length(alg)
    yf(i) = semilogy(N, doi(:,i), lsy{i}, 'Displayname', disname{i});
end
set(yf(:),'linewidth',3,'markersize',8);
set(gca,'fontsize',30);
hold off;
xlabel('$N$','fontsize',36);
ylabel('$\ | \ \omega_i - \omega_{0,i}\ |$','fontsize',36);
ylim([1e-15 1]);
xlim([300 1500]);
xticks(500:500:1500);
legend('location','northeast');
grid on;
%% Plot time
figure;
yf(1) = plot(N, tt(:,1), lsy{1}, 'Displayname', disname{1});
hold on;
for i = 2:length(alg)
    yf(i) = plot(N, tt(:,i), lsy{i}, 'Displayname', disname{i});
end
set(yf(:),'linewidth',3,'markersize',8);
set(gca,'fontsize',30);
hold off;
xlabel('$N$','fontsize',36);
ylabel('computational time (s)','fontsize',36);
legend('location','northwest');
xlim([300 1500]);
xticks(500:500:1500);
grid on;