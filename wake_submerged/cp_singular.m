clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end
%% Solver & Algorithm list
[method,alg,bflow,~,do_balancing,~,~,H,k,Fr2,Re,eps,c0,h,f] = pars_wake;
N = 100:100:1500;
eig_spectrum = 'max';
in_init = {N(1),H,k,h,Re,Fr2,bflow};
%% Run
tic;
case1 = wSubmerged(in_init{:});
case1.numMeth(method);
zL = case1.criticalH(c0);
addvar = struct('zL1',zL,'eps',eps);
for i = 1:length(N)
    case1.N = N(i);
    o = case1.solver(alg,'n', do_balancing, eig_spectrum, f, addvar);
    oi(i) = imag(o(1));
    ob = case1.solver(alg,'y', do_balancing, eig_spectrum, f, addvar);
    oib(i) = imag(ob(1));
end
doi = abs(oi-oi(end));
doib = abs(oib-oib(end));
toc;
%% Plot figure
figure;
semilogy(N, doi, '-ko','linewidth',3, 'Displayname', 'original');
hold on;
semilogy(N, doib, '--ro','linewidth',3, 'Displayname', 'de-singularized');
hold off;
xlabel('$N$','fontsize',36);
ylabel('$\ | \ \omega_i - \omega_{0,i}\ |$','fontsize',36);
legend('location','northeast');
grid on;
xlim([100 1500]);
xticks(500:500:1500);
ylim([5e-16 2e-5]);