clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,~,~,~,H,k,Fr2,Re,eps,c0,h,f] = pars_wake;
N = 300:100:2000;
eig_spectrum = 'max';
in_init = {N(1),H,k,h,Re,Fr2};
%% Run
tic;
case1 = wSubmerged(in_init{:});
case1.numMeth(method);
zL = case1.criticalH(c0);
addvar = struct('zL1',zL,'eps',eps);
for i = 1:length(N)
    case1.N = N(i);
    o = case1.solver(alg, de_singularize, 'n', eig_spectrum, f, addvar);
    oi(i) = imag(o(1));
    ob = case1.solver(alg, de_singularize, 'y', eig_spectrum, f, addvar);
    oib(i) = imag(ob(1));
end
doi = abs(oi-oi(end));
doib = abs(oib-oib(end));
toc;
%% Plot figure
figure;
semilogy(N, doi, '-ko','linewidth',3, 'Displayname', 'original');
hold on;
semilogy(N, doib, '--ro','linewidth',3, 'Displayname', 'balanced');
hold off;
xlabel('$N$');
ylabel('$\ | \ \omega_i(N_m) - \omega_i(N_{end})\ |$');
legend('location','northeast');
grid on;
xlim([300 2000]);