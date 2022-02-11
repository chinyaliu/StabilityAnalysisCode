clear all;
if ~contains(path,'wake_Zhang_orr\code_wake;')
    addpath('code_wake');
end 
%% Solver & Algorithm list
[method,alg,bflow,de_singularize,do_balancing,~,~,k,Fr2,Re,eps,c0,h,~] = pars_wake;
N = 300:100:2000;
eig_spectrum = 'max';
in_init = {N(1),k,h,Re,Fr2};
%% Run
tic;
case1 = wZhang_ddm(in_init{:});
case1.numMeth(method);
% zL = case1.criticalH(c0);
zL = 0.74708299;
addvar = struct('zL1',zL,'eps',eps);
for i = 1:length(N)
    case1.N = N(i);
    o = case1.solver(alg,de_singularize, do_balancing, eig_spectrum, wZhang_ddm.ddmtype(1), addvar);
    oi(i) = imag(o(1));
    o2 = case1.solver(alg,de_singularize, do_balancing, eig_spectrum, wZhang_ddm.ddmtype(2), addvar);
    oi2(i) = imag(o2(1));
    o4 = case1.solver(alg,de_singularize, do_balancing, eig_spectrum, wZhang_ddm.ddmtype(44), addvar);
    oi4(i) = imag(o4(1));
end
doi = abs(oi-oi(end));
doi2 = abs(oi2-oi2(end));
doi4 = abs(oi4-oi4(end));
toc;
%% Plot figure
figure;
semilogy(N, doi, '-ko','linewidth',3, 'Displayname', 'original');
hold on;
semilogy(N, doi2, '--bo','linewidth',3, 'Displayname', '2 domains');
semilogy(N, doi4, '-.ro','linewidth',3, 'Displayname', '4 domains');
hold off;
xlabel('$N$','fontsize',36);
ylabel('$\ | \ \omega_i - \omega_{0,i}\ |$','fontsize',36);
legend('location','northeast');
grid on;
xlim([300 2000]);
ylim([1e-17 1]);