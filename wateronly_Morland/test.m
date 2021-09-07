clear;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig", "invB"];
alg = solveGEPmeth(4);
baseflowlist = ["exponential", "error function"];
bflow = baseflowlist(1);
% Inputs
de_singularize = 'y';
do_balancing = 'y';
eig_spectrum = 'max';
N = 1000;
ud_nd = 2;
delta_nd = 0.291;
lambda_nd = [0.817, 0.4, 1.4];
ddm_number = 44;
f = wMorland.ddmtype(ddm_number);
sol_vars = {alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.1)};

%% Run solver
h_list = linspace(0.5,5,30);
c_list = cell(1,length(lambda_nd));
tic;
for j = 1:length(lambda_nd)
    wcase(j) = wMorland(N,0,ud_nd,delta_nd,lambda_nd(j),method,bflow);
    c_list{j} = NaN(1,length(h_list));
end
for i = 1:length(h_list)
    fprintf('h = %.2f times wave length.\n',h_list(i));
    for j = 1:length(lambda_nd)
        wcase(j).h = h_list(i)*lambda_nd(j);
        c = wcase(j).solvers(sol_vars{:});
        c_list{j}(i) = c;
    end
end
toc;

%% Plot diff(c_i) vs h
figspec = {'-ko', '--bo', '-.ro'};
figure;
hold on;
for j = 1:length(lambda_nd)
    dc = abs(diff(imag(c_list{j})));
    nam = ['$\lambda = ' sprintf('%.3f$',lambda_nd(j))];
    plot(h_list(1:end-1),dc,figspec{j},'Displayname',nam);
end
hold off; box on; grid on;
set(gca,'YScale','log');
xlabel('$n\lambda$');
ylabel('$\| c_i(m)-c_i(m+1) \|$');
legend;
