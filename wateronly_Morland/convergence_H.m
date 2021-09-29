clear;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,~,ud_nd,delta_nd,lambda_nd,c0,~,f] = pars_Morland(2);
init_var = {1000,0,ud_nd,delta_nd,lambda_nd,method,bflow};
sol_vars = {alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',delta_nd,'eps',0.1)};

%% Run solver
h_list = linspace(0.5,5,30);
N_list = 300:100:1500;
c_list = NaN(1,length(h_list));
h_Nlist = NaN(1,length(h_list));
tic;
case1 = wMorland(init_var{:});
addvar.zL1 = -case1.criticalH(c0);
for i = 1:length(h_list)
    fprintf('h = %.2f times wave length.\n',h_list(i));
    case1.h = h_list(i)*lambda_nd;
    c_temp = 0;
    for N = N_list
        fprintf('N = %3d\n', N);
        case1.N = N;
        c = case1.solvers(sol_vars{:});
        if abs(c-c_temp)<1e-8 %&& ~isnan(c_temp)
            break;
        elseif ~isnan(c)
            c_temp = c;
        end
        if ~isnan(case1.zc)
            sol_vars{6}.zL1 = -case1.zc;
        end
    end
    h_Nlist(i) = N;
    c_list(i) = c;
end
toc;

%% Plot N vs h
figure;
plot(h_list,h_Nlist,'-bo');
xlabel('$n\lambda$');
ylabel('$N$');
grid on;

%% Plot c_i vs h
figure;
plot(h_list,imag(c_list),'-bo');
xlabel('$n\lambda$');
ylabel('$c_i$');
grid on;

%% Plot diff(c_i) vs h
dc = abs(diff(imag(c_list)));
figure;
semilogy(h_list(1:end-1),dc,'-bo');
xlabel('$n\lambda$');
ylabel('$\| c_i(m)-c_i(m+1) \|$');
grid on;

%% Plot figure
ci = imag(c_list);
dci = abs(ci - ci(end).*ones(1,length(ci)))./ci(end);
figure;
semilogy(h_list, dci, '-bo');
xlabel('$n\lambda$');
ylabel('$\ | \ c_i(N_m) - c_i(N_{end})\ |/c_i(N_{end})$');
grid on;