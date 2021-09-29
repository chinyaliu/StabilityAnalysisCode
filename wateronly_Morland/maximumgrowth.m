clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Inputs
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,~,~,~,~,f] = pars_Morland(2);
%% change delta
delta_list = linspace(0,0.8,100);
lambda_list = NaN(1,length(delta_list));
c_list = NaN(1,length(delta_list));
t1 = tic;
for i = 1:length(delta_list)
    fprintf('delta = %.2f\n', delta_list(i));
    init_var = {N,6,ud_nd,delta_list(i),1,method,bflow};
    sol_vars = {alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',0.1,'eps',0.1)};
    [lambda,c] = findmaxgrowth(init_var, sol_vars, 10);
    if ~isnan(lambda)
        lambda_list(i) = lambda;
        c_list(i) = c;
    end
end
toc(t1);
o_list = 2*pi*c_list./lambda_list;

%% plot growth rate
figure('position',[50 50 720 640]);
plot(delta_list,imag(o_list),'-o');
xlabel('$\Delta / \lambda_m$');
ylabel('$(\lambda_m /c_m)\omega_i$');
xlim([0 0.8]);
ylim([0 1.5]);
grid on; axis square;

%% plot phase speed
figure('position',[50 50 720 640]);
plot(delta_list,real(c_list),'-o');
xlabel('$\Delta / \lambda_m$');
ylabel('$c_r / c_m$');
xlim([0 0.8]);
ylim([0 0.5]);
grid on; axis square;

%% plot lambda
figure('position',[50 50 720 640]);
plot(delta_list,lambda_list,'-o');
xlabel('$\Delta / \lambda_m$');
ylabel('$\lambda / \lambda_m$');
xlim([0 0.8]);
ylim([0 2.5]);
grid on; axis square;