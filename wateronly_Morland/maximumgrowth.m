clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Inputs
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
baseflowlist = ["exponential", "error function"];
bflow = baseflowlist(1);
do_balancing = 'y';
eig_spectrum = 'max';
N = 600;
ud_nd = 2;
ddm_number = 4;
f = wMorland.ddmtype(ddm_number);

%% change delta
delta_list = linspace(0.11,0.8,70);
lambda_list = NaN(1,length(delta_list));
c_list = NaN(1,length(delta_list));
t1 = tic;
inilam = 0.4;
for i = 1:length(delta_list)
    fprintf('delta = %.2f\n', delta_list(i));
    init_var = {N,6,ud_nd,delta_list(i),1,method,bflow};
    sol_var = {alg, do_balancing, eig_spectrum, f, struct('zL1',delta_list(i),'eps',0.2)};
    [lambda,c] = findmaxgrowth(init_var, sol_var, 10, inilam, 0.2);
    if ~isnan(lambda)
        lambda_list(i) = lambda;
        c_list(i) = c;
        inilam = lambda;
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