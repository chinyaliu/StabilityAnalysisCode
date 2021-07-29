clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end    
%% Set Solver & Algorithm
diff_meth = ["Schimd", "Trefethen"];
method = diff_meth(1);
solveGEPmeth = ["qr", "qz", "eig"];
alg = solveGEPmeth(1);
baseflowlist = ["exponential", "error function"];
bflow = baseflowlist(1);
% Inputs
do_balancing = 'y';
eig_spectrum = 'max';
N = 600;
ud_nd = 2;
delta_nd = 0.291;
h = @(x) 2*x;
ddm_number = 2;
f = wMorland.ddmtype(ddm_number);

%% Run solver
lambda_list = linspace(0.1,5,400);
c_list = NaN(1,length(lambda_list));
cutz_list = NaN(1,length(lambda_list));
t1 = tic;
flow1 = wMorland(N,h,ud_nd,delta_nd,1,method,bflow);
addvar.zL1 = delta_nd;
for i = 1:length(lambda_list)
    lam = lambda_list(i);
    fprintf('wavelength = %.2f\n', lam);
    flow1.k = 2*pi/lam;
    flow1.h = 3*delta_nd;
%     flow1.h = h(lam);
    c_list(i) = flow1.solver(alg, do_balancing, eig_spectrum, f, addvar);
    if isnan(flow1.zc)
        cutz_list(i) = NaN;
        addvar.zL1 = delta_nd;
    else
        cutz_list(i) = -flow1.zc;
        addvar.zL1 = cutz_list(i);
    end
end
toc(t1);
o_list = 2*pi*c_list./lambda_list;

%% Plot oi vs lambda
fig2 = figure('position',[50,0,1000,720]);
plot(lambda_list,imag(o_list),'k.','Markersize',6);
hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
xlabel('$\tilde{\lambda}$','fontsize',30);
ylabel('$\tilde{\omega _i}$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot critical height vs lambda
fig3 = figure('position',[50,0,1000,720]);
plot(lambda_list,-cutz_list,'k.','Markersize',6);
xlabel('$\tilde{\lambda}$','fontsize',30);
ylabel('$\tilde{z_c}$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
%% Plot ci vs cr
fig1 = figure('position',[50,0,1000,720]);
viscircles([ud_nd/2,0],ud_nd/2,'color','#898989');
hold on; 
yline(0,'linewidth',1.5,'color','#898989');
plot(real(c_list),imag(c_list),'k.','Markersize',6);
hold off;
xlabel('$\tilde{c _r}$','fontsize',30);
ylabel('$\tilde{c _i}$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');