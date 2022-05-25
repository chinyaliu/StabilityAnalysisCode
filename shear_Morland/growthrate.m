clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,delta_nd,lamtemp,c0,htemp,f,epss,Re] = pars_Morland;
% if strcmpi(bflow,'exponential')
%     lmaxmin = findneutral(ud_nd,delta_nd);
%     lambda_nd = flip(linspace(lmaxmin{1}(2),lmaxmin{1}(1),100));
% else
%     lambda_nd = flip(linspace(0.2,2,100));
% end
lambda_nd = flip(linspace(0.6,1.3,100));
h = lambda_nd*htemp/lamtemp;
in_init = {N,h(1),ud_nd,delta_nd,lambda_nd(1),method,bflow,Re};

%% Run solver
tic;
p1 = wMorland(in_init{:});
zL = p1.invbf(c0);
cutz = NaN(1,length(lambda_nd));
cutz(1) = zL(1);
addvar = struct('zL1',zL(1),'eps',epss);
o = NaN(1,length(lambda_nd));
m = 1;
for i = 1:length(lambda_nd)
    p1.setprop('lambda',lambda_nd(i),'h',h(i));
%     p1.setprop('k',p1.k-0.01i);
    o(i) = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    if isnan(p1.zc)
        addvar.zL1 = cutz(m);
    else
        cutz(i) = p1.zc;
        m = i;
    end
    fprintf('lambda = %.3f\n', lambda_nd(i));
end
indgood = ~isnan(cutz);
minnan = find(indgood,1,'first');
maxnan = find(indgood,1,'last');
midnan = boolean([zeros(1,minnan-1) ones(1,maxnan-minnan+1) zeros(1,length(cutz)-maxnan)]);
cutznew = interp1(lambda_nd(indgood), cutz(indgood), lambda_nd(midnan&~indgood));
indqq = find(midnan&~indgood);
for i = 1:length(indqq)
    p1.setprop('lambda',lambda_nd(indqq(i)),'h',h(indqq(i)));
    addvar.zL1 = cutznew(i);
    o(indqq(i)) = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    if ~isnan(p1.zc)
        cutz(indqq(i)) = p1.zc;
    end
    fprintf('lambda = %.3f\n', lambda_nd(indqq(i)));
end
toc;
k = 2*pi./lambda_nd;
c = o./k;

%% Plot oi vs or
fig1 = figure('position',[50,0,1000,720]);
plot(real(o),imag(o),'k.','Markersize',6);
hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
 set(gca,'XMinorTick','on','YMinorTick','on');
 
%% Plot oi vs lambda
oi = imag(o);
fig2 = figure;
% plot(lambda_nd,oi,'k.','Markersize',6);
plot(lambda_nd,oi,'k');
hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
xlabel('$\lambda$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
 set(gca,'XMinorTick','on','YMinorTick','on');
 
%% Plot critical height vs lambda
fig3 = figure('position',[50,0,1000,720]);
plot(lambda_nd,cutz,'k.','Markersize',6);
xlabel('$\lambda$','fontsize',30);
ylabel('$z_c$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
set(gca,'XMinorTick','on','YMinorTick','on');
 
%% Plot ci vs cr
fig1 = figure('position',[50,0,1000,720]);
plot(real(c),imag(c),'k.','Markersize',6);
hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
set(gca,'XMinorTick','on','YMinorTick','on');
xlabel('$c_r$','fontsize',30);
ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');