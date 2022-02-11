clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,H,k,Fr2,Re,eps,c0,h,f] = pars_wake(3);
inflec_pt = -0.74708299-H;
zL = (0.74708299+H)*ones(length(k),1);
cutz = NaN(1,length(k));
cutz(1) = -inflec_pt;
in_init = {N,H,k(1),h(1),Re,Fr2};
addvar = struct('zL1',zL(1),'eps',eps);

%% Run solver
tic;
p1 = wSubmerged(in_init{:});
p1.numMeth(method);
o = NaN(1,length(k));
j = 1;
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i);
%     addvar.zL1 = zL(i);
    o(i) = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    if isnan(p1.zc)
        addvar.zL1 = cutz(j);
    else
        cutz(i)=p1.zc;
        j = i;
        addvar.zL1 = p1.zc;
    end
    fprintf('k = %.2f\n', k(i));
end
toc;

%% Plot oi vs or
fig1 = figure('position',[50,0,1000,720]);
plot(real(o),imag(o),'k.','Markersize',6);
hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot oi vs k
fig2 = figure('position',[50,0,1000,720]);
o(isnan(o))=0;
% plot(k,imag(o)./k,'-b','linewidth',3);
plot(k,imag(o),'-b','linewidth',3);
hold on; yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off'); hold off;
xlim([0 max(k)]);
xlabel('$k$','fontsize',30);
% ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot critical height vs k
[zc,ind] = rmmissing(cutz);
ind2 = k>0.87;
fig3 = figure('position',[50,0,1000,720]);
plot(k(~(ind&ind2)),-cutz(~(ind&ind2)),'k','linewidth',3);
% plot(k,-cutz,'k.','Markersize',6);
hold on; 
yline(inflec_pt,'linewidth',1.5,'color','b');
xlim([0 max(k)]);
xlabel('$k$','fontsize',30);
ylabel('$z_c$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot ci vs cr
c = o./k;
fig1 = figure('position',[50,0,1000,720]);
plot(real(c),imag(c),'k.','Markersize',6);
hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
xlabel('$c _r$','fontsize',30);
ylabel('$c _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');