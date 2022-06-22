clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,~,N,ud_nd,delta_nd,lamtemp,c0,htemp,f,epss,Re] = pars_Morland(1,'inv');
eig_spectrum = 'all';
[lmin, lmax] = findneutral(ud_nd,delta_nd);
% if strcmpi(bflow,'exponential')
%     lambda_nd = flip(linspace(lmin,lmax,100));
%     lambda_nd = linspace(lmin,lmax,400);
% else
%     lambda_nd = flip(linspace(0.2,2,100));
% end
lambda_nd = linspace(lmin,lmax,100);
lambda_nd = flip([linspace(0.04,lmin,30) lambda_nd logspace(log(lmax),2,30)]);
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
o2 = NaN(1,length(lambda_nd));
m = 1;
for i = 1:length(lambda_nd)
    p1.setprop('lambda',lambda_nd(i),'h',h(i));
%     p1.setprop('k',p1.k-0.01i);
    oall = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    [~,ind] = min(real(oall));
    if real(oall(ind)) < 0
        o(i) = oall(ind);
    elseif imag(oall(1))>0
        o(i) = oall(1);
    end
    [~,ind] = max(real(oall));
    o2(i) = oall(ind);
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
    oall = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    [~,ind] = min(real(oall));
    if real(oall(ind)) < 0
        o(indqq(i)) = oall(ind);
    elseif imag(oall(1))>0
        o(indqq(i)) = oall(1);
    end
%     if imag(oall(1))>0
%         o(indqq(i)) = oall(1);
%     else
%         [~,ind] = min(real(oall));
%         if real(oall(ind)) <= 0
%             o(indqq(i)) = oall(ind);
%         end
%     end
    [~,ind] = max(real(oall));
    o2(indqq(i)) = oall(ind);
    if ~isnan(p1.zc)
        cutz(indqq(i)) = p1.zc;
    end
    fprintf('lambda = %.3f\n', lambda_nd(indqq(i)));
end
toc;
k = 2*pi./lambda_nd;
c = o./k;
c2 = o2./k;

% %% Plot oi vs or
% fig1 = figure('position',[50,0,1000,720]);
% plot(real(o),imag(o),'k.','Markersize',6);
% hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
% xlabel('$\omega _r$','fontsize',30);
% ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
%  set(gca,'XMinorTick','on','YMinorTick','on');
%  
% %% Plot oi vs lambda
% oi = imag(o);
% fig2 = figure;
% % plot(lambda_nd,oi,'k.','Markersize',6);
% plot(lambda_nd,oi,'k');
% hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
% xlabel('$\lambda$','fontsize',30);
% ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
%  set(gca,'XMinorTick','on','YMinorTick','on');
%  
% %% Plot critical height vs lambda
% fig3 = figure('position',[50,0,1000,720]);
% plot(lambda_nd,cutz,'k.','Markersize',6);
% xlabel('$\lambda$','fontsize',30);
% ylabel('$z_c$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% set(gca,'XMinorTick','on','YMinorTick','on');
%  
% %% Plot ci vs cr
% fig4 = figure('position',[50,0,1000,720]);
% plot(real(c),imag(c),'k.','Markersize',6);
% hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
% set(gca,'XMinorTick','on','YMinorTick','on');
% xlabel('$c_r$','fontsize',30);
% ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% %% Plot ci vs cr for both modes
% fig5 = figure('position',[50,0,1000,720]);
% plot(real(c),imag(c),'k.','Markersize',6);
% hold on;
% plot(real(c2),imag(c2),'k.','Markersize',6); 
% yline(0,'linewidth',1.5,'color','#898989'); 
% hold off;
% set(gca,'XMinorTick','on','YMinorTick','on');
% xlabel('$c_r$','fontsize',30);
% ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% %% Plot cr-k (mode 1)
% fig6 = figure('position',[50,0,1000,720]);
% plot(lambda_nd,real(c),'k.','Markersize',6);
% hold on; yline(0,'--','linewidth',1.5,'color','#898989'); hold off;
% set(gca,'XMinorTick','on','YMinorTick','on');
% xlabel('$\lambda$','fontsize',30);
% ylabel('$c_r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% %% Plot cr-k (mode 2)
% fig7 = figure('position',[50,0,1000,720]);
% plot(lambda_nd,real(c2),'k.','Markersize',6);
% cmin = min(c2);
% hold on; yline(cmin,'--','linewidth',1.5,'color','#898989'); hold off;
% set(gca,'XMinorTick','on','YMinorTick','on');
% xlabel('$\lambda$','fontsize',30);
% ylabel('$c_r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
%% Plot c vs lambda
ln = {'-','--',':'};
col = [0.16,0.44,1.00;0.93,0.00,0.00;0.00,0.57,0.00];
cc = [c;c2];
for j = 1:2
    fig = figure('position',[50,50,1620,540]);
    hh(j,1) = subplot(1,2,1);
    yline(0,'linewidth',1.5,'color','#898989','Handlevisibility','off');
    hold on;
    plot(lambda_nd,real(cc(j,:)),'linewidth',3,'color',col(j,:));
    hold off; box on;
    xlabel('$\lambda$','fontsize',30);
    ylabel('$c _r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    xticks(0:2:10);
    ax = gca;
    set(ax,'XMinorTick','on','YMinorTick','on');
    
    hh(j,2) = subplot(1,2,2);
    yline(0,'linewidth',1.5,'color','#898989','Handlevisibility','off');
    hold on;
    plot(lambda_nd,imag(cc(j,:)),'linewidth',3,'color',col(j,:));
    hold off; box on;
    xlabel('$\lambda$','fontsize',30);
    ylabel('$c _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    ax = gca;
    set(ax,'XMinorTick','on','YMinorTick','on');
    xticks(0:2:10);
end

%% Plot c_r
[cmin,ind] = min(c2);
cdrift = ud_nd + sqrt(0.5./lambda_nd);
csw = sqrt(0.5*lambda_nd);

hcr(1) = figure('position',[50,0,810,540]);
plot(lambda_nd,real(c2),'linewidth',3,'color',col(2,:),'Handlevisibility','off');
hold on;
plot(lambda_nd,cdrift,':k','linewidth',3,'Displayname','$U_0+\frac{1}{\sqrt{2\lambda}}$');
hold off; box on;
xlim([0 lambda_nd(ind)]);
ax(1) = gca;

hcr(2) = figure('position',[50,0,810,540]);
plot(lambda_nd,real(c2),'linewidth',3,'color',col(2,:),'Handlevisibility','off');
hold on;
plot(lambda_nd,csw,'--k','linewidth',3,'Displayname','$\sqrt{\frac{\lambda}{2}}$');
hold off; box on;
xlim([lambda_nd(ind) max(lambda_nd)-10]);
ax(2) = gca;

hcr(3) = figure('position',[50,0,810,540]);
plot(lambda_nd,real(c),'linewidth',3,'color',col(1,:),'Handlevisibility','off');
hold on;
plot(lambda_nd,-csw,'--k','linewidth',3,'Displayname','$-\sqrt{\frac{\lambda}{2}}$');
hold off; box on;
xlim([lambda_nd(ind) max(lambda_nd)-10]);
ax(3) = gca;

for i = 1:3
    set(ax(i),'XMinorTick','on','YMinorTick','on');
    set(get(ax(i),'xlabel'),'string','$\lambda$','fontsize',30);
    set(get(ax(i),'ylabel'),'string','$c _r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    legend(ax(i));
end