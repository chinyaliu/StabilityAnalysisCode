clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,H,k,Fr2,Re,eps,c0,h,f] = pars_wake;

%% Run solver
t1 = tic;
case1 = wSubmerged(N,H,k,h,Re,Fr2,bflow);
case1.numMeth(method);
zL = 0.74708299+H;
[o, an] = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',zL,'eps',eps));
o = o(real(o)>-50);
c = o/k;
zc = case1.zc;
toc(t1);

%% Choose eigenvalue
a = 1:length(c); 
crange = (real(c)>0) & (real(c)<0.5) & (imag(c)<0);
aa = a(crange);
abch = isoutlier(imag(c(aa)),'movmedian',5);
aa = [a(imag(c)>0) aa(abch)];
o_chosen = o(aa); c_chosen = c(aa);
an_c = an(:,aa);

%% Plot eigenvalue spectrum ci_cr
% figure;
% plot(real(c),imag(c),'ok');
% hold on; 
% yline(0,'linewidth',1.5,'color','#898989'); 
% scatter(real(c_chosen),imag(c_chosen),'b','filled');
% hold off;
% xlabel('$c_r$','fontsize',30);
% ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
% ymax = max(abs(imag(c)));
% ylim([-1.5*ymax 1.5*ymax]);
% titext = sprintf('$k=%.2f%+.2fi$',real(k),imag(k));
% title(titext);

%% Find velocity and wave amplitude
modd = 1;
[z, phi] = case1.getprop('modeshape',an_c(:,modd));
if sum(phi(:,1)) < 0
    phi = -phi;
end
up = phi(:,2);
wp = -1i*k*phi(:,1);
qp = an_c(end,1);

x = linspace(0,4*pi,100);
% x_phase = exp(1i*k*x);
T = linspace(0,4*pi,100);
x_phase = exp(1i*T);
upx = real(up*x_phase);
wpx = real(wp*x_phase);
qpx = real(qp*x_phase);

vel_unit = max(sqrt(upx.^2+wpx.^2),[],'all');
% vel_unit = max(max(upx,[],'all'),max(wpx,[],'all'));
% vel_unit = 2*max(qpx);
upx = upx/vel_unit;
wpx = wpx/vel_unit;
qpx = qpx/vel_unit;

%% Plot contour
% plotcont(x,z,upx,qpx,zc,H);
% title('$u''$');
% plotcont(x,z,wpx,qpx,zc,H);
% title('$w''$');

%% Plot velocity field
if H > 0
    yL = -2*H;
else
    yL = -6;
end
[~,zind] = min(abs(z-yL));
ptx = 30;
ptz = 30;
x2 = linspace(0,x(end),ptx);
z2 = -linspace(0,-z(zind),ptz);
[zu,ia] = unique(z,'stable');
xp2 = exp(1i*linspace(0,4*pi,ptx));
u2 = real((interp1(zu,up(ia),z2).')*xp2);
w2 = real((interp1(zu,wp(ia),z2).')*xp2);
[X, Z] = meshgrid(x2,z2);
figure;
quiver(X, Z, u2, w2, 'k', 'linewidth', 1);
hold on; 
plot(x,qpx,'b');
yline(-zc,':','linewidth',1.5,'color','#606060');
yline(case1.zc-2*H, '--', 'linewidth', 1.5,'color','#606060');
yline(-H,'-','linewidth',1.5,'color','#606060');
hold off;
ylim([-5*max(qpx) max(qpx)]);
xlim([0 x(end)]);
set(gca,'XTick',0:pi:4*pi,...
    'XTickLabel',{'0','$\pi$','$2\pi$','$3\pi$','$4\pi$'},...
    'XMinorTick', 'on', 'YMinorTick', 'on');
xlabel('$T$');
ylabel('$z$','rotation',0, 'HorizontalAlignment','right');

%%
figtitle = ["$\phi$", "$\phi_ z$", "$\phi_ {zz}$"];
for i = 1:3
    f = figure('position',[0 0 360 720]);
    plotvar = {abs(phi(:,i)),unwrap(angle(phi(:,i))),real(phi(:,i)),imag(phi(:,i))};
    xline(0,'--','linewidth',1.5,'color','#606060','HandleVisibility','off');
    hold on;
    if ~isnan(case1.zc)
        yline(-case1.zc, ':k','linewidth',1.5,'HandleVisibility','off');
        yline(case1.zc-2*H, ':k','linewidth',1.5,'HandleVisibility','off');
    end
%     arr = -case1.getprop('cut');
%     for k = 2:length(arr)-1
%         yline(arr(k), '--k','linewidth',1.5,'HandleVisibility','off');
%     end
%     yline(-zL, '--b','linewidth',1.5,'HandleVisibility','off');
    plot(real(phi(:,i)),z,'b','DisplayName','real', 'linewidth', 2);
    plot(imag(phi(:,i)),z,'--r','DisplayName','imag', 'linewidth', 2);
    hold off;
    ylabel('$z$','rotation',0, 'HorizontalAlignment','right');
    xticks(0);
    hold off;
    if H > 0
        ylim([-2*H 0]);
    else
        ylim([-6 0]);
    end
    grid on; box on;
    legend('location','southeast');
    title(figtitle(i));
end

function cmap = mycmap
    len =  size(get(gcf,'colormap'),1);
    cmap(1:2*len,:) = [zeros(len,1),zeros(len,1),linspace(147,255,len)';
    linspace(8,230,len)',linspace(8,230,len)',255*ones(len,1)]./255;
    cmap(2*len+1,:) = [1,1,1];
    cmap(2*len+2:4*len+1,:) = [255*ones(len,1),linspace(228,0,len)',linspace(228,0,len)';
    linspace(251,147,len)',zeros(len,1),zeros(len,1)]./255;
end

function plotcont(x,z,vel,q,zc,H)
    [X,Z] = meshgrid(x,z);
    lstep = max(vel,[],'all')/100;
    figure;
    contourf(X,Z,vel,'LineStyle','none','LevelStep',lstep);
    hold on; 
    plot(x,q,'k');
    yline(-zc,'--k','linewidth',1.5);
    hold off;
    if H > 0
        ylim([-2*H 1]);
    else
        ylim([-6 1]);
    end
    xlabel('$\widetilde{x}$');
    ylabel('$\widetilde{z}$','Rotation',0,'HorizontalAlignment','right');
    set(gca,'XTick',0:pi:4*pi,...
        'XTickLabel',{'0','$\pi$','$2\pi$','$3\pi$','$4\pi$'},...
        'XMinorTick', 'on', 'YMinorTick', 'on');
    colormap(mycmap); colorbar; caxis([-1 1]);
end