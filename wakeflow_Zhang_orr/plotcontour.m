clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,k,Fr2,Re,eps,c0,h,f] = pars_wake;

%% Run solver
t1 = tic;
case1 = wZhang_ddm(N,k,h,Re,Fr2);
case1.numMeth(method);
zL = 0.74708299;
[o, an] = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',zL,'eps',eps));
c = o/k;
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
[z, phi] = case1.findmodeshape(an_c(:,1));
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

vel_unit = max(max(upx,[],'all'),max(wpx,[],'all'));
upx = upx/vel_unit;
wpx = wpx/vel_unit;
qpx = qpx/vel_unit;

%% Plot contour
zc = case1.zc;
plotcont(x,z,upx,qpx,zc);
title('$u''$');
plotcont(x,z,wpx,qpx,zc);
title('$w''$');

function cmap = mycmap
len =  size(get(gcf,'colormap'),1);
cmap(1:2*len,:) = [zeros(len,1),zeros(len,1),linspace(147,255,len)';
    linspace(8,230,len)',linspace(8,230,len)',255*ones(len,1)]./255;
cmap(2*len+1,:) = [1,1,1];
cmap(2*len+2:4*len+1,:) = [255*ones(len,1),linspace(228,0,len)',linspace(228,0,len)';
    linspace(251,147,len)',zeros(len,1),zeros(len,1)]./255;
end

function plotcont(x,z,vel,q,zc)
[X,Z] = meshgrid(x,z);
figure;
contourf(X,Z,vel,'LineStyle','none','LevelStep',0.005);
hold on; 
plot(x,q/max(q),'k');
yline(zc,'--k','linewidth',1.5);
hold off;
ylim([-6 1]);
xlabel('$\widetilde{x}$');
ylabel('$\widetilde{z}$','Rotation',0,'HorizontalAlignment','right');
set(gca,'XTick',0:pi:4*pi,...
    'XTickLabel',{'0','$\pi$','$2\pi$','$3\pi$','$4\pi$'},...
	'XMinorTick', 'on', 'YMinorTick', 'on');
colormap(mycmap); colorbar; caxis([-1 1]);
end