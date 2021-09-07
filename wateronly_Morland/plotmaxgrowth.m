close all; clear;
global expdata erfdata
expdata = load('max_ud2_exp.mat');
erfdata = load('max_ud2_erf.mat');

%% plot growth rate
xf = @(x) x.('delta_list');
yf = @(x) imag(x.('o_list'));
f(1) = plotresult('Morland_growth.bmp',[0,0.8],[1.5,0],xf,yf);
xlabel('$\Delta$');
ylabel('$\omega_i$','rotation',0);
legend;

%% plot phase speed
xf = @(x) x.('delta_list');
yf = @(x) real(x.('c_list'));
f(2) = plotresult('Morland_cr.bmp',[0,0.8],[0.5,0],xf,yf);
xlabel('$\Delta$');
ylabel('$c_r$','rotation',0);
legend('location','southeast');

%% plot lambda
xf = @(x) x.('delta_list');
yf = @(x) x.('lambda_list');
f(3) = plotresult('Morland_lam.bmp',[0,0.8],[2.5,0],xf,yf);
xlabel('$\Delta$');
ylabel('$\lambda$','rotation',0);
legend('location','southeast');

%% Plot lambda with neutral curve
xf = @(x) x.('delta_list');
yf = @(x) x.('lambda_list');
f(3) = plotresult('Morland_lam.bmp',[0,0.8],[2.5,0],xf,yf);
dlt_list = linspace(0,expdata.delta_list(end),5000);
fx = @(x,y) [4*x^2+pi*y 0 -4*pi*y*x^4+4*x^2+2*pi*y 0 pi*y];
hold on;
for dlt = dlt_list
    polylam = fx(expdata.ud_nd,dlt);
    lam = roots(polylam);
    lam = real(lam(imag(lam)<eps));
    plot(dlt,lam,'.','color','#bfbfbf','Markersize',4,'HandleVisibility','off');
end
hold off;
xlabel('$\Delta$');
ylabel('$\lambda$','rotation',0);
legend('location','northwest','FontSize',14);

function f = plotresult(fignam,xL,yL,x,y)
% Read Morland's results
f = figure('position',[150 100 720 640]);
im1 = imread(fignam);
imagesc(xL,yL,im1);
set(gca,'YDir','normal');
% Plot results
global expdata erfdata
expdisp = {'bx','Displayname','exponential','Markersize',8};
erfdisp = {'ro','Displayname','error function'};
hold on;
plot(-1,-1,'-k','Displayname','Morland piecewise-linear');
plot(-1,-1,'--k','Displayname','Morland error function');
plot(x(erfdata),y(erfdata),erfdisp{:});
plot(x(expdata),y(expdata),expdisp{:});
hold off;
grid on; axis square;
end

