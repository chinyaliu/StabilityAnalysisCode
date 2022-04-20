% clear all;
ud = 2;
dlt_list = linspace(0,0.8,1000);
fx = @(x,y) [4*x^2+pi*y 0 -4*pi*y*x^4+4*x^2+2*pi*y 0 pi*y];
lam_list = cell(1,length(dlt_list));
for i = 1:length(dlt_list)
    polylam = fx(ud,dlt_list(i));
    lamm = roots(polylam);
    lam_list{i} = real(lamm(imag(lamm)<eps));
end

%% Plot Morland's neutral curve
c0 = linspace(0,1.5,400);
lam_m = 0.5*pi*(1-c0.^4)./c0.^2;
f = figure('position',[150 100 720 640]);
plot(lam_m,c0,'k-'); 
grid on;axis square;
ylim([0 1.5]);
xlim([0 9]);
xlabel('$\tilde{\lambda}/\tilde{\Delta}$');
ylabel('$\tilde{c_0}/\tilde{u_d}\quad$','rotation',0);

%% Plot fast
ud = 2;
dlt_list = linspace(0,0.8,1000);
fx = @(x,y) [4*x^2+pi*y 0 -4*pi*y*x^4+4*x^2+2*pi*y 0 pi*y];
% f = figure('position',[150 100 720 640]);
hold on;
for dlt = dlt_list
    polylam = fx(ud,dlt);
    lamm = roots(polylam);
    lamm = real(lamm(imag(lamm)<eps & real(lamm)>-eps));
    plot(dlt,lamm,'r.','Markersize',4);
end
hold off; 
% grid on;axis square;
ylim([0 2.5]);
xlabel('$\Delta$');
ylabel('$\lambda\quad$','rotation',0);
title(sprintf('$u_d=%1.1f$',ud));