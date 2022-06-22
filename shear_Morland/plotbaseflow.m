clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end

%% Plot the base flow of this profile
clear;
[method,~,~,~,~,~,~,ud_nd,delta_nd,~,~,~,~,~,~] = pars_Morland(1,'inv');
h = 3;
z = linspace(-h, 0, 100).';
bflow = ["error function","exponential"];
case1 = wMorland(0,h,ud_nd,delta_nd,0,method,bflow(1),0);
case2 = wMorland(0,h,ud_nd,delta_nd,0,method,bflow(2),0);
U = case1.baseflow(z);
U2 = case2.baseflow(z);

%% plot U
figure('position',[200 100 360 720]);
plot(U(:,1),z,'-k','linewidth',3);
hold on;
plot(U2(:,1),z,':r','linewidth',3);
xline(0, '--', 'color', '#606060','HandleVisibility','off');
hold off;
ylim([-h 0]);
xlim([0 ud_nd]);
xticks([0 ud_nd]);
xticklabels({'0','$U_0$'});
xlabel('$U$');
ylabel('$z$','rotation',0,'HorizontalAlignment','right');
legend(bflow,'fontsize',24,'location','southeast');

%% plot U'
figure('position',[200 100 360 720]);
plot(U(:,2),z,'-k','linewidth',3);
hold on;
plot(U2(:,2),z,':r','linewidth',3);
xline(0, '--', 'color', '#606060','HandleVisibility','off');
hold off;
ylim([-h 0]);
% xlim([0 ud_nd]);
% yticks([]);
% xticks([]);
xlabel('$U_z$');
ylabel('$z$','rotation',0,'HorizontalAlignment','right');
legend(bflow,'fontsize',24,'location','southeast');

%% plot U''
figure('position',[200 100 360 720]);
plot(U(:,3),z,'-k','linewidth',3);
hold on;
plot(U2(:,3),z,':r','linewidth',3);
xline(0, '--', 'color', '#606060','HandleVisibility','off');
hold off;
ylim([-h 0]);
% xlim([0 ud_nd]);
% yticks([]);
% xticks([]);
xlabel('$U_{zz}$');
ylabel('$z$','rotation',0,'HorizontalAlignment','right');
legend(bflow,'fontsize',24,'location','southeast');