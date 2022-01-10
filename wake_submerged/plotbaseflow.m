% Plot tha base flow of this profile
clear;
[~,~,~,~,~,~,N,H,~,~,Re,~,~,h,~] = pars_wake;
case1 = subRay(N,0,h,struct('H',H,'k',[],'dm',[]));
z = case1.z;
U = case1.U;
inflec_pt = [-H-0.74708299 -H+0.74708299];

%% plot U
figure('position',[200 100 360 720]);
plot(U(:,1),z,'k','linewidth',3);
hold on;
yline(inflec_pt(1), '-b','linewidth',2); 
yline(inflec_pt(2), '-b','linewidth',2); 
yline(-H, '--', 'color', '#46bf12','linewidth',2); 
xline(0, '--', 'color', '#606060');
hold off;
if (h > 6)
    blim = -6-2*H;
else
    blim = fix(-h);
end
ylim([blim 0]);
xlim([0 1.2]);
yticks([]);
xticks([0 1]);
xlabel('$U$');
ylabel('$z$','rotation',0,'HorizontalAlignment','right');

%% plot dU/dz
figure('position',[600 100 360 720]);
plot(U(:,2),z,'k','linewidth',3);
hold on;
yline(inflec_pt(1), '-b','linewidth',2); 
yline(inflec_pt(2), '-b','linewidth',2); 
yline(-H, '--', 'color', '#46bf12','linewidth',2); 
xline(0, '--', 'color', '#606060');
hold off;
if (h > 6)
    blim = -6-2*H;
else
    blim = fix(-h);
end
ylim([blim 0]);
xlim([-0.9 0.1]);
yticks([]);
xticks([0 1]);
xlabel('$U_z$');
ylabel('$z$','rotation',0,'HorizontalAlignment','right');

%% plot (d^2/dz^2)U
figure('position',[1000 100 360 720]);
plot(U(:,3),z,'k','linewidth',3);
hold on;
yline(inflec_pt(1), '-b','linewidth',2); 
yline(inflec_pt(2), '-b','linewidth',2); 
yline(-H, '--', 'color', '#46bf12','linewidth',2); 
xline(0, '--', 'color', '#606060');
hold off;
if (h > 6)
    blim = -6-2*H;
else
    blim = fix(-h);
end
ylim([blim 0]);
xlim([-0.6 1.6]);
yticks([]);
% xticks([0 1]);
xlabel('$U_{zz}$');
ylabel('$z$','rotation',0,'HorizontalAlignment','right');

%% plot (d^3/dz^3)U
c1 = 0.9988; c2 = 0.8814;
U3 = 8*c1*c2^3*tanh(c2*(z+H)).*(sech(c2*(z+H))).^2.*((tanh(c2*(z+H))).^2-2*(sech(c2*(z+H))).^2);
figure('position',[1400 100 360 720]);
plot(U3,z,'k','linewidth',3);
hold on;
yline(inflec_pt(1), '-b','linewidth',2); 
yline(inflec_pt(2), '-b','linewidth',2); 
yline(-H, '--', 'color', '#46bf12','linewidth',2); 
xline(0, '--', 'color', '#606060');
hold off;
if (h > 6)
    blim = -6-2*H;
else
    blim = fix(-h);
end
ylim([blim 0]);
xlim([-0.5 3]);
yticks([]);
% xticks([0 1]);
xlabel('$U_{zzz}$');
ylabel('$z$','rotation',0,'HorizontalAlignment','right');
