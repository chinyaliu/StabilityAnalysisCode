% Plot tha base flow of this profile
clear;
[~,~,~,~,~,~,N,H,~,~,Re,~,~,h,~] = pars_wake;
case1 = subRay(N,0,h,struct('H',H,'k',[],'dm',[]));
z = case1.z;
U = case1.U;
inflec_pt = [-H-0.74708299 -H+0.74708299];

%% plot
figure('position',[200 100 360 720]);
plot(U(:,3),z,'k','linewidth',3);
hold on;
yline(inflec_pt(1),'--r','linewidth',2); 
yline(inflec_pt(2),'--r','linewidth',2); 
yline(-H,'--k','linewidth',2); 
hold off;
ylim([-h 0]);
xlim([0 1.2]);
yticks([]);
xticks([0 1]);
xlabel('$U$');
ylabel('$z$','rotation',0,'HorizontalAlignment','right');