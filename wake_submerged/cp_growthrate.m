clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end
%% Set Solver & Algorithm
numofcase = 2;
[method{1},alg{1},bflow,de_singularize{1},do_balancing{1},~,N{1},H,k,Fr2,Re{1},eps,c0,h,f{1}] = pars_wake(3,'inv');
in_init{1} = {N{1},H,k(1),h(1),Re{1},Fr2,bflow};
[method{2},alg{2},~,de_singularize{2},do_balancing{2},~,N{2},~,~,~,Re{2},~,~,~,f{2}] = pars_wake(3,'vis');
in_init{2} = {N{2},H,k(1),h(1),Re{2},Fr2,bflow};
% [method{2},alg{2},~,de_singularize{2},do_balancing{2},~,N{2},~,~,~,Re{2},~,~,~,f{2}] = pars_wake(3,'inv');
% in_init{2} = in_init{1};
f{2} = wSubmerged.ddmtype(1);
inflec_pt = -0.74708299-H;
zL = (0.74708299+H)*ones(length(k),1);
eig_spectrum = 'max';

%% Run solver
tic;
o = NaN(numofcase,length(k));
cutz = NaN(numofcase,length(k));
c01 = c0(1);
parfor i = 1:numofcase
    p{i} = wSubmerged(in_init{i}{:});
    p{i}.numMeth(method{i});
    addvar = struct('zL1',p{i}.invbf(c01),'eps',eps);
    [o(i,:), cutz(i,:)] = findmode(p{i},k,h,alg{i}, de_singularize{i}, do_balancing{i}, f{i}, addvar);
end
toc;

%% Plot oi vs or
lsy = {'-k','--r','--b'};
fig1 = figure('position',[50,0,1000,720]);hold on;
for i = 1:length(p)
    plot(real(o(i,:)),imag(o(i,:)),lsy{i}, 'Displayname', num2str(Re{i}));
end
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off'); 
hold off; box on; grid on;
xlabel('$\omega _r$');
ylabel('$\omega _i$','rotation',0, 'HorizontalAlignment','right');
leg = legend('location','northeast');
title(leg,'$Re$');

%% Plot oi vs k
o(isnan(o))=0;
fig2 = figure('position',[50,0,1000,720]);hold on;
for i = 1:length(p)
    plot(k,imag(o(i,:)),lsy{i}, 'Displayname', sprintf('$%.0f$',Re{i}));
end
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off; box on; grid on;
xlim([0 max(k)]);
xlabel('$k$');
ylabel('$\omega _i$','rotation',0, 'HorizontalAlignment','right');
leg = legend('location','northeast');
title(leg,'$Re$');
% ylim([0 0.03]);

%% Plot oi vs k
o(isnan(o))=0;
lsy = {'-k',':r','--b'};
nam = {'Original','DDM'};
fig2 = figure('position',[50,0,1000,720]);hold on;
for i = 1:length(p)
    plot(k,imag(o(i,:)),lsy{i}, 'Displayname', nam{i});
end
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off; box on; grid on;
xlim([0 max(k)]);
xlabel('$k$');
ylabel('$\omega _i$','rotation',0, 'HorizontalAlignment','right');
leg = legend('location','northeast');
ylim([0 0.03]);

%% Plot critical height vs k
ind2 = k>0.87;
fig3 = figure('position',[50,0,1000,720]);hold on;
for i = 1:length(p)
    cutemp = cutz(i,:);
    [~,ind] = rmmissing(cutemp);
    plot(k(~(ind&ind2)),-cutemp(~(ind&ind2)),lsy{i}, 'Displayname', num2str(Re{i}));
end
yline(inflec_pt,'linewidth',1.5,'color','#898989','HandleVisibility','off');
xlim([0 max(k)]);
xlabel('$k$');
ylabel('$z_c$','rotation',0, 'HorizontalAlignment','right');
leg = legend('location','northeast');
title(leg,'$Re$');

% %% Plot ci vs cr
% c = o./k;
% fig4 = figure('position',[50,0,1000,720]);
% plot(real(c),imag(c),'k.','Markersize',6);
% hold on; yline(0,'linewidth',1.5,'color','#898989'); hold off;
% xlabel('$c _r$','fontsize',30);
% ylabel('$c _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

function [o, cutz] = findmode(p1,k,h,alg, de_singularize, do_balancing, f, addvar)
    cutz = NaN(1,length(k));
    cutz(1) = addvar.zL1;
    o = NaN(1,length(k));
    j = 1;
    for i = 1:length(k)
        p1.k = k(i); p1.h = h(i);
        o(i) = p1.solver(alg, de_singularize, do_balancing,'max', f, addvar);
        if isnan(p1.zc)
            addvar.zL1 = cutz(j);
        else
            cutz(i)=p1.zc;
            j = i;
            addvar.zL1 = p1.zc;
        end
        fprintf('k = %.2f\n', k(i));
    end
end