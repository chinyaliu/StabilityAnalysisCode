clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,H,k,Fr2,Re,eps,c0,h,~] = pars_wake(3);
zL = (0.74708299+H)*ones(length(k),1);
in_init = {N,H,k(1),h(1),Re,Fr2,bflow};
addvar = struct('zL1',zL(1),'eps',eps);

%% Run solver
tic;
o = NaN(2,length(k));
parfor j = 1:2
    if j == 1
        f = wSubmerged.ddmtype(1);
    else
        f = wSubmerged.ddmtype(45);
    end
    smeth = {alg, de_singularize, do_balancing, eig_spectrum, f, addvar};
    o(j,:) = findmode(k,h,in_init,method,smeth);
end
toc;

%% Plot oi vs k
% fig2 = figure('position',[50,0,720,540]);
fig2 = figure('position',[50,0,810,540]);
% figure;
o(isnan(o))=0;
hold on;
for i = 1:2
    plot(k,imag(o(i,:)),'linewidth',3);
end
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off'); 
hold off; box on; grid on;
ax = gca;
set(ax,'XMinorTick','on','YMinorTick','on');
xlim([0 max(k)]);
ylim([0 0.04]);
yticks(0:0.01:0.04);
xlabel('$k$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

function o = findmode(k,h,in_init,method,smeth)
p1 = wSubmerged(in_init{:});
p1.numMeth(method);
o = NaN(1,length(k));
c0 = min(0.99, 1./sqrt(real(k(1))*in_init{6}));
cutz = NaN(1,length(k));
cutz(1) = -p1.invbf(c0);
smeth{end}.zL1 = cutz(1);
j = 1;
for i = 1:length(k)
    p1.k = k(i); p1.h = h(i);
    o(i) = p1.solver(smeth{:});
    if isnan(p1.zc)
        smeth{end}.zL1 = cutz(j);
    else
        cutz(i)=p1.zc;
        j = i;
        smeth{end}.zL1 = p1.zc;
    end
    fprintf('k = %.2f\n', k(i));
end
end