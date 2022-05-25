clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
cod = colororder;

%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,k,Fr2,Re,eps,c0,h,f] = pars_wake;
pertE = logspace(-12, -2, 6);
loopnum = 40;

%% Eigenvalue
case1 = wZhang_ddm(N,k,h,Re,Fr2);
case1.numMeth(method);
zL = 0.74708299;
[o, an] = case1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, struct('zL1',zL,'eps',eps));
if isnan(-case1.zc)
    zc = zL;
else
    zc = -case1.zc;
end
c = o(real(o)>-50)/k;
c = c(imag(c)>-5);

%%
tic;
[A,B] = case1.getAB(f,struct('zL1',zL,'eps',eps));
if strcmpi(de_singularize,'y')
    er = -200i;
    sinB = all(B<eps,2);
    B(sinB,:) = A(sinB,:)/er;
end
% BA = B\A;
% KI = k*eye(size(BA,2));
xx = linspace(-1,1.5,50);
yy = linspace(-3,0.2,50);
% xx = linspace(-0.5,1,50);
% yy = linspace(-4,0.5,50);
lny = length(yy);
[X,Y] = meshgrid(xx,yy);
Z = NaN(length(xx),length(yy));
z = complex(X,Y);
parfor i = 1:length(xx)
    for j = 1:lny
%         Z(i,j) = norm(inv(BA-z(i,j)*KI))^(-1);
%         Z(i,j) = norm(inv(z(i,j)*KI-BA))^(-1);
        Z(i,j) = norm(inv(A-z(i,j)*B))^(-1);
%         Z(i,j) = norm(inv(A-z(i,j)*k*B))^(-1);
    end
end
toc;

%%
figure;
scatter(real(c),imag(c),40,'r','filled','HandleVisibility','off');
hold on;
[cc,hh] = contour(X,Y,log10(Z),-20:2:0,'ShowText','on','linewidth',2);
clabel(cc,hh,'FontSize',18,'LabelSpacing',200);
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold off; box on; grid on;
xlim([-0.5 1]);
ylim([-4 0.5]);
aa = colorbar; 
aa.Label.String = '$log(\sigma_\epsilon)$';
aa.Label.Interpreter = 'latex';
caxis([-10 0]);
title(sprintf('$k=%.2f$',real(k)));
xlabel('$c _r$','fontsize',30);
ylabel('$c _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
set(gca,'XMinorTick','on','YMinorTick','on','Fontsize',28);
yticks([-4:1:1]);

%% Pseudoeigenvalue
t1 = tic;
figure;
yline(0,'linewidth',1.5,'color','#898989','HandleVisibility','off');
hold on;
for i = 1:length(pertE)
    cpall = case1.solverE(alg, de_singularize, do_balancing, f, struct('zL1',zc,'eps',eps), pertE(i) ,loopnum);
    [cpBx, cpBy] = findB(c,cpall,pertE(i));
%     scatter(real(cpall),imag(cpall),5,cod(i,:),'filled','HandleVisibility','off');
    plot(cpBx,cpBy,'color',cod(i,:),'linewidth',2,'Displayname',sprintf('%1.0e',pertE(i)));
end
scatter(real(c),imag(c),30,'k','filled','HandleVisibility','off');
hold off;
toc(t1);

%%
xlim([-0.5 1]);
ylim([-4 1]);
xlabel('$c _r$','fontsize',30);
ylabel('$c _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
titext = sprintf('$k=%.2f$',real(k));
title(titext);
legend('location','northeastoutside');

function [cpr,cpi] = findB(c,cp,pertE) 
    % Find boundary
    C1 = repmat(c,1,length(cp));
    C2 = repmat(cp.',length(c),1);
    CC = abs(C1-C2);
    c_min = min(CC,[],1);
    cpall = cp(c_min>pertE*10);
    cpr = real(cpall); cpi = imag(cpall);
    j = boundary(cpr,cpi,1);
    cpr = cpr(j); cpi = cpi(j);
end