clear all;
if ~contains(path,'code_Morland;')
    addpath('code_Morland');
end
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,ud_nd,~,~,~,~,f,epss] = pars_Morland;
dt = linspace(0.01,0.8,40);
lnx = length(dt);
lam = linspace(0.01,2.5,40);
lny = length(lam);

%% Run solver
c = nan(lnx,lny);
k = repmat(2*pi./lam,lnx,1);
tic;
parfor j = 1:lny
    h =2*lam(j);
    c0 = min(pi+pi./lam(j)^2, 0.99*ud_nd);
    dtj = dt;
    for i = 1:lnx
        case1 = wMorland(N,h,ud_nd,dtj(i),lam(j),method,bflow);
        addvar = struct('zL1',case1.invbf(c0),'eps',epss);
        c1 = case1.solver_RE(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        c(i,j) = c1(1);
    end
end
o = c.*k;
toc;

% %% Plot growth rate
[X,Y] = meshgrid(dt,lam);
% f = figure;
% contourf(X,Y,imag(o).',20,'LineColor','none')
% [mo, I] = max(imag(o),[],2);
% indi = (I>2)&(mo>1e-4);
% hold on;
% plot(dt(indi),lam(I(indi)),'ko','markersize',6);
% pltneutral(ud_nd,linspace(dt(1),dt(end),1000),f);
% hold off;
% xlim([0 0.8]);
% ylim([0 2.5]); 

% %% Plot phase velocity
% figure;
% rc = real(c);
% rc(imag(o)<eps) = 0;
% contourf(X,Y,rc.',20,'LineColor','none')

%%
lamlim = findneutral(ud_nd,dt);
ustab = false(lnx,lny);
for i = 1:lnx
    if length(lamlim{i}) > 1
        maxlam = max(lamlim{i});
        minlam = min(lamlim{i});
        ulam = lam > minlam & lam < maxlam;
        ustab(i,ulam) = true;
    end
end
rcc = real(c);
rcc(imag(o)<eps) = 0;
rcc(imag(o)<eps & ustab) = nan;
F = fillmissing(rcc,'linear',1,'SamplePoints',lam);
F(F<=0) = nan;
% figure;
% contourf(X,Y,F.',20,'LineColor','none');
dt2 = linspace(0.01,0.8,100);
lam2 = linspace(0.01,2.5,100);
[X2,Y2] = meshgrid(dt2,lam2);
F2 = interp2(X,Y,F,X2,Y2);
lnx2 = length(dt2);
lny2 = length(lam2);
zL1 = nan(lnx2,lny2);
for i = 1:lnx2
    for j = 1:lny2
        if ~isnan(F2(i,j))
            h =2*lam2(j);
            case1 = wMorland(N,h,ud_nd,dt2(i),lam2(j),method,bflow);
            zL1(i,j) = case1.invbf(F2(i,j));
        end
    end
end

%% Run solver with given initial guess
c2 = nan(lnx2,lny2);
k2 = repmat(2*pi./lam2,lnx2,1);
tic;
parfor j = 1:lny2
    h =2*lam2(j);
    dtj = dt2;
    for i = 1:lnx2
        if ~isnan(zL1(i,j))
            case1 = wMorland(N,h,ud_nd,dtj(i),lam2(j),method,bflow);
            addvar = struct('zL1',zL1(i,j),'eps',epss);
            c1 = case1.solver_RE(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
            c2(i,j) = c1(1);
        end
    end
end
o2 = c2.*k2;
toc;

%% Plot growth rate
f = figure;
% contourf(X2,Y2,imag(o2).',20,'LineColor','none');
pcolor(dt2,lam2,ioo2.');
shading flat;
colorbar;
[mo, I] = max(imag(o2),[],2);
indi = (I>2)&(mo>1e-4);
hold on;
% plot(dt2(indi),lam2(I(indi)),'-k.','markersize',20);
pltneutral(ud_nd,linspace(dt2(1),dt2(end),4000),f);
hold off;
xlim([0 0.8]);
ylim([0 2.5]); 
xlabel('$\Delta$');
ylabel('$\lambda$','rotation',0, 'HorizontalAlignment','right');