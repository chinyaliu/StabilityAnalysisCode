clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,k,Fr2,Re,eps,c0,h,f] = pars_wake(3);
% h = 6*ones(1,length(k));
inflec_pt = -0.74708299;
zL = 0.74708299;
cutz = NaN(1,length(k));
cutz(1) = -inflec_pt;
in_init = {N,k(1),h(1),Re,Fr2};
%% Run solver
tic;
p1 = wZhang_ddm(in_init{:});
p1.numMeth(method);
oall = cell(1,3);
R = [1e3,inf];
% R = [1e2,1e3,inf];
for j = 1:length(R)
    fprintf('Re = %4d\n', R(j));
    p1.chgRe(R(j),method);
    o = NaN(1,length(k));
    addvar = struct('zL1',zL,'eps',eps);
    for i = 1:length(k)
        p1.k = k(i); p1.h = h(i);
        o(i) = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
        if isnan(p1.zc)
            cutz(i)=cutz(1);
        else
            cutz(i)=-p1.zc;
        end
        addvar.zL1 = cutz(i);
        fprintf('k = %.2f\n', k(i));
    end
    oall{j} = o;
end
toc;

%% Read Zhang's results & Plot oi vs k
dn = {'100','1000','$\infty$'};
imagedata = imread('Zhang.bmp');
im2 = imbinarize(imagedata).*255;
im2(im2==0) = 100;
im2 = cast(im2,'uint8');
fig = figure('position',[50,0,1000,720]);
imagesc([0,4],[0.04,0],im2);
set(gca,'YDir','normal');
xticks(0:1:4);
yticks(0:0.01:0.04);
hold on;
for i = 1:3
    hh(i) = plot(k,imag(oall{i}),'.','markersize',6, 'DisplayName', dn{i});
end
hold off;
xlabel('$k$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
hCopy = copyobj(hh, ax); 
for i = 1:3
    set(hCopy(i),'XData', NaN', 'YData', NaN);
    hCopy(i).MarkerSize = 15; 
end
hleg = legend(hCopy,'location','east');
title(hleg,'$\textbf{Re}$','FontSize',20);

%% Plot real omega
mk = {'o','+','x'};
ln = {'--','-.',':'};
figure; 
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for i = 1:3
    or = real(oall{i});
    hh2(i) = plot(k,or,ln{i},'linewidth',1.5, 'DisplayName', dn{i});
%     hh2(i) = plot(k,or,mk{i},'markersize',3, 'DisplayName', dn{i});
end
hold off; box on;
xlabel('$k$');
xticks(0:1:4);
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('$\omega_r$,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
hCopy = copyobj(hh2, ax); 
for i = 1:3
    set(hCopy(i),'XData', NaN', 'YData', NaN);
    hCopy(i).MarkerSize = 15; 
end
hleg = legend(hCopy,'location','northwest');
title(hleg,'$\textbf{Re}$','FontSize',20);

%% Plot real c
mk = {'o','+','x'};
ln = {'--','-.',':'};
figure; 
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for i = 1:3
    cr = real(oall{i})./k;
    hh2(i) = plot(k,cr,ln{i},'linewidth',1.5, 'DisplayName', dn{i});
%     hh2(i) = plot(k,or,mk{i},'markersize',3, 'DisplayName', dn{i});
end
hold off;
xlabel('$k$','fontsize',30);
xticks(0:1:4);
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('$c_r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
ax = gca;
hCopy = copyobj(hh2, ax); 
for i = 1:3
    set(hCopy(i),'XData', NaN', 'YData', NaN);
    hCopy(i).MarkerSize = 15; 
end
hleg = legend(hCopy,'location','southeast');
title(hleg,'$\textbf{Re}$','FontSize',20);