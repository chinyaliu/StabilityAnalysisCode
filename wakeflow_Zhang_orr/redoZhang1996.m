clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,k,Fr2,Re,eps,c0,h,f] = pars_wake(3);
% h = 6*ones(2*length(k),1);
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
R = [1e2,1e3,inf];
for j = 1:3
    fprintf('Re = %4d\n', R(j));
    p1.chgRe(R(j),method);
    o = NaN(1,length(k));
    addvar = struct('zL1',zL,'eps',eps);
    if j == 3
        m = 1;
        k2 = [linspace(0.01,1,40) linspace(1.02,3,50) linspace(3.05,4,58)];
        h = 3*2*pi./real(k2);
        o = NaN(1,length(k2));
        for i = 1:length(k2)
            p1.k = k2(i); p1.h = h(i);
            o(i) = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, wZhang_ddm.ddmtype(44), addvar);
            if isnan(p1.zc)
                cutz(i)=cutz(m);
            else
                cutz(i)=-p1.zc;
                m = i;
            end
            addvar.zL1 = cutz(i);
            fprintf('k = %.2f\n', k2(i));
        end
        o = [o(1:90) o(91:3:end)];
    else
        for i = 1:length(k)
            p1.k = k(i); p1.h = h(i);
            o(i) = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, wZhang_ddm.ddmtype(1), addvar);
            fprintf('k = %.2f\n', k(i));
        end
    end
    oall{j} = o;
end
toc;

%% Read Zhang's results & Plot oi vs k
col = [0.16,0.44,1.00;0.93,0.00,0.00;0.00,0.57,0.00];
dn = {'100','1000','$\infty$'};
imagedata = imread('Zhang.bmp');
im2 = imbinarize(imagedata).*255;
im2(im2==0) = 100;
im2 = cast(im2,'uint8');
fig = figure('position',[50,0,720,540]);
imagesc([0,4],[0.04,0],im2);
set(gca,'YDir','normal');
xticks(0:1:4);
yticks(0:0.01:0.04);
hold on;
for i = 1:3
%     hh(i) = plot(k,imag(oall{i}),'.','markersize',10, 'DisplayName', dn{i},'color',col(i,:));
    hh(i) = plot(k,imag(oall{i}),':','linewidth',3, 'DisplayName', dn{i},'color',col(i,:));
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
set(gca,'XMinorTick','on','YMinorTick','on');
ylabel('$\omega_r$','rotation',0, 'HorizontalAlignment','right');
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