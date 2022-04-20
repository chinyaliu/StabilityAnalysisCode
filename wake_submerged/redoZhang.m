clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,~,N,H,k,Fr2,Re,eps,c0,h,f] = pars_wake(3);
inflec_pt = -0.74708299;
zL = -inflec_pt;
in_init = {N,H,k(1),h(1),Re,Fr2,bflow};
addvar = struct('zL1',zL,'eps',eps);
smeth = {alg, de_singularize, do_balancing, 'all', f, addvar};
%% Select specific modes to observe
% 100
cs(1,1) = -5.69935554321657 + 0.155435654161913i;
cs(1,2) = 0.0943324372553924 - 0.359774643534542i;
cs(1,3) = 7.64805342765973 - 0.0885146263928389i;
% 1000
cs(2,1) = -5.69042132200112 + 0.0157451056762317i;
cs(2,2) = -0.0903635749346843 - 0.0288947409210721i;
cs(2,3) = 7.64522803387940 - 0.00878703878220784i;
% inf
cs(3,1) = -5.69036042186261;
cs(3,2) = 0.0246719832072093 + 0.0344079342992104i;
cs(3,3) = 7.64521877450795;
%% Run solver
tic;
% p1 = wSubmerged(in_init{:});
% p1.numMeth(method);
% oall = cell(1,3);
R = [1e2,1e3,inf];
% css = cell(3,1);
% oss = cell(3,1);
k1 = [linspace(0.01,0.2,20) linspace(0.21,4,95)];
% k1 = [linspace(0.01,0.2,20) linspace(0.21,1,50) linspace(1.01,4,100)];
% parfor j = 1:length(R)
parfor j = 1:2
    in_init = {N,H,k(1),h(1),Re,Fr2,bflow};
    addvar = struct('zL1',zL,'eps',eps);
    smeth = {alg, de_singularize, do_balancing, 'all', f, addvar};
    in_init{5} = R(j);
    if isinf(R(j))
        [css{j},oss{j}] = findmodeRE(k1,h,in_init,method,smeth,cs(j,:));
    else
        [css{j},oss{j}] = findmode(k,h,in_init,method,smeth,cs(j,:));
    end
end
toc;

%% parameters
mk = {'o','+','x'};
ln = {'-','--',':'};
% col = cbrewer('qual','Set1',3);
col = [0.16,0.44,1.00;0.93,0.00,0.00;0.00,0.57,0.00];
dn = {'100','1000','$\infty$'};

%% Read Zhang's results & Plot oi vs k
imagedata = imread('Zhang.bmp');
im2 = imbinarize(imagedata).*255;
im2(im2==0) = 100;
im2 = cast(im2,'uint8');
fig = figure('position',[50,0,720,540]);
% fig = figure('position',[50,0,1000,720]);
imagesc([0,4],[0.04,0],im2);
set(gca,'YDir','normal');
xticks(0:1:4);
yticks(0:0.01:0.04);
hold on;
for i = 1:3
    [~,ind] = max(imag(oss{i}),[],1,'linear');
    osmax = oss{i}(ind);
    if i == 3
%         plot(k1,imag(osmax),':','linewidth',3,'color',col(i,:),'Displayname',dn{i});
        plot(k1,imag(osmax),'.','markersize',10,'color',col(i,:),'Displayname',dn{i});
    else
%         plot(k,imag(osmax),':','linewidth',3,'color',col(i,:),'Displayname',dn{i});
        plot(k,imag(osmax),'.','markersize',10,'color',col(i,:),'Displayname',dn{i});
    end
end
hold off;
xlabel('$k$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
lg = legend('location','east','fontsize',24);
title(lg,'$Re$')

%% Read Zhang's results & Plot oi vs k (by mode)
imagedata = imread('Zhang.bmp');
im2 = imbinarize(imagedata).*255;
im2(im2==0) = 100;
im2 = cast(im2,'uint8');
fig = figure('position',[50,0,720,540]);
% fig = figure('position',[50,0,1000,720]);
imagesc([0,4],[0.04,0],im2);
set(gca,'YDir','normal');
xticks(0:1:4);
yticks(0:0.01:0.04);
hold on;
for i = 1:3
    [~,ind] = max(imag(oss{i}),[],1,'linear');
    osmax = oss{i}(ind);
    if i == 3
        for j = 1:3
        	plot(k1,imag(oss{i}(j,:)),ln{i},'linewidth',3,'color',col(j,:),'Handlevisibility','off');
        end
    else
        for j = 1:3
        	plot(k,imag(oss{i}(j,:)),ln{i},'linewidth',3,'color',col(j,:),'Handlevisibility','off');
        end
    end
end
plot(nan,nan,'.','markersize',25,'color',col(2,:),'displayname','mode 1');
plot(nan,nan,'.','markersize',25,'color',col(3,:),'displayname','mode 2');
plot(nan,nan,'.','markersize',25,'color',col(1,:),'displayname','mode 3');
hold off;
xlabel('$k$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
lg = legend('location','east','fontsize',24);

%% Plot o vs k
for j = 1:3
    fig = figure('position',[50,50,1620,540]);
    subplot(1,2,1);
    yline(0,'linewidth',1.5,'color','#898989','Handlevisibility','off');
    hold on;
    for i = 1:2
        plot(k,real(oss{i}(j,:)),ln{i},'linewidth',3,'color',col(j,:),'Displayname',dn{i});
    end
    plot(k1,real(oss{3}(j,:)),ln{3},'linewidth',3,'color',col(j,:),'Displayname',dn{3});
    hold off; box on;
    xlabel('$k$','fontsize',30);
    ylabel('$\omega _r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    xlim([0 4]);
    xticks(0:1:4);
    ylim([-1.5 1.5]);
    yticks(-1.5:0.5:1.5);
    ax = gca;
    set(ax,'XMinorTick','on','YMinorTick','on');
    lg = legend('location','southeast','fontsize',24);
    title(lg,'$Re$')
    
    subplot(1,2,2);
    yline(0,'linewidth',1.5,'color','#898989','Handlevisibility','off');
    hold on;
    for i = 1:2
        plot(k,imag(oss{i}(j,:)),ln{i},'linewidth',3,'color',col(j,:),'Displayname',dn{i});
    end
    plot(k1,imag(oss{3}(j,:)),ln{3},'linewidth',3,'color',col(j,:),'Displayname',dn{3});
    hold off; box on;
    xlabel('$k$','fontsize',30);
    ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    xlim([0 4]);
    ylim([-0.03 0.03]);
    xticks(0:1:4);
    yticks(-0.03:0.01:0.03);
    ax = gca;
    set(ax,'XMinorTick','on','YMinorTick','on');
    ax.YAxis.Exponent = -2;
    lg = legend('location','northeast','fontsize',24);
    title(lg,'$Re$')
end

%% Plot c real vs imag
figure; 
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for i = 1:3
    for j = 1:3
        plot(real(css{i}(j,:)),imag(css{i}(j,:)),'.','Markersize',12,'color',col(j,:));
    end
end
hold off; box on;
xlabel('$c_r$','fontsize',30);
xticks(0:1:4);
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('$c_i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

%% Plot real c
figure; 
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for i = 1:3
    for j = 1:3
        plot(k,real(css{i}(j,:)),ln{i},'linewidth',1.5,'color',col(j,:));
    end
end
hold off;
xlabel('$k$','fontsize',30);
xticks(0:1:4);
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('$c_r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

function o = maxeig(ev)
    a = 1:length(ev);
    ev_ind = a(abs(ev)<1e+3 & abs(ev)>1e-5 & abs(imag(ev)) > 1e-10);
    [~, ind] = sort(imag(ev(ev_ind)),'descend');
    ev_ind = ev_ind(ind);
    w = ev(ev_ind);
    if isempty(w)
        o = NaN;
    else
        o = w(1);
    end
end

function [css,oss] = findmode(k,h,in_init,method,smeth,cs)
    p1 = wSubmerged(in_init{:});
    p1.numMeth(method);
    cutz = smeth{end}.zL1;
    for i = 1:length(k)
        p1.k = k(i); p1.h = h(i);
        oall = p1.solver(smeth{:});
        o = maxeig(oall);
        if (real(o) > 0) && (real(o/k(i)) < 1)
            cutz = p1.invbf(real(o)/k(i));
        end
        smeth{end}.zL1 = cutz;
        
        call = oall/k(i);
        if (i~=1)
            cmatn = repmat(calln,1,length(call));
            cmat = repmat(call.',length(calln),1);
            [~,ind] = min(abs(cmatn-cmat),[],2);
            cnind = call(setdiff(1:end,ind));
            calln = call(ind);
        else
            cnind = call(imag(call)>-2);
            c1 = repmat(cs.',1,length(cnind));
            c2 = repmat(cnind.',length(cs),1);
            [~,ind] = min(abs(c1-c2),[],2);
            calln = cnind(setdiff(1:end,ind));
        end

        % Choose selected eigenvalues
        if ~isempty(cnind)
            for m = 1:size(cs,2)
                [~,ind] = min(abs(cnind-cs(m)));
                c = cnind(ind);
                [~,ind2] = min(abs(c-cs));
                if ind2 == m
                    cs(m) = c;
                    css(m,i) = c;
                    oss(m,i) = c*k(i);
                end
            end
        end
    end
end

function [css,oss] = findmodeRE(k,h,in_init,method,smeth,cs)
    smeth{5} = wSubmerged.ddmtype(44);
    p1 = wSubmerged(in_init{:});
    p1.numMeth(method);
    cutz = smeth{end}.zL1;
    p1.N = 600;
    for i = 1:length(k)
        p1.k = k(i); p1.h = h(i);
        oall = p1.solver_RE(smeth{:});
        o = oall(1);
        if (real(o) > 0) && (real(o/k(i)) < 1)
            cutz = p1.invbf(real(o)/k(i));
        end
        smeth{end}.zL1 = cutz;
        cnind = oall/k(i);

        % Choose selected eigenvalues
        if ~isempty(cnind)
            for m = 1:size(cs,2)
                [~,ind] = min(abs(cnind-cs(m)));
                c = cnind(ind);
                [~,ind2] = min(abs(c-cs));
                if ind2 == m
                    cs(m) = c;
                    css(m,i) = c;
                    oss(m,i) = c*k(i);
                end
            end
        end
    end
end