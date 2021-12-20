clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,~,N,H,k,Fr2,Re,eps,c0,h,f] = pars_wake(3);
eig_spectrum = 'all';
inflec_pt = -0.74708299-H;
zL = (0.74708299+H)*ones(length(k),1);
cutz = -inflec_pt;
addvar = struct('zL1',zL(1),'eps',eps);
in_init = {N,H,k(1),h(1),Re,Fr2};
%% Select specific modes to observe
switch(Fr2)
    % k = 0.05
    case 1.5^2
        cs(1) = -2.5;
        cs(2) = 0.099+0.096i;
        cs(3) = 4;
        cs(4) = 0.9778+0.0672i;
%     % k = 0.04
%     case 0.5^2
%         cs(1) = -9.08;
%         cs(2) = 0.0852+0.0894i;
%         cs(3) = 10.92;
%         cs(4) = 0.9742+0.1101i;
    % k = 0.05, Re = 1000
    case 0.5^2
        cs(1) = -8;
        cs(2) = 0.065+0.046i;
        cs(3) = 10;
        cs(4) = 0.9933+0.1053i;
    % k = 0.04 H = 2
    case 2^2
        cs(1) = -1.678;
        cs(2) = 0.1129+0.1018i;
        cs(3) = 3.36;
        cs(4) = 0.8769+0.00313i;
    case 3^2
        cs(1) = -0.8;
        cs(2) = 0.0117+0.099i;
        cs(3) = 2;
        cs(4) = 0.978+0.068i;
    case 5^2
        cs(1) = -10;
        cs(2) = 0.156877+0.08418i;
        cs(3) = 1.92632;
        cs(4) = 0.973055+0.111415i;
    % k = 0.04, H = 1.16, Re = 180
    case 0.3^2
        cs(1) = -15.8-0.003i;
        cs(2) = 0.066-0.049i;
        cs(3) = 17;
        cs(4) = 0.635-0.254i;
%     % k = 0.04, H = 0
%     case 10^2
%         cs(1) = -0.05;
%         cs(2) = 1.41;
%     % k = 0.01, H = 10
%     case 10^2
%         cs(1) = -0.157;
%         cs(2) = 0.1095+0.0645i;
%         cs(3) = 1.966;
%         cs(4) = 0.9965+0.0307i;
%     % k = 0.05, H = 5, Re=1000
    case 10^2
        cs(1) = -0.072+0.0253i;
        cs(2) = 0.558-0.02i;
        cs(3) = 1.316;
        cs(4) = 0.988+0.122i;
%     % k = 0.05, H = 5, Re=1000
%     case 2^2
%         cs(1) = -1.4;
%         cs(2) = 0.082+0.046i;
%         cs(3) = 3.11;
%         cs(4) = 0.9932+0.1059i;
    otherwise
        error('Initial mode value not given.');
end
%% Run solver
figure;
yline(0,'linewidth',1.5,'color','#898989');
xlabel('$\omega_r$');
ylabel('$\omega_i$','rotation',0, 'HorizontalAlignment','right');
xlim([-0.3 1.5]);
ylim([-0.15 0.25]);
xticks(-1:0.5:1.5);
yticks(-0.15:0.05:0.25);
hold on;

col = get(gca,'colororder');
css = nan(length(cs),length(k));
tic;
p1 = wSubmerged(in_init{:});
p1.numMeth(method);
for i = 1:length(k)
    fprintf('k = %.2f\n', k(i));
    p1.k = k(i); p1.h = h(i);
    oall = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    oall = oall(real(oall)>-50); % Remove the eigenvalues assigned by de-singularizing
    o = maxeig(oall);
    if real(o) > 0
        addvar.zL1=p1.criticalH(real(o)/k(i));
    end

    % Plot the eigenspectrum
    g(1) = scatter(real(oall),imag(oall),'ok');
    call = oall/k(i);
    
    if isinf(Re)
        a = 1:length(call); 
        crange = ((real(call)-0.0012>-1e-5) & (real(call)-1<=1e-5));
        dis = ((real(call)-1).^2 +imag(call).^2)>1e-5;
        aa = a(crange&dis);
        abch = isoutlier(imag(call(aa)),'movmedian',20);
        aa = [a(dis&~crange) aa(abch)];
        cnind = call(aa);
    else
        if (i~=1) && (abs(k(i)-0.8)>1e-5)
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
    end
    
    % Plot selected eigenvalues
    for j = 1:length(cs)
        [mc,ind] = min(abs(cnind-cs(j)));
        c = cnind(ind);
        [~,ind2] = min(abs(c-cs));
        if ind2 == j
            plot(real(c*k(i)),imag(c*k(i)),'.','Markersize',12,'color',col(j,:));
            cs(j) = c;
            css(j,i) = c;
        end
    end
    
    titext = sprintf('$k=%.2f$',k(i));
    title(titext);
    F(i) = getframe(gcf);
    delete(g);
end
hold off;
close(gcf);
toc;

%%
writerObj = VideoWriter('reinf_fr15_H08.mp4','MPEG-4');
writerObj.FrameRate = 40;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);

%% Plot oi vs or
dm = {'mode 1','mode 2','mode 3','mode 4'};
col([3 1 2 4],:) = col([1 2 3 4],:);
css([3 1 2 4],:) = css([1 2 3 4],:);
figure;
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for j = 1:length(cs)
%     hh(j) = plot(real(css(j,:)),imag(css(j,:)),'o','markersize',5,'linewidth',2, 'DisplayName', dm{j},'color',col(j,:));
    hh(j) = plot(real(k.*css(j,:)),imag(k.*css(j,:)),'o','markersize',5,'linewidth',2, 'DisplayName', dm{j},'color',col(j,:));
end
hold off; grid on;
% xlabel('$k$','fontsize',30);
xlabel('$\omega _r$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
xlim([-0.5 1.5]);
ylim([-0.05 0.15]);
ax = gca;
hCopy = copyobj(hh, ax); 
for i = 1:length(cs)
    set(hCopy(i),'XData', NaN', 'YData', NaN);
    hCopy(i).MarkerSize = 15; 
end
hleg = legend(hCopy,'location','east');
title(hleg,'$\textbf{Mode}$','FontSize',20);

%% Plot oi vs k
% dm = {'mode 3','mode 1','mode 2','mode 4'};
% col = {"#77AC30",[1 0 0],[0 0 1],[0 0 0]};
figure;
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for j = 1:length(cs)
    hh(j) = plot(k,imag(k.*css(j,:)),'o','markersize',5, 'DisplayName', dm{j},'color',col(j,:));
%     hh(j) = plot(real(k.*css(j,:)),imag(k.*css(j,:)),'o','markersize',5,'linewidth',2, 'DisplayName', dm{j},'color',col(j,:));
end
hold off; grid on;
% xlabel('$k$','fontsize',30);
xlabel('$k$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
xlim([-0.5 2.5]);
ylim([-0.05 0.15]);
ax = gca;
hCopy = copyobj(hh, ax); 
for i = 1:length(cs)
    set(hCopy(i),'XData', NaN', 'YData', NaN);
    hCopy(i).MarkerSize = 15; 
end
hleg = legend(hCopy,'location','east');
title(hleg,'$\textbf{Mode}$','FontSize',20);

%% Plot oi vs k (seperate)
% dm = {'mode 1','mode 2','mode 3'};
figure('position',[50,0,960,1080]);
col = get(gca,'colororder');
for j = 1:length(cs)
    subplot(length(cs),1,j);
    yline(0,'linewidth',1.5,'color','#898989');
    hold on;
    hh(j) = plot(k,imag(k.*css(j,:)),'linewidth',2, 'color',col(j,:));
    hold off; grid on;
    xlim([0 2.5]);
    ylim([-0.01 0.03]);
    if j == 3
        xlabel('$k$');
    end
    ylabel('$\omega _i$','rotation',0, 'HorizontalAlignment','right');
    title(dm{j},'FontSize',20);
end

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