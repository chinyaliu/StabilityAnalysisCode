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
%         cs(4) = 0.46+0.008i;
        cs(4) = 0.9778+0.0672i;
    % k = 0.04
    case 0.5^2
        cs(1) = -9.08;
        cs(2) = 0.0852+0.0894i;
        cs(3) = 10.92;
        cs(4) = 0.9742+0.1101i;
    case 5^2
        cs(1) = -10;
        cs(2) = 0.156877+0.08418i;
        cs(3) = 1.92632;
        cs(4) = 0.973055+0.111415i;
end
%% Run solver
col = get(gca,'colororder');
css = nan(length(cs),length(k));
tic;
% figure;
% yline(0,'linewidth',1.5,'color','#898989');
% xlabel('$\omega_r$');
% ylabel('$\omega_i$','rotation',0, 'HorizontalAlignment','right');
% xlim([-0.3 1.5]);
% ylim([-0.15 0.15]);
% xticks(-1:0.5:1.5);
% yticks(-0.15:0.05:0.15);
% hold on;
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
%     g(1) = scatter(real(oall),imag(oall),'ok');
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
%             plot(real(c*k(i)),imag(c*k(i)),'.','Markersize',12,'color',col(j,:));
            cs(j) = c;
            css(j,i) = c;
        end
    end
    
%     titext = sprintf('$k=%.2f$',k(i));
%     title(titext);
%     F(i) = getframe(gcf);
%     delete(g);
end
% hold off;
% close(gcf);
toc;

%%
% writerObj = VideoWriter('reinf_fr5_H5.mp4','MPEG-4');
% writerObj.FrameRate = 40;
% writerObj.Quality = 100;
% open(writerObj);
% for i = 1:length(F)
%     frame = F(i);
%     writeVideo(writerObj, frame);
% end
% close(writerObj);

%% Plot oi vs or
% dm = {'mode 3','mode 1','mode 2','mode 4'};
% col = {"#77AC30",[1 0 0],[0 0 1],[0 0 0]};
dm = {'mode 1','mode 2','mode 3','mode 4'};
css([3 1 2 4],:) = css([1 2 3 4],:);
col([3 1 2 4],:) = col([1 2 3 4],:);
% col([3 1 2 4]) = col([1 2 3 4]);
% figure('position',[50,50,1000,720]);
figure;
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for j = 1:4
%     hh(j) = plot(k,imag(k.*css(j,:)),'linewidth',2, 'DisplayName', dm{j},'color',col(j,:));
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
for i = 1:4
    set(hCopy(i),'XData', NaN', 'YData', NaN);
    hCopy(i).MarkerSize = 15; 
end
hleg = legend(hCopy,'location','east');
title(hleg,'$\textbf{Mode}$','FontSize',20);

%% Plot oi vs k (seperate)
% dm = {'mode 1','mode 2','mode 3'};
figure('position',[50,0,960,1080]);
col = get(gca,'colororder');
for j = 1:3
    subplot(3,1,j);
    yline(0,'linewidth',1.5,'color','#898989');
    hold on;
    hh(j) = plot(k,imag(k.*css(j,:)),'linewidth',2, 'color',col(j,:));
    hold off; grid on;
    xlim([0 4]);
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