clear all;
if ~contains(path,'code_wake;')
    addpath('code_wake');
end 
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,eig_spectrum,N,H,k,~,Re,eps,c0,h,f] = pars_wake;
fr2_list = linspace(0.5,10,100).^2;
inflec_pt = -0.74708299-H;
zL = (0.74708299+H)*ones(length(fr2_list),1);
cutz = -inflec_pt;
addvar = struct('zL1',zL(1),'eps',eps);
in_init = {N,H,k,h,Re,fr2_list(1)};
%% Select specific modes to observe
cs(1) = -9.08;
cs(2) = 0.0852+0.0894i;
cs(3) = 10.92;
cs(4) = 0.9742+0.1101i;
%% Run solver
col = get(gca,'colororder');
css = nan(length(cs),length(fr2_list));
tic;
p1 = wSubmerged(in_init{:});
p1.numMeth(method);
for i = 1:length(fr2_list)
    fprintf('fr = %.2f\n', sqrt(fr2_list(i)));
    p1.Fr2 = fr2_list(i);
    oall = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    oall = oall(real(oall)>-50); % Remove the eigenvalues assigned by de-singularizing
    o = maxeig(oall);
    if real(o) > 0
        addvar.zL1=p1.criticalH(real(o)/k);
    end
    call = oall/k;
    
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
            cs(j) = c;
            css(j,i) = c;
        end
    end
end
toc;

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
% xlim([-0.5 1.5]);
% ylim([-0.05 0.15]);
ax = gca;
hCopy = copyobj(hh, ax); 
for i = 1:4
    set(hCopy(i),'XData', NaN', 'YData', NaN);
    hCopy(i).MarkerSize = 15; 
end
hleg = legend(hCopy,'location','east');
title(hleg,'$\textbf{Mode}$','FontSize',20);

%% Plot oi vs fr (seperate)
% dm = {'mode 1','mode 2','mode 3'};
fr = (fr2_list).^(0.5);
figure('position',[50,0,960,1080]);
col = get(gca,'colororder');
for j = 1:4
    subplot(4,1,j);
    yline(0,'linewidth',1.5,'color','#898989');
    hold on;
    hh(j) = plot(fr,imag(k.*css(j,:)),'linewidth',2, 'color',col(j,:));
    hold off; grid on;
    xlim([0.5 2.5]);
%     ylim([-0.01 0.03]);
    if j == 4
        xlabel('$Fr$');
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