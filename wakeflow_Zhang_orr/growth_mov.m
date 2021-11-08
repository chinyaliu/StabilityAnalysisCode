clear all;
if ~contains(path,'./code_wake;')
    addpath('./code_wake');
end 
%% Set Solver & Algorithm
[method,alg,bflow,de_singularize,do_balancing,~,N,k,Fr2,Re,eps,c0,h,f] = pars_wake(3);
eig_spectrum = 'all';
inflec_pt = -0.74708299;
zL = 0.74708299*ones(length(k),1);
cutz = NaN(1,length(k));
cutz(1) = -inflec_pt;
addvar = struct('zL1',zL(1),'eps',eps);
in_init = {N,k(1),h(1),Re,Fr2};
%% Select specific modes to observe
switch(Re)
%     % k = 0.1
%     case 100
%         cs(1) = -1.35243512853880 + 0.0741596105496452i;
%         cs(2) = 0.0964176159568175 - 0.106243617763662i;
%         cs(3) = 2.92650940472364 - 0.0194903151119640i;
%     case 1000
%         cs(1) = -1.34599172404507 + 0.00743722414759224i;
%         cs(2) = 0.123010290800953 + 0.0623048508068197i;
%         cs(3) = 2.92617260767161 - 0.00193739708710329i;
%     case inf
%         cs(1) = -1.34595057448545;
%         cs(2) = 0.137230193360875 + 0.108602511575087i;
%         cs(3) = 2.92617533571192;
    % k = 0.01
    case 100
        cs(1) = -5.69935554321657 + 0.155435654161913i;
        cs(2) = 0.0943324372553924 - 0.359774643534542i;
        cs(3) = 7.64805342765973 - 0.0885146263928389i;
    case 1000
        cs(1) = -5.69042132200112 + 0.0157451056762317i;
        cs(2) = -0.0903635749346843 - 0.0288947409210721i;
        cs(3) = 7.64522803387940 - 0.00878703878220784i;
    case inf
        cs(1) = -5.69036042186261;
        cs(2) = 0.0246719832072093 + 0.0344079342992104i;
        cs(3) = 7.64521877450795;
end
%% Run solver
col = get(gca,'colororder');
css = nan(3,length(k));
tic;
yline(0,'linewidth',1.5,'color','#898989');
% xlabel('$c_r$');
% ylabel('$c_i$','rotation',0, 'HorizontalAlignment','right');
% xlim([-6 8]);
% ylim([-0.15 0.15]);

xlabel('$\omega_r$');
ylabel('$\omega_i$','rotation',0, 'HorizontalAlignment','right');
xlim([-0.3 1.5]);
ylim([-0.15 0.15]);
xticks(-1:0.5:1.5);
yticks(-0.15:0.05:0.15);
hold on;
p1 = wZhang_ddm(in_init{:});
p1.numMeth(method);
for i = 1:length(k)
    fprintf('k = %.2f\n', k(i));
    p1.k = k(i); p1.h = h(i);
    oall = p1.solver(alg, de_singularize, do_balancing, eig_spectrum, f, addvar);
    o = maxeig(oall);
    if real(o) > 0
        cutz(i)=-p1.criticalH(real(o)/k(i));
    else
        cutz(i)=cutz(1);
    end
    addvar.zL1 = cutz(i);

    % Plot the eigenspectrum
    call = oall/k(i);
%     g(1) = scatter(real(call),imag(call),'ok');
    g(1) = scatter(real(oall),imag(oall),'ok');
    
    if isinf(Re)
        a = 1:length(call); 
        crange = ((real(call)-0.0012>-1e-5) & (real(call)-1<=1e-5));
        dis = ((real(call)-1).^2 +imag(call).^2)>1e-5;
        aa = a(crange&dis);
        abch = isoutlier(imag(call(aa)),'movmedian',5);
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
%         if ~(((mc>0.1) && (length(cnind)<4))|| isnan(cs(j)))
        if ind2 == j
%             plot(real(c),imag(c),'.','Markersize',12,'color',col(j,:));
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
writerObj = VideoWriter('reinf.mp4','MPEG-4');
writerObj.FrameRate = 40;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);

%% Plot oi vs k
dm = {'mode 3','mode 1','mode 2'};
% figure('position',[50,50,1000,720]);
figure;
col = get(gca,'colororder');
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for j = 1:3
    hh(j) = plot(k,imag(k.*css(j,:)),'linewidth',2, 'DisplayName', dm{j},'color',col(j,:));
end
hold off; grid on;
xlabel('$k$','fontsize',30);
ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
%     xlim([0 4]);
%     ylim([0 0.03]);
ax = gca;
hCopy = copyobj(hh, ax); 
for i = 1:3
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

%% Plot oi vs or
dm = {'mode 1','mode 2','mode 3'};
css([3 1 2],:) = css([1 2 3],:);
col([3 1 2],:) = col([1 2 3],:);
figure;
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for j = 1:3
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
for i = 1:3
    set(hCopy(i),'XData', NaN', 'YData', NaN);
    hCopy(i).MarkerSize = 15; 
end
hleg = legend(hCopy,'location','east');
title(hleg,'$\textbf{Mode}$','FontSize',20);

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