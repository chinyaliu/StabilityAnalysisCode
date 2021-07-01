close all;
%% Load results from redoZhang1996_2
% load('save06102114.mat');
%% Parameters
mk = {'o','+','x'};
ln = {'-','-.','--'};
h1 = figure;
col = get(gca,'colororder');
close(h1);
dn = {'100','1000','$\infty$'};
dm = {'mode 3','mode 1','mode 2'};
%% Read Zhang's results & Plot oi vs k
figure('position',[50,50,1000,720]);
imagedata = imread('Zhang.bmp');
im2 = imbinarize(imagedata).*255;
im2(im2==0) = 100;
im2 = cast(im2,'uint8');
imagesc([0,4],[0.04,0],im2);
set(gca,'YDir','normal');
xticks(0:1:4);
yticks(0:0.01:0.04);
hold on;
for i = 1:3
    [~,ind] = max(imag(oss{i}),[],1,'linear');
    osmax = oss{i}(ind);
    hh(i) = plot(k,imag(osmax),'.','markersize',8, 'DisplayName', dn{i});
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

%% Read Zhang's results & Plot oi vs k (Mode)
figure('position',[50,50,1000,720]);
imagedata = imread('Zhang.bmp');
im2 = imbinarize(imagedata).*255;
im2(im2==0) = 100;
im2 = cast(im2,'uint8');
imagesc([0,4],[0.04,0],im2);
set(gca,'YDir','normal');
xticks(0:1:4);
yticks(0:0.01:0.04);
hold on;
for i = 1:3
    for j = 1:3
        if i == 1
            hh(j) = plot(k,imag(oss{i}(j,:)),'.','markersize',8, 'DisplayName', dm{j},'color',col(j,:));
        else
            plot(k,imag(oss{i}(j,:)),'.','markersize',8, 'HandleVisibility', 'off','color',col(j,:));
        end
    end
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
title(hleg,'$\textbf{Mode}$','FontSize',20);

%% Plot ci vs cr (Mode)
for i = 1:3
    figure;
    yline(0,'linewidth',1.5,'color','#898989');
    hold on;
    for j = 1:3
        hh(j) = plot(real(css{i}(j,:)),imag(css{i}(j,:)),'.','markersize',8, 'DisplayName', dm{j},'color',col(j,:));
    end
    hold off;
    xlabel('$c _r$','fontsize',30);
    ylabel('$c _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
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
end

%% Plot oi & or vs k
for i = 1:3
    fig(i) = figure('position',[50,50,1800,720]);
    
    subplot('position',[0.1 0.15 0.35 0.7]);
    yline(0,'linewidth',1.5,'color','#898989', 'HandleVisibility', 'off');
    hold on;
    for j = 1:3
        hh(j) = plot(k,imag(oss{j}(i,:)),ln{j},'linewidth',2, 'DisplayName', dn{j});
    end
    hold off;
    xlabel('$k$','fontsize',30);
    ylabel('$\omega _i$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    xlim([0 4]);
    ylim([-0.08 0.03]);
    xticks(0:1:4);
    yticks(-0.08:0.02:0.03);
    legend('location','northeast');
    
    subplot('position',[0.55 0.15 0.35 0.7]);
    yline(0,'linewidth',1.5,'color','#898989', 'HandleVisibility', 'off');
    hold on;
    for j = 1:3
        hh(j) = plot(k,real(oss{j}(i,:)),ln{j},'linewidth',2, 'DisplayName', dn{j});
    end
    hold off;
    xlabel('$k$','fontsize',30);
    ylabel('$\omega _r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    xlim([0 4]);
    ylim([-1.5 1.5]);
    xticks(0:1:4);
    yticks(-1.5:0.5:1.5);
    legend('location','east');

    sgtitle(dm{i},'fontsize',30);
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

%% Plot real omega
figure; 
yline(0,'linewidth',1.5,'color','#898989');
hold on;
for i = 1:3
    for j = 1:3
        plot(k,real(oss{i}(j,:)),ln{i},'linewidth',1.5,'color',col(j,:));
    end
end
hold off; box on;
xlabel('$k$','fontsize',30);
xticks(0:1:4);
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('$\omega_r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

% %% Plot real c
% figure; 
% yline(0,'linewidth',1.5,'color','#898989');
% hold on;
% for i = 1:3
%     for j = 1:3
%         plot(k,real(css{i}(j,:)),ln{i},'linewidth',1.5,'color',col(j,:));
%     end
% end
% hold off;
% xlabel('$k$','fontsize',30);
% xticks(0:1:4);
% set(gca,'XMinorTick','on','YMinorTick','on')
% ylabel('$c_r$','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

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