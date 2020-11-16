close all; clear all;% clc
load('diffk.mat');
k = k(k<=4);
o = o(k<=4);
o1 = imag(o);
o1r = real(o);
o1(isnan(o1)) = nan;
o1r(imag(o1)<0) = nan;
inflec_pt = -0.74708299;
for i = 1:length(k)
    [maxval, maxpos] = max(abs(phi{i}(:,1)));
    if (real(phi{i}(maxpos,1)) < 0)
        phi{i} = -phi{i}/maxval;
    else
        phi{i} = phi{i}/maxval;
    end
    ang(:,i) = angle(phi{i}(:,1));
end

%% Plot phi
figure('position',[0 0 1920 1080]);
for i = 1:length(k)
    yt = sort([-6:1:0 z_c(i)]);
    ind = find(yt==z_c(i));
    
    subplot('Position',[0.8 0.7 0.18 0.25]);
    plot(k,o1,'linewidth',2);
    ylim([0 0.04]);
    xlim([0 4]);
    hold on; xline(k(i), '--r', 'linewidth', 2); hold off;
    set(gca,'fontsize',20);
    xlabel('$\tilde{k}$', 'Interpreter', 'LaTeX','fontsize',30);
    ylabel('$\tilde{\omega_i}$', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    ax = gca;
    ax.YAxis.Exponent = -2;
    grid on;
    
    subplot(6,5,(6:5:30));
    plot(abs(phi{i}(:,1)),z(:,i),'-o','linewidth',1.5,'markersize',3);
    xlim([0 1]);
    ylim([-6 0]);
    if (~isnan(z_c(i)))
        hold on; yline(z_c(i), '-.r', 'linewidth', 1.5); 
        yline(inflec_pt, '--k', 'linewidth', 1.5); 
        hold off;
        yticks(yt);
        ax = gca;
        ax.YTickLabel{ind} = ['\color{red}' ax.YTickLabel{ind}];
    end
    xlabel('$||\phi||$','FontSize',30, 'Interpreter', 'LaTeX');
    ylabel('z','FontSize',30, 'Interpreter', 'LaTeX','rotation',0);
    grid on;
    set(gca,'fontsize',24);
    
    subplot(6,5,(7:5:30));
    plot(ang(:,i),z(:,i),'-o','linewidth',1.5,'markersize',3);
    if (~isnan(z_c(i)))
        hold on; yline(z_c(i), '-.r', 'linewidth', 1.5); 
        yline(inflec_pt, '--k', 'linewidth', 1.5); 
        hold off;
        yticks(-6:1:0);
    end
    xlabel('$angle(\phi)$','FontSize',30, 'Interpreter', 'LaTeX');
    xlim([-pi pi]);
    xticks(-pi:pi:pi);
%     ylabel('z','FontSize',30, 'Interpreter', 'LaTeX','rotation',0);
    ylim([-6 0]);
    set(gca,'fontsize',24,'xticklabels',{'-\pi','0','\pi'});
    grid on;
    
    subplot(6,5,(8:5:30));
    plot(abs(phi{i}(:,2)),z(:,i),'-o','linewidth',1.5,'markersize',3);
    if (~isnan(z_c(i)))
        hold on; yline(z_c(i), '-.r', 'linewidth', 1.5); 
        yline(inflec_pt, '--k', 'linewidth', 1.5); 
        hold off;
        yticks(-6:1:0);
    end
    xlabel('$||\frac{d \phi }{dz}||$','FontSize',30, 'Interpreter', 'LaTeX');
    if (k(i) < 0.7)
        xlim([0 2]);
    else
        xlim([0 4]);
    end
%     ylabel('z','FontSize',30, 'Interpreter', 'LaTeX', 'rotation',0);
    ylim([-6 0]);
    set(gca,'fontsize',24);
    grid on;
    
    subplot(6,5,(9:5:30));
    plot(abs(phi{i}(:,3)),z(:,i),'-o','linewidth',1.5,'markersize',3);
    if (~isnan(z_c(i)))
        hold on; yline(z_c(i), '-.r', 'linewidth', 1.5); 
        yline(inflec_pt, '--k', 'linewidth', 1.5); 
        hold off;
        yticks(-6:1:0);
    end
    xlabel('$||\frac{d^2 \phi }{dz^2}||$','FontSize',30, 'Interpreter', 'LaTeX');
    if (k(i) < 0.7)
        xlim([0 2.1]);
    else
        xlim([0 10]);
    end
%     ylabel('z','FontSize',30, 'Interpreter', 'LaTeX', 'rotation',0);
    ylim([-6 0]);
    set(gca,'fontsize',24);
    grid on;
    
    tet = sprintf('k = %.2f',k(i));
    sgtitle(tet,'FontSize',30);
    F(i) = getframe(gcf);
end
writerObj = VideoWriter('phi.mp4','MPEG-4');
writerObj.FrameRate = 10;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);