close all; clear all;% clc
load('modeshape.mat');
for i = 1:length(k)
    % Data process
    [maxval, maxpos] = max(abs(mk(i).phi));
    if (real(mk(i).phi(maxpos)) < 0)
        mk(i).phi = -mk(i).phi/maxval;
        mk(i).up = -mk(i).up/maxval;
        mk(i).wp = -mk(i).wp/maxval;
    else
        mk(i).phi = mk(i).phi/maxval;
        mk(i).up = mk(i).up/maxval;
        mk(i).wp = mk(i).wp/maxval;
    end
end
%% Plot phi
figure('position',[0 0 1440 1280]);
for i = 1:length(k)

    subplot(1,2,1);
    plot(real(mk(i).phi),z,'-k','linewidth',1);
    % hold on; plot(real(phi(minpos)),z(minpos),'ro'); hold off;
    xlabel('real(\phi)','FontSize',30);
    ylabel('z','FontSize',30);
    if (k(i) < 0.7)
        xlim([-0.5 1]);
    else
        xlim([0 1]);
    end
    ylim([-h(k(i)) 0]);
    title('real','FontSize',30);
    grid on;
    set(gca,'fontsize',24);
    subplot(1,2,2);
    plot(imag(mk(i).phi),z,'-k','linewidth',1);
    % hold on; plot(imag(phi(minpos)),z(minpos),'ro'); hold off;
    xlabel('imag(\phi)','FontSize',30);
    if (k(i) < 0.7)
        xlim([-0.2 0.15]);
    else
        xlim([0 0.11]);
        xticks(0:0.025:0.1);
    end
    ylabel('z','FontSize',30);
    ylim([-h(k(i)) 0]);
    title('imag','FontSize',30);
    set(gca,'fontsize',24);
    grid on;
    tet = sprintf('\\phi, k = %.2f',k(i));
    sgtitle(tet,'FontSize',30);
    F(i) = getframe(gcf);
end
writerObj = VideoWriter('phiT.mp4','MPEG-4');
writerObj.FrameRate = 10;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);
%% Plot u
figure('position',[0 0 1440 1280]);
for i = 1:length(k)
    subplot(1,2,1);
    plot(real(mk(i).up),z,'.-k','linewidth',1);
    % hold on; plot(real(up(minpos)),z(minpos),'ro'); hold off;
    xlabel('real(u'')','FontSize',30);
    if (k(i) < 0.7)
        xlim([-5 2]);
    else
        xlim([0 1]);
    end
    ylabel('z','FontSize',30);
    ylim([-h(k(i)) 0]);
    title('real','FontSize',30);
    grid on;
    set(gca,'fontsize',24);
    subplot(1,2,2);
    plot(imag(mk(i).up),z,'.-k','linewidth',1);
    % hold on; plot(imag(up(minpos)),z(minpos),'ro'); hold off;
    xlabel('imag(u'')','FontSize',30);
    if (k(i) < 0.7)
        xlim([-1.1 1]);
    else
        xlim([-0.7 0.3]);
    end
    ylabel('z','FontSize',30);
    ylim([-h(k(i)) 0]);
    title('imag','FontSize',30);
    set(gca,'fontsize',24);
    grid on;
    tet = sprintf('u'', k = %.2f',k(i));
    sgtitle(tet,'FontSize',30);
    F(i) = getframe(gcf);
end
writerObj = VideoWriter('up.mp4','MPEG-4');
writerObj.FrameRate = 6;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);
%% Plot w
figure('position',[0 0 1440 1280]);
for i = 1:length(k)
    subplot(1,2,1);
    plot(real(mk(i).wp),z,'k','linewidth',1);
    % hold on; plot(real(wp(minpos)),z(minpos),'ro'); hold off;
    xlabel('real(w'')','FontSize',30);
    ylabel('z','FontSize',30);
    ylim([-h(k(i)) 0]);
    title('real','FontSize',30);
    set(gca,'fontsize',24);
    grid on;
    subplot(1,2,2);
    plot(imag(mk(i).wp),z,'k','linewidth',1);
    % hold on; plot(imag(wp(minpos)),z(minpos),'ro'); hold off;
    xlabel('imag(w'')','FontSize',30);
    ylabel('z','FontSize',30);
    ylim([-h(k(i)) 0]);
    title('imag','FontSize',30);
    set(gca,'fontsize',24);
    grid on;
    sgtitle('w''','FontSize',30);
    F(i) = getframe(gcf);
end
writerObj = VideoWriter('wp.mp4','MPEG-4');
writerObj.FrameRate = 6;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);