close all; clear all;% clc
N = 401;
load('neutral_pts.mat');
c_r = real(omax)./kmax;

fig1 = figure('position',[0 0 1920 1280]);
for i = 1:length(kmax)
%     fig1 = figure('position',[0 0 1920 1280]);
    [~,minpos] = min(abs(U(1:fix(N+1)/2)-c_r(i)));
    subplot(1,4,1);
    plot(U(1:fix(N+1)/2),z(1:fix(N+1)/2),'k','linewidth',1);
    hold on; 
    plot(U(minpos),z(minpos),'ro','markersize',10,'markerfacecolor','r'); 
    hold off;
    ylim([0 1]);
    set(gca,'fontsize',20);
    grid on;
    xlabel('$\tilde{U}$', 'Interpreter', 'LaTeX','fontsize',30);
    ylabel('$\tilde{z}\ $', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');

    subplot(1,4,[2 3 4]);
    plot(Re0,k0,'o','color','#737975');
    xlim([1e+3 1e+9]);
    ylim([0.1 1.1]);
    set(gca,'fontsize',20,'Xscale','log');
    xlabel('$Re$', 'Interpreter', 'LaTeX','fontsize',30);
    ylabel('$\tilde{k}\ $', 'Interpreter', 'LaTeX','fontsize',30,'rotation',0, 'HorizontalAlignment','right');
    hold on
    plot(Rmax,kmax,'k','linewidth',1);
    plot(Rmax(i),kmax(i),'ro','markersize',10,'markerfacecolor','r');
    hold off
    grid on
    F(i) = getframe(gcf);
end

writerObj = VideoWriter('zc.mp4','MPEG-4');
writerObj.FrameRate = 6;
writerObj.Quality = 100;
open(writerObj);
for i = 1:length(F)
    frame = F(i);
    writeVideo(writerObj, frame);
end
close(writerObj);