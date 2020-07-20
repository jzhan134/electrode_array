clear
clc
set(0,'DefaultFigureWindowStyle','normal')
addpath('./Matlab');
addpath('./Matlab/Image_Analysis_v0626');
figure(10)
set(gcf,'position',[10 10 900 900])
movie = false;
% data = load('library/configurations/3600_quench2.txt');
% pnum = size(data,1);
data = load('traj/xyz0.dat');
% pnum = size(data,1);
pnum = max(data(:,1))+1;
a = 0.8;
frame = size(data,1)/pnum;
if movie
    v = VideoWriter('traj_demo.avi');
    v.FrameRate = 50;
    open(v);
end
% C6 = [];
for f = 1:5:frame
    CF = data((f-1)*pnum+1:f*pnum,:);
%     [~, ~, ~, ~, G_C6, ~, ~,~] = IMAGE_ANALYSIS(CF);
%     G_C6
%     C6 = cat(1,C6, G_C6);
    CF_plot = CF*a;
    clf
    hold on
    for dx = -90:60:90
        plot([dx dx],[-90 90],'k-')
        plot([-90 90],[dx dx],'k-')
    end
%     plot_electrodes
    plot(CF_plot(:,2),CF_plot(:,3),'ko',...
        'markersize', 6,...
        'MarkerFaceColor','k');
    display_cluster_sizes;
    xticks('')
    yticks('')
    axis([-90 90 -90 90])
    pbaspect([1,1,1])
    sb = subplot(1, 1, 1);
    set(sb,'position',[0 0 1 1])
    box on
    hold off
    annotation('textbox', ...
        [.9 .0 .1 .1], ...
    'String', num2str(f/10), ...
    'EdgeColor', 'none', ...
    'BackgroundColor','w',...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom',...
    'FontSize', 16);
    drawnow;
%     if movie
%         this_frame = getframe(gcf);
%         writeVideo(v,this_frame.cdata);
%     end
end
if movie
    close(v)
end
% figure(1)
% plot(C6)
% axis([0 inf 0 1])