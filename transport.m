clear
clc
set(0,'DefaultFigureWindowStyle','normal')
addpath('./Matlab');
addpath('./Matlab/Image_Analysis_v0626');
figure(10)
set(gcf,'position',[10 10 900 900])
movie = false;
% data = load('library/configurations/single_cluster.txt');
% pnum = size(data,1);
data = load('traj/xyz0.dat');
pnum = max(data(:,1))+1;
a = 0.8;
frame = size(data,1)/pnum;
if movie
    v = VideoWriter('traj_demo.avi');
    v.FrameRate = 2;
    open(v);
end
for f = 1:10:frame
    CF = data((f-1)*pnum+1:f*pnum,:);
    CF_plot = CF*a;
    clf
    hold on
    plot_electrodes(0);
    plot(CF_plot(:,2),CF_plot(:,3),'ko',...
        'markersize', 6,...
        'MarkerFaceColor','k');
%     display_cluster_sizes;
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
    if movie
        this_frame = getframe(gcf);
        writeVideo(v,this_frame.cdata);
    end
end
if movie
    close(v)
end
% figure(1)
% plot(res)
% axis([0 inf 0 1])