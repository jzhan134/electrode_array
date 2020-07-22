% [col,row], left: col = -1; 
% bot: row = -1

clear
clc
set(0,'DefaultFigureWindowStyle','normal')
a = 0.8;
data = load('library/configurations/test2.txt');
CF = data;
% N0 = [0, 0];
% N1 = [0,-1];
N0 = [0, 0];
N1 = [0, 0];
out = preprocess(CF, N0, N1, 3);
res = out;
% N0 = [1,1];
% N1 = [0,-1];
% out2 = preprocess(CF, N0, N1, 3);
% res = [out;out2];
res(:,1) = (0:size(res,1)-1)';

figure(11)
clf
set(gcf,'position',[10 10 900 900])
hold on
for dx = -90:60:90
    plot([dx dx],[-90 90],'k-',[-90 90],[dx dx],'k-')
end
plot(res(:,2)*a,res(:,3)*a,'ko','markersize', 6,'MarkerFaceColor','k');
xticks('')
yticks('')
axis([-90 90 -90 90])
pbaspect([1,1,1])
sb = subplot(1, 1, 1);
set(sb,'position',[0 0 1 1])
box on
hold off
drawnow;
save('library/configurations/single_cluster.txt','res','-ascii')
