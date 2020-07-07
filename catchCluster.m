clear
clc
set(0,'DefaultFigureWindowStyle','normal')
figure(10)
clf
set(gcf,'position',[10 10 900 900])
movie = false;
CF = load('library/configurations/3600_quench.txt');
a = 0.8;
idx = round(CF(:,2)*a/60)+ round(-CF(:,3)*a/60)*3;

row = [0,0];
col = [0,1];
curr = [];
for i = 1:length(row)
    r = row(i);
    c = col(i);
    idx0 = r+c*3;
    curr = cat(1,curr,CF(idx == idx0,[1:3]));
end
plot(curr(:,2)*a,curr(:,3)*a,'ko',...
    'markersize', 4,...
    'MarkerFaceColor','k');
axis([-90 90 -90 90])
pbaspect([1,1,1])
save('library/configurations/test.txt','curr','-ascii')
