clear
clc
file ="./library/fields/array/3x3_frame.txt";
set(0,'DefaultFigureWindowStyle','normal')
x = -100:0.25:100;
[xx,yy] = meshgrid(x,x);
kT = 1.38e-23 * 298;
epsilon = 80 * 8.85e-12;
a = 1.5e-6;
fcm = 0.4667;
pref = 2*pi*epsilon*a^3*fcm/kT;

%% read electric field file
data = load(file);
Ex = reshape(data(:,3),[length(x),length(x)]);
Ey = reshape(data(:,4),[length(x),length(x)]);
Fx = reshape(data(:,5),[length(x),length(x)]);
Fy = reshape(data(:,6),[length(x),length(x)]);
% Ex = imgaussfilt(Ex,10);
% Ey = imgaussfilt(Ey,10);
E = sqrt(Ex.^2 + Ey.^2);
E = E./(1e4);
for i = 1:801
    for j = 1:801
        if (E(i,j) > 4)
            E(i,j) = 4;
        end
    end
end
            
figure(12)
clf
set(gcf,'position',[80 50 900 900])
sb = subplot(1,1,1);
set(sb,'position',[0 0 1 1])
colormap((jet))
hold on 
contourf(xx,yy, E',50,'linecolor','none');
pbaspect([1 1 1])
axis([-90 90 -90 90])
xticks('')
% axis([-30 30 -30 30])
yticks('')
% caxis([0 5])
box on
% [qx, qy] = meshgrid(-100:10:100);
% qEx = Ex(1:40:end,1:40:end);
% qEy = Ey(1:40:end,1:40:end);
% qExn = qEx./sqrt(qEx.^2+qEy.^2);
% qEyn = qEy./sqrt(qEx.^2+qEy.^2);
% hold on
% h = quiver(qy,qx,qExn,qEyn,'color','w','LineStyle','none','AutoScale','on',...
%     'AutoScaleFactor',1,'MaxHeadSize',0.1,'linewidth',2);
% h.Head.LineStyle = 'solid';  
