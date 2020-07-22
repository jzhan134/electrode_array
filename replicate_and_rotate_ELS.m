clear
clc
x = -100:0.25:100;
[xx,yy] = meshgrid(x,x);
kT = 1.38e-23 * 298;
epsilon = 80 * 8.85e-12;
a = 1.5e-6;
fcm = 0.4667;
pref = 2*pi*epsilon*a^3*fcm/kT;

%% read electric field file
data = load("./library/fields/5v3_112.txt");
Ex = reshape(data(:,3),[801,801]);
Ey = reshape(data(:,4),[801,801]);
Fx = reshape(data(:,5),[801,801]);
Fy = reshape(data(:,6),[801,801]);
Ex = Ex(281:521,281:521);
Ey = Ey(281:521,281:521);
Fx = Fx(281:521,281:521);
Fy = Fy(281:521,281:521);
Ex1 = imgaussfilt(Ex,1);
Ey1 = imgaussfilt(Ey,1);
Fx1 = imgaussfilt(Fx,1);
Fy1 = imgaussfilt(Fy,1);

data = load("./library/fields/5v3_67.txt");
Ex = reshape(data(:,3),[801,801]);
Ey = reshape(data(:,4),[801,801]);
Fx = reshape(data(:,5),[801,801]);
Fy = reshape(data(:,6),[801,801]);
Ex = Ex(281:521,281:521);
Ey = Ey(281:521,281:521);
Fx = Fx(281:521,281:521);
Fy = Fy(281:521,281:521);
Ex2 = imgaussfilt(Ex,1);
Ey2 = imgaussfilt(Ey,1);
Fx2 = imgaussfilt(Fx,1);
Fy2 = imgaussfilt(Fy,1);


GB_Ex = zeros(801,801);
GB_Ey = zeros(801,801);
GB_Fx = zeros(801,801);
GB_Fy = zeros(801,801);
grid = [41,281,521];
for x = 1:3
    for y = 1:3
        i = grid(x);
        j = grid(y);
        
        if mod(x+y,2) == 1
            GB_Ex(i:i+240,j:j+240) = Ex1;
            GB_Fx(i:i+240,j:j+240) = Fx1;
            GB_Ey(i:i+240,j:j+240) = Ey1;
            GB_Fy(i:i+240,j:j+240) = Fy1;
        else
            GB_Ex(i:i+240,j:j+240) = Ex2;
            GB_Fx(i:i+240,j:j+240) = Fx2;
            GB_Ey(i:i+240,j:j+240) = Ey2;
            GB_Fy(i:i+240,j:j+240) = Fy2;
        end
    end
end
E = sqrt(GB_Ex.^2 + GB_Ey.^2);
E = E./(1e4);
figure(12)
clf
set(gcf,'position',[300 -1000 600 600])
sb = subplot(1,1,1);
set(sb,'position',[0 0 1 1])
colormap((jet))
hold on 
contourf(xx,yy, E',50,'linecolor','none');
pbaspect([1 1 1])
axis([-90 90 -90 90])
xticks('')
yticks('')
box on
% colorbar

newData = data;
newData(:,3) = reshape(GB_Ex,[801*801,1]);
newData(:,4) = reshape(GB_Ey,[801*801,1]);
newData(:,5) = reshape(GB_Fx,[801*801,1]);
newData(:,6) = reshape(GB_Fy,[801*801,1]);
save("./library/fields/array/5v3_II.txt", "newData",'-ascii')