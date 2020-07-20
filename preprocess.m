% from particles in frame CF from cluster indexed N0 to N1, with rotation
% as the count of 90 counter-clockwise rotaton.
function out = preprocess(CF, N0, N1, rot)
x0 = N0(1);
y0 = N0(2);
x1 = N1(1);
y1 = N1(2);
a = 0.8;
CF = CF*a;
idx = ...
    round(CF(:,2)/60) == x0 & ...
    round(CF(:,3)/60) == y0;
temp = CF(idx,2:3);
meanX = mean(temp);
temp = temp - mean(temp);
while rot ~= 0
    temp = [-temp(:,2),temp(:,1)];
    rot = rot - 1;
end
temp(:,1) = temp(:,1) + x1*60;
temp(:,2) = temp(:,2) + y1*60;
out = [(0:size(temp,1)-1)',temp/a];
end