function offDiagonal = minGry(X,t)
R = [cos(t*pi/180), -sin(t*pi/180); sin(t*pi/180), cos(t*pi/180)];
gyr = (X*R)'*(X*R);
offDiagonal = abs(gyr(1,2));
end