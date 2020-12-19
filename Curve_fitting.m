% Curve Fitting
x = 1:0.8:10;
y = x.^2;

xx1 = [1 2 3 4 5];
yy1 = [0.5 1.7 3.5 5.7 8.4];

%xx = [3 4 5 7 8 9 11 12];
%yy = [1.6 3.6 4.4 3.4 2.2 2.8 3.8 4.6];

xx = [0 1 2 3 4 5];
yy = [2.1 7.7 13.6 27.2 40.9 61.1];

least_squarefit(xx,yy)
log_leastfit(xx,yy)
cubic_leastfit(xx,yy)
quad_leastfit(xx,yy)

function least_squarefit(x,y)

f1 = figure(1);
plot(x, y, 'ko', 'Markersize', 10,...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');hold on;
xlabel('x');
ylabel('y');
set(gca, 'FontSize', 20);

sum_x = sum(x);
sum_y = sum(y);
sum_xy = sum(x.*y);
sum_x2 = sum(x.^2);
mean_x = mean(x);
mean_y = mean(y);

N = length(x);

a1 = (N*sum_xy -sum_x*sum_y)/(N*(sum_x2) - (sum_x).^2);
disp('Slope of Straight Line = '); disp(a1)
a0 = mean_y - a1*mean_x;
disp('Intercept with Y-axis = '); disp(a0)
y_new = a0 + a1*x;

plot(x, y_new, 'r-', 'linewidth', 2);
set(gca,'FontSize',20);
hold off;
end

function log_leastfit(xx1,yy1)
f1 = figure(1);
plot(xx1,yy1,'b-','LineWidth',2.0,'Marker','o','MarkerSize',10); 
xlabel('xx');
ylabel('yy');
set(gca, 'FontSize', 20);
zz1 = log(yy1); 

% y = a0 + a1*x % slope = a1% intercept with Y-axis = a0;
sum_xx1 = sum(xx1);
sum_zz1 = sum(zz1);
sum_xxzz = sum(xx1.*zz1);
sum_xx2 = sum(xx1.^2);
mean_xx1 = mean(xx1);
mean_zz1 = mean(zz1);

NN = length(xx1);

a11 = (NN*sum_xxzz -sum_xx1*sum_zz1)/(NN*(sum_xx2) - (sum_xx1).^2);
disp('Slope of Straight Line = '); disp(a11)
a00 = mean_zz1 - a11*mean_xx1;
disp('Intercept with Y-axis = '); disp(a00)
y_new1 = a00 + a11*xx1;

f2 = figure(2);
plot(xx1,zz1,'o','MarkerSize',10);  hold on;
plot(xx1, y_new1, 'r-', 'linewidth', 2); 
xlabel('xx'); ylabel('ln(yy)');
set(gca,'FontSize',20);
end

function cubic_leastfit(xx,yy)
f1 = figure(1);
plot(xx,yy,'ko','MarkerSize',12);hold on;
xlabel('X');ylabel('Y');

set(gca,'FontSize',20);

n = length(xx);

sum_x = sum(xx);
sum_x2 = sum(xx.^2);
sum_x3 = sum(xx.^3);
sum_x4 = sum(xx.^4);
sum_x5 = sum(xx.^5);
sum_x6 = sum(xx.^6);

sum_y = sum(yy);
sum_xy = sum(xx.*yy);
sum_x2y = sum(xx.^2.*yy);
sum_x3y = sum(xx.^3.*yy);

A1 = [n sum_x sum_x2 sum_x3; sum_x sum_x2 sum_x3 sum_x4; sum_x2 sum_x3 sum_x4 sum_x5; sum_x3 sum_x4 sum_x5 sum_x6]
b1 = [sum_y sum_xy sum_x2y sum_x3y]'

aa = A1\b1;

ynew = aa(1) + aa(2)*xx + aa(3)*xx.^2 + aa(4)*xx.^3;
plot(xx,ynew,'b-','LineWidth',2);
l1 = legend('Observations','Quadratic fit');
l1.Box = 'off';
l1.Location = 'best';
end

function quad_leastfit(xx,yy)
f1 = figure(1);
plot(xx,yy,'ko','MarkerSize',12); hold on;
xlabel('X'); ylabel('Y');
set(gca,'FontSize',20);
N = length(xx);

sum_xx   = sum(xx);
sum_xx2  = sum(xx.^2); 
sum_xx3  = sum(xx.^3); 
sum_xx4  = sum(xx.^4); 
sum_yy   = sum(yy);
sum_xx_yy = sum(xx.*yy);
sum_xx2_yy= sum(xx.^2.*yy);

A1 = [N sum_xx sum_xx2; sum_xx sum_xx2 sum_xx3; sum_xx2 sum_xx3 sum_xx4]
b1 = [sum_yy sum_xx_yy sum_xx2_yy]'
aa = A1\b1
ynew = aa(1) + aa(2)*xx + aa(3)*xx.^2;
plot(xx,ynew,'b-','LineWidth',2);
l1 = legend('Observations','Quadratic fit');
l1.Box = 'off';
l1.Location = 'best';
end