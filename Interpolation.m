%% interpolation

x = 1:1:10
y = x.^2

lagrange(2.6,x,y)
% Lagrange 
function y=lagrange(x,pointx,pointy)
n=size(pointx,2);
L=ones(n,size(x,2));
if (size(pointx,2)~=size(pointy,2))
   fprintf(1,'\nERROR!\nPOINTX and POINTY must have the same number of elements\n');
   y=NaN;
else
   for i=1:n
      for j=1:n
         if (i~=j)
            L(i,:)=L(i,:).*(x-pointx(j))/(pointx(i)-pointx(j));
         end
      end
   end
   y=0;
   for i=1:n
      y=y+pointy(i)*L(i,:);
   end
end
end

function y_inter = f_newton_poly(x,y,x_inter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_newton_poly: Newton interpolating polynomial

% yint = f_newton_poly(x,y,x_inter)

% Uses an (n - 1)-order Newton interpolating polynomial based on n data points (x, y)
% to determine a value of the dependent variable (yint) at a given value of the independent 
% variable, x_inter.

% input:
%x       = independent variable
%y       = dependent variable
%x_inter = value of independent variable at which interpolation is calculated

% output:
%y_inter = interpolated value of dependent variable

% compute the finite divided differences in the form of a difference table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify that we have n dependent variables for n independent variables
n = length(x);
if length(y)~=n
    error('x and y must be same length'); 
end

% Making a matrix that contains the difference table
b = zeros(n,n);

% Assign dependent variables to the first column of b.
b(:,1) = y(:); % the (:) ensures that y is a column vector.

% Filling different orders of the derivatives
for jj = 2:n
    for ii = 1:n-jj+1
        b(ii,jj) = (b(ii+1,jj-1)-b(ii,jj-1))/(x(ii+jj-1)-x(ii));
    end
end

% Use the finite divided differences to interpolate

xt = 1;
y_inter = b(1,1);

for jj = 1:n-1
xt      = xt*(x_inter-x(jj));
y_inter = y_inter+b(1,jj+1)*xt;
end
end

function [a,b,c]=Quadratic(x,y)
p1=x(2)-x(3);
p2=x(3)-x(1);
p3=x(1)-x(2);
p4=x(3)^2-x(2)^2;
p5=x(1)^2-x(3)^2;
p6=x(2)^2-x(1)^2;
p7=x(2)^2*x(3)-x(2)*x(3)^2;
p8=x(1)*x(3)^2-x(1)^2*x(3);
p9=x(1)^2*x(2)-x(1)*x(2)^2;
delta=x(1)^2*(x(2)-x(3))-x(1)*(x(2)^2-x(3)^2)+1*(x(2)^2*x(3)-x(2)*x(3)^2);
a=(1/delta)*((x(2)-x(3))*y(1)+(x(3)-x(1))*y(2)+(x(1)-x(2))*y(3));
b=(1/delta)*(p4*y(1)+p5*y(2)+p6*y(3));
c=(1/delta)*(p7*y(1)+p8*y(2)+p9*y(3));
end