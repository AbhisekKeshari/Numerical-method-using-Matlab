%% Newton Cotter - - -
%  grid based integration
%  Richardsosn extrapollation


fun1 = @(x) 0.2 + 25*x - 200*x.^2 + 675*x.^3 -900*x.^4 + 400*x.^5; 
x_low = 0;
x_upper = 0.8;
true_value = 1.640533;
N_iter = 100;

Trapazoidal(x_low,x_upper,true_value,fun1);
Multiple_trapazoidal(x_low,x_upper,true_value,fun1,N_iter);
Simpson(x_low,x_upper,true_value,fun1);
Multiple_simp(x_low,x_upper,true_value,fun1,N_iter);
simp38(x_low,x_upper,true_value,fun1);
Multiple_simp38(x_low,x_upper,true_value,fun1,N_iter);

%%

function I_tz1 = Trapazoidal(x_low,x_upper,true_value,fun1)
a      = x_low;      b = x_upper;    % lower and upper limits of the integral 
I_true = true_value; 

I_tz1 = (b - a)*( fun1(a) + fun1(b) )/2; 
Error_true_tz1 = abs( (I_true - I_tz1)/I_true ) * 100;
fprintf('Single application trapezoidal rule true relative perc error = %f \n', Error_true_tz1); 
end

function I_m_trap = Multiple_trapazoidal(x_low,x_upper,true_value,fun1,N_iter)
N = N_iter;
a      = x_low;      b = x_upper;    % lower and upper limits of the integral 
I_true = true_value; 

h = (b - a)/N;
%x = x_lower:h:x_upper;
f_lower = fun1(a);
f_upper = fun1(b);

f_mid = 0;         % Initiation of sum of integrations for non-boundary points

for ii = 1 : N-1 
  f_mid = f_mid + fun1(a + ii*h);
end

I_m_trap = (h/2)*(f_lower + 2*f_mid + f_upper);
fprintf('Multiple application trapezoidal rule result = %f \n', I_m_trap); 
Error_true_tzn = abs( (I_true - I_m_trap)/I_true ) * 100;
fprintf('Multiple application trapezoidal rule true relative perc error = %f \n', Error_true_tzn); 
end

function I_simp_13 = Simpson(x_low,x_upper,true_value,fun1)

a      = x_low;      b = x_upper;    % lower and upper limits of the integral 
I_true = true_value; 

x_mid = a + (b-a)/2;   % The mid point
I_simp_13 = (b - a)*( fun1(a) + 4*fun1(x_mid) + fun1(b) )/6; 
fprintf('Single application Simpson 1/3 rule result = %f \n', I_simp_13); 
Error_true_s13 = abs( (I_true - I_simp_13)/I_true ) * 100;
fprintf('Multiple application trapezoidal rule true relative perc error = %f \n', Error_true_s13); 
end

function I_M_simp13 = Multiple_simp(x_low,x_upper,true_value,fun1,N_iter)
N = N_iter;
a      = x_low;      b = x_upper;    % lower and upper limits of the integral 
I_true = true_value; 

h = (b - a)/N;
simp13 = fun1(a) + fun1(b);

% Calculating odd numbered terms  
for ii = 1:2:N-1
  simp13 = simp13 + 4*fun1(a+ii*h);
end

% Calculating even numbered terms  
for ii = 2:2:N-2
  simp13 = simp13 + 2*fun1(a+ii*h);
end

I_M_simp13 = (h/3)*simp13 ;

fprintf('MULTIPLE APPLICATION SIMPSON 1/3 RULE RESULT = %f \n', I_M_simp13); 
Error_true_s13n = abs( (I_true - I_M_simp13)/I_true ) * 100;
fprintf('Multiple application trapezoidal rule true relative perc error = %f \n', Error_true_s13n);
end

function I_simp_38 = simp38(x_low,x_upper,true_value,fun1)

a      = x_low;      b = x_upper;    % lower and upper limits of the integral 
I_true = true_value; 

x_m1 = a + (b-a)/3;     % The inside point 1
x_m2 = a + 2*(b-a)/3;   % The inside point 2

I_simp_38 = (b - a)*( fun1(a) + 3*fun1(x_m1) + 3*fun1(x_m2) + fun1(b) )/8; 
fprintf('Single application Simpson 1/3 rule result = %f \n', I_simp_38);
Error_true_s38 = abs( (I_true - I_simp_38)/I_true ) * 100;
fprintf('Multiple application trapezoidal rule true relative perc error = %f \n', Error_true_s38);
end

function m_simp38 = Multiple_simp38(x_low,x_upper,true_value,fun1,N_iter)
N = N_iter;
a      = x_low;      b = x_upper;    % lower and upper limits of the integral 
I_true = true_value; 

h = (b - a)/N;
% Sample Vector 
x = a:h:b;
% Function at the endpoints
f_lower = fun1(x(1));
f_upper = fun1(x(end));
    
simp3 = 0;
simp2 = 0;

for ii = 2:3:N -1
  simp3 = simp3 + 3*(fun1(x(ii)) + fun1(x(ii+1)));
end

for ii = 4:3:N -2
  simp2 = simp2 + 2*fun1(x(ii));
end

% Bringing everything together
m_simp38 = (3/8)*h*(f_lower + simp3 + simp2 + f_upper);
fprintf('MULTIPLE APPLICATION SIMPSON 3/8 RESULT = %f\n',m_simp38);
end

