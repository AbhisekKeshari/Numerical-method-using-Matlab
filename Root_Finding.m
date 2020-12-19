%% Input parameters Abhisek-Keshari : 2018UME0126

% x_low             ===>    Lower limit At the iteration i
% x_high            ===>    Higher limit At the iteration i
% N_iter            ===>    Number of iterations
% x_actual          ===>    Actual Root
% f_x               ===>    Function defination
% tolerance         ===>    Stopping Criteria
% f_derivative      ===>    Function Derivative
x_high = 1;
f_x =@(x) x^2 - 4/9;
f_partial = @(x) x^2 - 4/9 + x;
tolerance = 1e-5;
N_iter = 100;
x_actual = 0.66;
% x_guess           ===>    Guessed Root


x_low = 0;
f_derivative = @(x) 2*x;
x_guess = 0.8;

bisection_method(x_low,x_high,N_iter,x_actual,f_x,tolerance)
False_position(x_low,x_high,N_iter,x_actual,f_x,tolerance)
Newtons(x_guess,f_x,f_derivative,N_iter,x_actual,tolerance)
Fixes_point(x_guess,f_partial,N_iter,x_actual,tolerance)
Secant(x_guess,x_low,x_high,f_x,N_iter,x_actual,tolerance)
modified_secant(x_guess,f_x,N_iter,x_actual,tolerance)

%% Functions
% Bisection Method
function root = bisection_method(x_low,x_high,N_iter,x_actual,f_x,tolerance)
    for i = 1:N_iter
        root = (x_low+x_high)/2;
        
        if(f_x(x_low)*f_x(root)<0)
           x_high = root; 
        end
        
        if(f_x(x_high)*f_x(root)<0)
           x_low = root; 
        end
        
        if((abs(root - x_actual)/x_actual)<tolerance)
           break 
        end
    end
end

% False Position Method
function root = False_position(x_low,x_high,N_iter,x_actual,f_x,tolerance)
    
    for i = 1:N_iter
       root = x_high -(f_x(x_high)*(x_low - x_high)/(f_x(x_low) - f_x(x_high)));
       
       if(f_x(root)*f_x(x_high)< 0 )
          x_low = root; 
       end
       
       if((abs(root - x_actual)/x_actual)<tolerance)
           break 
       end
    end
end

% Newtons Method 
function root = Newtons(x_guess,f_x,f_derivative,N_iter,x_actual,tolerance)
    for i = 1:N_iter
       root = x_guess;
       x_new = root - f_x(root)/f_derivative(root);
       
       if((abs(root - x_actual)/x_actual) > tolerance)
           x_guess = x_new;
       else
           break
       end
    end
end

% Fixed Point Method
function x_root = Fixes_point(x_guess,f_x,N_iter,x_actual,tolerance)
    for i = 1:N_iter
       x_root = x_guess;
       x_new = f_x(x_root);
       
       if((abs(x_new - x_actual)/x_actual)>tolerance)
           x_guess = x_new;
       else
           x_root = x_new;
           break
       end
    end
end

% Secant Method 
function root = Secant(x_guess,x_low,x_high,f_x,N_iter,x_actual,tolerance)
    
    %x_1 = x_min;
    %x_2 = x_max;
    
    for i = 1:N_iter
       
        x_new = x_high -(f_x(x_high)*(x_low - x_high)/(f_x(x_low) - f_x(x_high)));
        
        if((abs(x_new - x_guess)/x_actual)>tolerance)
            x_guess = x_new;
            x_low = x_high;
            x_high = x_new;
        else
            root = x_new;
            break
        end        
    end
end

%% - - - Non linear equation - - -
%% "modified secant"      - - -moderate speed- - -

function root = modified_secant(x_guess,f_x,N_iter,x_actual,tolerance)
    delta = 0.01;
    x_low = x_guess;

    for i = 1:N_iter
        x_new = x_low - delta*x_low*f_x(x_low)/(f_x(x_low + delta*x_low) - f_x(x_low));
        
        if((abs(x_new - x_guess)/x_actual)>tolerance)
            x_guess = x_new;
            x_low = x_new;
        else
            root = x_new;
            break
        end 
    end
end
