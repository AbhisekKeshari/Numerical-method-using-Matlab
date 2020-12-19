%% Input parameters
Coeff_matrix = [3,    -0.1,    -0.2, 0.1; 
                0.1,     7,    -0.3, 2.5; 
                0.3,  -0.2,      10, 3;
                0.12,   1,  2, 4 ];

Const_matrix = [7.85, -19.3,   71.4 15];

Size = size(Const_matrix,2);

tollerance = 1e-5;

gauss_jordan(Coeff_matrix ,Const_matrix)
gauss_elemination(Coeff_matrix ,Const_matrix)
gauss_seidel(Coeff_matrix , Const_matrix, tollerance, Size)
Jacobi_method(Coeff_matrix , Const_matrix, tollerance, Size)

%% Function Block

%gauss_elemination method
function x_vector = gauss_elemination(A ,b)
    N = size(A,1);
    
    for i = 1:N-1
        
       for j = i+1:N
           factor = A(j,i)/A(i,i);
           %A(j,:) = A(j,:) - (A(j,i)/A(i,i)).*A(i,:);
           for k = i:N
               A(j,k) = A(j,k) - (factor).*A(i,k);
           end
           b(j) = b(j) - (factor)*b(i);
       end
    end
    
    x_vector = zeros(1,N);
    x_vector(N) = b(N)/A(N,N);

    for i = N-1:-1:1
        sum = b(i);
        
        for j = i+1:N
           sum = sum -A(i,j)*x_vector(j);
        end
        
        x_vector(i) = sum/A(i,i);
    end

end

%gauss_jordan method
function x_vector = gauss_jordan(A ,b)
    N = size(A,1);
    
    for i = 1:N-1
       for j = i+1:N
           factor = A(j,i)/A(i,i);
           %A(j,:) = A(j,:) - (A(j,i)/A(i,i)).*A(i,:);
           for k = i:N
               A(j,k) = A(j,k) - (factor).*A(i,k);
           end
           b(j) = b(j) - (factor)*b(i);
       end
    end

    b(N) = b(N)/A(N,N);
    A(N,N) = 1;
    
    for i = N:-1:2
        for j = i-1:-1:1
           factor = A(j,i)/A(i,i);
           %A(j,:) = A(j,:) - (A(j,i)/A(i,i)).*A(i,:);
           for k = N:-1:i
               A(j,k) = A(j,k) - (factor).*A(i,k);
           end
           b(j) = b(j) - (factor)*b(i);
        end
        b(i-1) =   b(i-1)/A(i-1,i-1);
        A(i-1,i-1) = 1;
    end

    x_vector = b;
end

%gauss_seidel method
function x_vector = gauss_seidel(A , b, tolerance, SIZE)
    norm_val=Inf;
    
    x_vector = zeros(1,SIZE);
    
    while norm_val>tolerance
        x_old = x_vector;
       
        for i=1:SIZE
        
            SUM=0;
        
            for j=1:i-1
                SUM=SUM+A(i,j)*x_vector(j);
            end
        
            for j=i+1:SIZE
                SUM=SUM+A(i,j)*x_old(j);
            end
        
            x_vector(i)=(1/A(i,i))*(b(i)-SUM);
         
        end
        
        norm_val=norm(x_old-x_vector);
    end
end

%Jacobi_method
function x_vector = Jacobi_method(A , b, tolerance, SIZE)
    norm_val=Inf;
    
    x_vector = zeros(1,SIZE);
    
    while norm_val>tolerance
        x_old = x_vector;
       
        for i=1:SIZE
        
            SUM=0;
        
            for j=1:SIZE
                if j~=i
                    SUM=SUM+A(i,j)*x_vector(j);
                end
            end
        
            x_vector(i)=(1/A(i,i))*(b(i)-SUM);
         
        end
        
        norm_val=norm(x_old-x_vector);
    end
end

%% Relaxation - - - 30 sept - - -