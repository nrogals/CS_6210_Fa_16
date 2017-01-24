function [ x, fmin ] = final_p3solve( A, M )
%Return the minimizing vector x and the minimum value fmin

n=max(size(A)); 

C_column_one = ones ( n , 1 );
C = zeros ( n , n ); 
C(: , 1 ) = C_column_one; 

C_trans = transpose(C); 

display(C_trans); 

L = chol(M,'lower') ; 

L_trans = transpose(L); 

%D_trans is constraint matrix
D_trans = C_trans / L_trans ; 

%D_other_way = C_trans * inv(L_trans); 
%display(D_other_way); 

D= transpose(D_trans); 

[Q_full,R_full] = qr(D) ; 

Q_2 = Q_full(:, 2 : end); 
%Since it is fixed that C has only one constraint, W is basis for solution
%they say. y_p is equal to Q_1 * inv (transpose(R_1)) d where d = 0 and so
%by linearity yp=0; 




%Form new matrix over lower dim space for optimization 
%transpose(Q_2)*inv(L)*A*inv(L_trans) * Q_2

B=A/L_trans ; 
middle_matrix = L\B ; 
%display(middle_matrix); 


%middle_matrix_diff = inv(L) * A *inv(L_trans); 
%Numerically unstable
%{
display(middle_matrix); 
display(middle_matrix_diff); 
assert(1<0);
 %}

%display(size(transpose(Q_2))); 
%display(size(middle_matrix)); 
%display(size(Q_2)); 


opt_matrix = transpose(Q_2) *middle_matrix*Q_2 ; 

%assert(1<0);
%if matrix is always upper triangular read off eigenvalues from the
%diagonal...save time. 

eval_fun = @(x)  (transpose(x) * opt_matrix * x ) / (transpose(x) * x) ; 

minimizing_z = find_minimizing_eigvec ( opt_matrix, eval_fun) ; 
%Immediately equivalent to minimizing y since yp = ; 



minimizing_y = Q_2 * minimizing_z ; 

minimizing_x = L_trans \ minimizing_y ; 

display(minimizing_x); 



f_min = (transpose(minimizing_x) * A * minimizing_x ) ./ (transpose(minimizing_x) *M * minimizing_x ) ; 
fmin=f_min ; 
x = minimizing_x ; 

end



function minimizing_z = find_minimizing_eigvec (opt, eval_fun)

%Opt may be upper triangular, always ???

%Eigenvalues, eigenvecs are stationary points. 


[V,D] = eig(opt) ; 

[num_rows, num_columns] = size(V); 

current_num = inf ; 
minimizing_vec = NaN; 

for i = 1 : num_columns
    trial_vec = V(:, i) ; 
    num = eval_fun(trial_vec) ; 
    
    if num < current_num 
        current_num = num; 
        minimizing_vec = trial_vec ; 
    end
end 

minimizing_z = minimizing_vec ; 

end




