function [ Q , T ] = final_p4solve( U, V )
%Compute the Schur factorization of an involution with orthonormal subspace
%bases U (for eigenvalue 1) and V (for eigenvalue -1)

[num_rows_pos , num_columns_pos]= size(U); 
[num_rows_neg, num_columns_neg]=size(V); 

n = num_rows_pos ; 

total_num_columns = num_columns_pos + num_columns_neg ; 
%Should partition space
assert(n==total_num_columns); 
%I = eye(num_rows_pos); 
%lhs = ; 
%Unstable
T = V - U*(U\V) ; 

[W , R] = qr(T , 0);

X_int_1 = -2 * ( U*( U \ V));

X_int_2= X_int_1/R;

X = U\X_int_2; 

Y = [U W]; 

T = [ eye(num_columns_pos) X ; zeros(num_columns_neg, num_columns_pos) , -1*eye(num_columns_neg)];

end

