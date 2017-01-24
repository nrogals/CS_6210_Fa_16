%Test.m

%{
% Test Problem 2 

n=4;
x = rand(n , 1); 

row = rand(n , 1); 
column = rand(n, 1); 

row(1)=1; 
column(n)=1; 

A = toeplitz(column , row);

display(A);

b= A*x; 

b_calc = hw1toeplitz(column,row,x) ;

%display(b_calc); 
%display(b); 

display ********** ; 


A = hankel(column, row); 
display(A); 
display(x); 
display(A*x); 

b=A*x; 
b_calc= hankel_multiply(transpose(row) , transpose(column) , x ); 

display(b); 
display(b_calc);

%}


%{
%Final Problem 3 Main 


%Matrix M is spd
%A is symmetric


M = [2 -1 0 ; -1 2 -1 ; 0 -1 2 ] ; 

B = [1 3 5; 7 8 -3 ; 1 8 6]; 
B_trans = transpose(B); 

A = B_trans * B ; 
normA = norm(A , 'fro') ; 
A = A ./ normA ; 
%Scale nicely. 

fun = @(x) (transpose(x) * A * x) / (transpose(x) * M * x) ; 

x0 = [0.3634; -0.4266 ; 0.0632]; 

W = [1 1 1]; 
b = 0; 

x_actual = fmincon(fun,x0,[], [], W,b) ; 

val = fun(x_actual) ; 

[ x, fmin_solved ] = final_p3solve( A, M ) ; 

display(x); 
display(x_actual); 
display(fmin_solved); 
display(val); 

%}

%{

%Problem 4 Involution 


A = [ 3 , - 4 , 4 ; 
    0   , -1 , 0;
    -2 , 2, -3];

display(eig(A));


[V , J ] = jordan(A); 

Pos_Space = V(:, 1 ); 
Neg_Space = V(: , 2:end); 


[Q_pos , R_pos] = qr(Pos_Space, 0); 
[Q_neg , R_neg] = qr(Neg_Space, 0 ); 

display(Q_pos);

[U, T] = schur(A); 

[ W , T ] = final_p4solve( Q_pos, Q_neg ) ; 
%Returns W and X 


W_trans = transpose(W); 

display(W); 
Orthonormal = W_trans * W ; 
display(Orthonormal); 
display(T); 



calc_A = W * T* transpose(W); 


display(calc_A); 
display(A);

assert(1<0);


X_actual = T(1, 2:end); 
W_actual = U(: , 2:end); 

display ***************
display(X); 
display(W); 

display ***************
display(V); 
display(J);
display(T); 
display(X_actual); 
display(W_actual); 


%}

%{

%5 Main Test 

n=70; 

A = rand(n); 


b = rand(n,1); 


k=1000; 

[ x] = final_p5solve( A , b , k ) ; 

display(size(x));
display(size(A));


res = A*x-b; 


%Not in Kyrolv subspace...
x_actual = lsqlin(A,b, [] , [] , ones(1, n), 1);

sum_x_actual = sum(x_actual); 

display(sum_x_actual);

norm_error = norm(x_actual - x) ; 

display(norm_error);



relative_residual = norm(res) / norm (b) ; 
norm_res = norm(res); 
sum_val = sum(x); 

%}






