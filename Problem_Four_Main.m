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








%{
n=10; 

P = rand (n); 

[W,S,V] = svd(P); 



D = diag ([1 , 1 , 1 , 1 , -1 , -1 , -1 ,-1 , -1 , -1]); 

%display(D); 

B= rand (n); 

X= triu(B , 1);

Schur_T = X + D ; 



A = W * Schur_T * transpose(W); 

N= A*A ; 
display(N);

assert(1<0);


[U_Schur,T_Schur] = schur(A); 

C = rand (n); 

[R,I,E] = svd(C);
%Attain the mixer orthogonal matrix


U = W(:, 1:4);
orthog_V = W(:, 5:end); 


[V , J ] = jordan(A); 

display(V_jor); 

assert(1<0); 

vectors = V(: , 5:end); 

display(V); 



assert(1<0); 




[ Q , T ] = final_p4solve( U, mixed_V ) ; 

display(Q); 
display(orthog_V); 

display *************** 

display(T); 
display(X); 


%}




























