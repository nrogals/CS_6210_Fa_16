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


