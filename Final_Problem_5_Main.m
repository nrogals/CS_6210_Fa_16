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




