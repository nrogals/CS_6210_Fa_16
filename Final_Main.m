%Final Main Test Script 


n=10; 

time_vals = zeros(1, n); 
n_vals = zeros(1, n); 

for i = 1 : 4
    
u = rand (1, 2*n-1) ; 

u = u./(2*(2*n-1)); 
norm_one_u = norm(u,1) ; 

%display(norm_one_u); 

assert (norm_one_u<0.5); 

column = u(1:n) ; 

%Bottom Row of Hankel U
row = u(n:2*n-1); 

% hankel(c,r)
h= hankel(column, row); 

I = eye(n); 

A = h + I; 

b= rand(n , 1 ) ; 

x = A\b; 

tol = 1e-9; 

%final_p2solve

tic; 
[ x_final ] = final_p2solve( u , b, tol, A) ; 
time= toc; 

res = A*x_final - b ; 

norm_res = norm(res); 
norm_b = norm(b); 

relative_res_error = norm_res ./ norm_b ; 

display(x_final);

display(x); 

assert ( relative_res_error < tol); 

time_vals(i)=time; 

n=10*n; 

n_vals(i)=n; 


end

loglog(n_vals, time_vals ); 
grid on 

display(n_vals); 
display(time_vals); 

