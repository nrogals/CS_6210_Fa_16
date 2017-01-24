function [ x ] = final_p2solve( u , b, tol )
%Solve Ax=b to a relative residual of at most tol

n = max (size ( b)) ; 

%First column of Hankel U
column = u(1:n) ; 
%Bottom Row of Hankel U
row = u(n:2*n-1); 
I = speye ( n ) ; 
U_diag = sparse_diagonal (u , n); 

D = I + U_diag ; 



residual_condition_not_met =  1 ; 

%Intialize guess
x = rand (n , 1 ) ; 

normb = norm(b, Inf); 


while residual_condition_not_met 
    U_x = hankel_multiply( row , column , x) ; 
    
    r = b - (x + U_x); 
    
    normres = norm (r, Inf); 
    
    relative_residual_error = normres ./ normb ; 
    
    if relative_residual_error < tol
        residual_condition_not_met = 0 ;
        x = x +  D\r; 
    else
         x = x +  D\r; 
    end
    
end

x_final = x  ;

end


function sdiag_from_u = sparse_diagonal (u , n ) 


x_index_list = zeros (1 , n ); 
y_index_list = zeros (1 , n );  

value_list = zeros (1 , n );  

for i = 1 : n 
    val = u(2*i - 1) ;
    value_list (i) = val ; 
    x_index_list (i) = i ; 
    y_index_list (i) = i ; 
end

sdiag_from_u = sparse(x_index_list,y_index_list, value_list ,n,n) ; 


end 


function H_x = hankel_multiply(row , column , x )
%Form the Hankel Matrix given by H(row , column) and multiply it by x using
%Toeplitz multiply. 
row = transpose(row); 
column=transpose(column);

flip_column = flipud (column) ; 

x_flipped = hw1toeplitz ( flip_column, row , x) ; 

H_x = flipud(x_flipped); 

end



function [y] = hw1toeplitz(c,r,x)

Fa = fft([c; flipud(r(2:end))]);
xp = fft([x; zeros(length(r)-1,1)]);
yp = ifft(Fa .* xp);
y = yp(1:length(x));

end



