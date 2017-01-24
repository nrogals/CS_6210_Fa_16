%Test Toeplitz 

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

