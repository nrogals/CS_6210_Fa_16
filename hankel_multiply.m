function H_x = hankel_multiply(row , column , x )
%Form the Hankel Matrix given by H(row , column) and multiply it by x using
%Toeplitz multiply. 



row = transpose(row); 
column=transpose(column);

flip_column = flipud (column) ; 



x_flipped = hw1toeplitz ( flip_column, row , x) ; 

H_x = flipud(x_flipped); 

end


