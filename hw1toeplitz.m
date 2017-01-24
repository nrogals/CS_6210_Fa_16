
function [y] = hw1toeplitz(c,r,x)

Fa = fft([c; flipud(r(2:end))]);
xp = fft([x; zeros(length(r)-1,1)]);
yp = ifft(Fa .* xp);
y = yp(1:length(x));

end

