function y = Psifunc(i,j,x,l,lambda)
% Computes Psi_{i,j}(x,s) for j = 1,2 defined in Table 1.

if j == 1
    y = exp(lambda(i,1)*(x-l(i)));
%     Re = real(lambda(i,1)*(x-l(i)));
%     Im = imag(lambda(i,1)*(x-l(i)));
%     y = exp(Re)*(cos(Im) + 1i*sin(Im));
else
    if i == 1
        y = exp(lambda(i,2)*x);
%         Re = real(lambda(i,2)*x);
%         Im = imag(lambda(i,2)*x);
%         y = exp(Re)*(cos(Im) + 1i*sin(Im));
    else
        y = exp(lambda(i,2)*(x-l(i-1)));
%         Re = real(lambda(i,2)*(x-l(i-1)));
%         Im = imag(lambda(i,2)*(x-l(i-1)));
%         y = exp(Re)*(cos(Im) + 1i*sin(Im));
    end
end