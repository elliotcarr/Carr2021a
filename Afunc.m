function A = Afunc(i,x,Psi,m,D,theta,beta,lambda,l,outlet)
% Computes A_{i}(x,s) defined in Eqs (31)-(33) and Table 1.

aL = outlet{1}; bL = outlet{2};

if i == 1
    A = 1/(beta(1))*(lambda(1,2)*Psi(1,2,l(1))*Psi(1,1,x) - lambda(1,1)*Psi(1,2,x));
elseif i == m
    A = 1/(theta(m)*D(m)*beta(m))*((aL + bL*lambda(m,2))*Psi(m,2,l(m))*Psi(m,1,x) ...
        - (aL + bL*lambda(m,1))*Psi(m,2,x));
else
    A = 1/(theta(i)*D(i)*beta(i))*(lambda(i,2)*Psi(i,2,l(i))*Psi(i,1,x) - lambda(i,1)*Psi(i,2,x));
end