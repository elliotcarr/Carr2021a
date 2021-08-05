function B = Bfunc(i,x,Psi,m,D,theta,beta,lambda,l,inlet)
% Computes B_{i}(x,s) defined in Eqs (31)-(33) and Table 1.

a0 = inlet{1}; b0 = inlet{2};

if i == 1
    B = 1/(theta(1)*D(1)*beta(1))*((a0-b0*lambda(1,1))*Psi(1,1,0)*Psi(1,2,x)...
            - (a0-b0*lambda(1,2))*Psi(1,1,x));
elseif i == m
    B = 1/(beta(m))*(lambda(m,1)*Psi(m,1,l(m-1))*Psi(m,2,x) - lambda(m,2)*Psi(m,1,x));
else
    B = 1/(theta(i)*D(i)*beta(i))*(lambda(i,1)*Psi(i,1,l(i-1))*Psi(i,2,x) - lambda(i,2)*Psi(i,1,x));
end