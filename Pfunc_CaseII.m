function P = Pfunc_CaseII(i,s,x,lambda,alpha,beta,l,omega,xi,Psi,m,inlet,outlet)
% Computes P_{i}(x,s) defined in Eqs (42)-(44) and Table 1.

a0 = inlet{1}; aL = outlet{1};

P = (omega(i)/s + alpha(i))/lambda(i);
if i == 1
    P = (1 + a0/beta(1)*(xi(1,1)*Psi(1,2,x) - xi(1,2)*Psi(1,2,l(1))*Psi(1,1,x)))*P;
elseif i == m
    P = (1 + aL/beta(m)*(xi(m,2)*Psi(m,1,x) - xi(m,1)*Psi(m,1,l(m-1))*Psi(m,2,x)))*P;
end