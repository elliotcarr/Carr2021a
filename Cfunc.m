function C = Cfunc(s,x,v,D,theta,mu,R,gamma,l,f,inlet,outlet)
% Evaluates the Laplace-domain concentration from Eqs (31)-(33) given values of x and s.

m = length(D);
D = reshape(D,m,1); 
v = reshape(v,m,1); 
mu = reshape(mu,m,1); 
R = reshape(R,m,1); 
gamma = reshape(gamma,m,1); 
f = reshape(f,m,1); 

% Lambda values from Table 1
lambda(:,1) = (v+sqrt(v.^2+4*D.*(R*s+mu)))./(2*D);
lambda(:,2) = (v-sqrt(v.^2+4*D.*(R*s+mu)))./(2*D);

% Boundary condition parameters
a0 = inlet{1}; b0 = inlet{2}; G0 = inlet{4};
aL = outlet{1}; bL = outlet{2}; GL = outlet{4};

% Beta values from Table 1
beta = zeros(m,1);
beta(1) = (a0 - b0*lambda(1,1))*lambda(1,2)*exp(-(lambda(1,1)-lambda(1,2))*l(1)) ...
    - (a0 - b0*lambda(1,2))*lambda(1,1);
for i = 2:m-1
    beta(i) = lambda(i,1)*lambda(i,2)*(exp(-(lambda(i,1)-lambda(i,2))*(l(i)-l(i-1)))-1);
end
beta(m) = (aL + bL*lambda(m,2))*lambda(m,1)*exp(-(lambda(m,1)-lambda(m,2))*(l(m)-l(m-1))) ...
    - (aL + bL*lambda(m,1))*lambda(m,2);

% Functions from Table 1
Psi = @(i,j,x) Psifunc(i,j,x,l,lambda);
P = @(i,s,x) Pfunc(i,s,x,mu,R,gamma,f,l,beta,lambda,Psi,m,inlet,outlet);
A = @(i,x) Afunc(i,x,Psi,m,D,theta,beta,lambda,l,outlet); 
B = @(i,x) Bfunc(i,x,Psi,m,D,theta,beta,lambda,l,inlet);

if m == 2 % two layers
    
    Amat = B(1,l(1)) - A(2,l(1));
    bmat = P(2,s,l(1)) - P(1,s,l(1)) - A(1,l(1))*G0(s) + B(2,l(1))*GL(s);
    G = bmat/Amat;
    
else % three or more layers
    
    % Build linear system given in Eq (35)
    Amat = zeros(m-1,m-1); bvec = zeros(m-1,1);
    for i = 2:m-2
        Amat(i,i-1) = A(i,l(i));
        Amat(i,i) = B(i,l(i)) - A(i+1,l(i));
        Amat(i,i+1) = -B(i+1,l(i));
        bvec(i) = P(i+1,s,l(i)) - P(i,s,l(i));
    end
    
    Amat(1,1) = B(1,l(1)) - A(2,l(1));
    Amat(1,2) = -B(2,l(1));
    bvec(1) = P(2,s,l(1)) - P(1,s,l(1)) - A(1,l(1))*G0(s);
    
    Amat(m-1,m-1) = B(m-1,l(m-1)) - A(m,l(m-1));
    Amat(m-1,m-2) = A(m-1,l(m-1));
    bvec(m-1) = P(m,s,l(m-1)) - P(m-1,s,l(m-1)) + B(m,l(m-1))*GL(s);
    
    G = Amat\bvec;
    
end

% Determine layer where value of x is located
for i = m:-1:1
    if x <= l(i)
        layer = i;
    end
end

% Laplace-domain concentration from Eqs (31)-(33)
i = layer;
if i == 1
    C = P(1,s,x) + A(1,x)*G0(s) + B(1,x)*G(1);
elseif i == m
    C = P(m,s,x) + A(m,x)*G(m-1) + B(m,x)*GL(s);
else
    C = P(i,s,x) + A(i,x)*G(i-1) + B(i,x)*G(i);
end