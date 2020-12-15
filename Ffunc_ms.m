function F = Ffunc_ms(t,c,R,D,v,M,gamma,theta,m,h,layer,interface,interior,N,inlet,outlet,n)
% Right-hand side function appearing FVM discretisation: M dc/dt = F(c)

F = zeros(N,n);
c = reshape(c,N,n);

for j = 1:n
    
    % Treatment of inlet boundary condition
    a0 = inlet{j}{1};
    b0 = inlet{j}{2};
    g0 = inlet{j}{3};
    
    if b0 == 0
        F(1,j) = a0*c(1,j) - g0(t);
    else
        r = M(j,:)*c(1,:)';
        F(1,j) = D(1)*(c(2,j)-c(1,j))/h - v(1)*(c(1,j)+c(2,j))/2 + D(1)/b0*g0(t) + ...
            (v(1)-D(1)*a0/b0)*c(1,j) + h/2*(r+gamma(1,j));
        F(1,j) = F(1,j)/(h/2*R(1,j));
    end
    
    % Treatment of outlet boundary condition    
    aL = outlet{j}{1};
    bL = outlet{j}{2};
    gL = outlet{j}{3};
    
    if bL == 0
        F(N,j) = aL*c(N,j) - gL(t);
    else
        r = M(j,:)*c(N,:)';
        F(N,j) = D(m)/bL*gL(t) - (v(m) + D(m)*aL/bL)*c(N,j) - D(m)*(c(N,j)-c(N-1,j))/h + ...
            v(m)*(c(N-1,j)+c(N,j))/2 + h/2*(r+gamma(m,j));
        F(N,j) = F(N,j)/(h/2*R(m,j));
    end
    
    % Interface nodes
    for i = 1:m-1 % node located at interface between layers i and i+1        
        k = interface(i);
        r = M(j,:)*c(k,:)';
        F(k,j) = theta(i+1)*D(i+1)*(c(k+1,j)-c(k,j))/h - theta(i+1)*v(i+1)*(c(k,j)+c(k+1,j))/2 - ...
            theta(i)*D(i)*(c(k,j)-c(k-1,j))/h + v(i)*theta(i)*(c(k-1,j)+c(k,j))/2 + ...
            h/2*(theta(i)*(r+gamma(i,j))+theta(i+1)*(r+gamma(i+1,j)));
        F(k,j) = F(k,j)/(h/2*(theta(i)*R(i,j)+theta(i+1)*R(i+1,j)));
        
    end
    
    % Interior nodes
    for k = interior
        i = layer(k); % node located in interior of layer i
        r = M(j,:)*c(k,:)';
        F(k,j) = D(i)*(c(k+1,j)-c(k,j))/h - v(i)*(c(k,j)+c(k+1,j))/2 - D(i)*(c(k,j)-c(k-1,j))/h ...
            + v(i)*(c(k-1,j)+c(k,j))/2 + h*(r+gamma(i,j));
        F(k,j) = F(k,j)/(h*R(i,j));
    end
    
end

F = reshape(F,N*n,1);