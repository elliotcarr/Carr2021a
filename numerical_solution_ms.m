function c = numerical_solution_ms(tspan,R,D,v,M,gamma,theta,l,f,N,m,inlet,outlet,n)
% Computes a numerical solution to the multilayer transport model (1)-(6) using a finite volume
% spatial discretisation and ode15s time stepping.

% Node spacing
x = linspace(0,l(end),N)';

% Initial condition
layer = zeros(N,1);
cnt = 1;
c0 = zeros(n*N,1);
for j = 1:n
    for k = 1:N
        for i = m:-1:1
            if x(k) == l(i) && ismember(i,1:m-1)
                if j == 1, interface(cnt) = k; end
                cnt = cnt + 1;
                if (f(i,j) ~= 0) && (f(i+1,j) == 0)
                    c0((j-1)*N+k) = f(i,j);
                elseif (f(i,j) == 0) && (f(i+1,j) ~= 0)
                    c0((j-1)*N+k) = f(i+1,j);
                else 
                    c0((j-1)*N+k) = (f(i,j)+f(i+1,j))/2;
                end
            elseif x(k) < l(i)
                if j == 1, layer(k) = i; end
                c0((j-1)*N+k) = f(i,j);
            end
        end
    end
end
if length(interface) ~= m-1
    error('Choose N such that nodes are located on all interfaces');
end
h = x(2)-x(1); % uniform spacing
interior = setdiff(2:N-1,interface); % interior node indexes

% Mass matrix
Mt = sparse(n*N);  mass_singular = false;
for j = 1:n
    for k = 1:N
        Mt((j-1)*N+k,(j-1)*N+k) = 1;
    end
    b0 = inlet{j}{2};
    if b0 == 0
        Mt((j-1)*N+1,(j-1)*N+1) = 0;
        mass_singular = true;
    end
    bL = outlet{j}{2};
    if bL == 0
        Mt(j*N,j*N) = 0;
        mass_singular = true;
    end
end

% Jacobian sparsity pattern
Jp = sparse(n*N);
for j = 1:n
    % Treatment of inlet boundary condition
    b0 = inlet{j}{2};   
    if b0 == 0
        Jp((j-1)*N+1,(j-1)*N+1) = 1;
    else
        Jp((j-1)*N+1,([1:n]-1)*N+1) = 1;
        Jp((j-1)*N+1,(j-1)*N+2) = 1;        
    end
    % Treatment of outlet boundary condition
    bL = outlet{j}{2};
    if bL == 0
        Jp((j-1)*N+N,(j-1)*N+N) = 1;
    else
        Jp((j-1)*N+N,([1:n]-1)*N+N) = 1;
        Jp((j-1)*N+N,(j-1)*N+N-1) = 1;
    end
    % Interior and interface nodes
    for k = 2:N-1
        Jp((j-1)*N+k,(j-1)*N+k-1) = 1;
        Jp((j-1)*N+k,([1:n]-1)*N+k) = 1;
        Jp((j-1)*N+k,(j-1)*N+k+1) = 1;
    end
end
    
if mass_singular
    options = odeset('Mass',Mt,'MassSingular','yes','AbsTol',1e-8,'RelTol',1e-8,...
        'Jpattern',Jp);
else
    options = odeset('AbsTol',1e-8,'RelTol',1e-8,'Jpattern',Jp);
end

% Solve ODE system
F = @(t,c) Ffunc_ms(t,c,R,D,v,M,gamma,theta,m,h,layer,interface,interior,N,inlet,outlet,n);
[~,c] = ode15s(F,[0,tspan],c0,options);

% Remove initial solution
c = c(2:end,:);