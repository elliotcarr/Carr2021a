function C = Cfunc_wrap(s,x,v,D,theta,M,Rspecies,gamma,l,f,inlet,outlet,n,m)
% Returns vector of Laplace transformed solutions for all species at given
% values of x and s.

C = zeros(n,1); % Concentration in Laplace domain for each species
Ct = zeros(n,1); % Transformed concentration in Laplace domain for each species

%% Apply species-decoupling transformation
R = diag(Rspecies);
[Q,Lambda] = eig(M - s*R);
lambda = diag(Lambda);
alpha = -Q\(R*(f'));
omega = -Q\(gamma');

inlett = cell(n,1);
outlett = cell(n,1);
% ut = Q\eye(n,1);
% for i = 1:n
%     inlett{i} = {inlet{i}{1}, inlet{i}{2}, @(t) ut(i)*inlet{1}{3}(t), @(s) ut(i)*inlet{1}{4}(s)};
% end
Qinv = Q\eye(n);
for j = 1:n
    inlett{j} = cell(4,1);
    inlett{j}{1} = inlet{j}{1};
    inlett{j}{2} = inlet{j}{2};
    inlett{j}{4} = @(s) G0s_func(j,s,inlet,Qinv,n);
    outlett{j} = cell(4,1);
    outlett{j}{1} = outlet{j}{1};
    outlett{j}{2} = outlet{j}{2};
    outlett{j}{4} = @(s) G0s_func(j,s,outlet,Qinv,n);    
end

%% Solve transformed problem in Laplace domain
for j = 1:n
    Ct(j) = Cfunc_CaseII(s,x,v,D,theta,lambda(j)*ones(1,m),alpha(j,:),l,omega(j,:),inlett{j},outlett{j});
end

%% Reverse species-decoupling transformation
for i = 1:n
    for j = 1:n
        C(i) = C(i) + Q(i,j)*Ct(j);
    end
end