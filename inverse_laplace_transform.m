function invL = inverse_laplace_transform(g,t,N,c,z)
% Numerical inversion of the Laplace transform as described in Eq (36).

invL = 0.0;
for k = 1:2:N
    s = z(k)/t;
    invL = invL-c(k)*g(s);
end
invL = 2*real(invL)/t;
% [z,indx] = sort(z);
% c = c(indx);

% invL = 0.0;
% for k = 1:N
%     s = z(k)/t;
%     invL = invL-c(k)*g(s);
% end
% invL = real(invL)/t;