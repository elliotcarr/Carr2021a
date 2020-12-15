function g0t = g0t_func(j,t,inlet,Qinv,n)

g0t = 0;
for i = 1:n
    g0t = g0t + Qinv(j,i)*inlet{i}{3}(t);
end
