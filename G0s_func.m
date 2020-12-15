function G0t = G0s_func(j,s,inlet,Qinv,n)

G0t = 0;
for i = 1:n
    G0t = G0t + Qinv(j,i)*inlet{i}{4}(s);
end