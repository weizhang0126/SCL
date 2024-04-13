function x=encode(u)
N=length(u);
if N==1
    x=u;
else
    u1u2=mod(u(1:N/2)+u(N/2+1:N),2);
    u2=u(N/2+1:N);
   x=[encode(u1u2) encode(u2)];
end
end
