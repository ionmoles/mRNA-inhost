function F = AB_fun(t,a,alpha)

B = 0;
A = 0;
for k=1:3
    B = B + a(k)*exp(-alpha(k)*t);
    if alpha(k)==1
        A = A + a(k)*t*exp(-t);
    else
        A = A + a(k)/(1-alpha(k))*(exp(-alpha(k)*t)-exp(-t));
    end
end

F = B-A;

end