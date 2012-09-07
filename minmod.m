function s = minmod(a,b)

if(a*b < 0) 
    s=0;
elseif(abs(a)<abs(b))
    s=a;
else
    s=b;
end