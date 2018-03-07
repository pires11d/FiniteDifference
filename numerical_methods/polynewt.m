function p = polynewt(x,c,t)

n = max(size(x)) - 1;
p = c(n+1);

for i=n:-1:1
p = p.*(t-x(i)) + c(i);
end
