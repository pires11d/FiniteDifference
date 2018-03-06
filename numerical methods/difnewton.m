function c = difnewton(x,fx)

n = max(size(x)) - 1;
a = zeros(n+1,n+1);

for i=1:n+1
	a(i,1) = fx(i);
end

for j = 2:n+1
	for i=1:(n+1)-(j-1)
		a(i,j) = (a(i+1,j-1) - a(i,j-1))/(x(i+j-1) - x(i));
	end    
end
a
c = a(1,:);
end





