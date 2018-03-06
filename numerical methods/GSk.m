function x = GSk(A,b,tol,kmax)
m = size(A,1);
k = 0;
erro = 2*tol;

for i=1:m
x(i) = 0;
end

while k<kmax
k = k+1;			
	for i=1:m
	soma1(i) = 0;
		for j=1:m
			if (j~=i)
			soma1(i) = soma1(i) + A(i,j).*x(j);
			end
		end
	x(i) = (b(i)-soma1(i))./A(i,i);
	end
end

k
erro

end
