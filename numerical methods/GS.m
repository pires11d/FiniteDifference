function x = GS(A,b,tol,kmax)
m = size(A,1);
k = 0;
erro = 2*tol;

for i=1:m
xk(i) = 0;
end

while k<kmax
	while erro>=tol
	k = k+1;			
		for i=1:m
		soma1(i) = 0;
			for j=1:m
				if (j~=i)
				soma1(i) = soma1(i) + A(i,j).*xk(j);
				end
		x(i) = (b(i)-soma1(i))./A(i,i);
		d(i) = abs(x(i)-xk(i));	
			end
		end
	erro = sum(d);
	xk = x;
		if (erro<=tol)
		break
		end
	end
if (k<=kmax)
break
end
end

k
erro

end
