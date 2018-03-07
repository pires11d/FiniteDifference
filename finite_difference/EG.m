function x = EG(A,b)

n = size(A, 1);

%Etapa de Eliminação
for k=1:(n-1)
	for i=k+1:n
		m = A(i,k)/A(k,k);
		for j=k:n
			A(i,j) = A(i,j) - m*A(k,j);
		end			
		b(i) = b(i) - m*b(k);	
	end
end


%Etapa de Retrossubstituição
for i=n:-1:1
	soma = 0;
	for j=i+1:n
		soma = soma + A(i,j)*x(j);
	end
	x(i) = (b(i) - soma)/ A(i,i);
end
x = x';
