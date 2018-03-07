function [c,erro] = bissecao(a,b,tol)

fa = ele_ene(a); fb = ele_ene(b);
c = (a+b)/2; fc = ele_ene(c);
erro = abs((a-b)/2);

while(erro>tol)
	if fc == 0
		erro = 0;
		break;
	end	
	if (fa*fc < 0)
		b = c;
		fb = fc;
	else
		a = c;
		fa = fc;
	end
	c = (a+b)/2;
	fc = ele_ene(c);
	erro = abs((a-b)/2);
end
end
  
