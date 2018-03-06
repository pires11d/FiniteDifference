      function [ux]=super(xl,xu,n,u,v)
      dx=(xu-xl)/(n-1);
      delta=1.0e-05;
      if v >= 0.0
        for i=3:n-1
          if(abs(u(i)-u(i-1))<delta)
            phi(i)=0.0;
          else
            r(i)=(u(i+1)-u(i))/(u(i)-u(i-1));
            phi(i)=max3(0.0,min(2.0*r(i),1.0),min(r(i),2.0));
          end
          if(abs(u(i-1)-u(i-2))<delta)
            phi(i-1)=0.0;
          else
            r(i-1)=(u(i)-u(i-1))/(u(i-1)-u(i-2));
            phi(i-1)=max3(0.0,min(2.0*r(i-1),1.0),min(r(i-1),2.0));
          end
          flux2=u(i  )+(u(i  )-u(i-1))*phi(i  )/2.0;
          flux1=u(i-1)+(u(i-1)-u(i-2))*phi(i-1)/2.0;
          ux(i)=(flux2-flux1)/dx;
        end
        ux(1)=(-u(1)+u(2))/dx;
        ux(2)=(-u(1)+u(2))/dx;
        ux(n)=(u(n)-u(n-1))/dx;
      end
      if v < 0.0
        for i=2:n-2
          if(abs(u(i)-u(i+1))<delta)
            phi(i)=0.0;
          else
            r(i)=(u(i-1)-u(i))/(u(i)-u(i+1));
            phi(i)=max3(0.0,min(2.0*r(i),1.0),min(r(i),2.0));
          end
          if(abs(u(i+1)-u(i+2))<delta)
            phi(i+1)=0.0;
          else
            r(i+1)=(u(i)-u(i+1))/(u(i+1)-u(i+2));
            phi(i+1)=max3(0.0,min(2.0*r(i+1),1.0),min(r(i+1),2.0));
          end
          flux2=u(i  )+(u(i  )-u(i+1))*phi(i  )/2.0;
          flux1=u(i+1)+(u(i+1)-u(i+2))*phi(i+1)/2.0;
          ux(i)=-(flux2-flux1)/dx;
        end
      ux(1)=(u(2)-u(1))/dx;
      ux(n-1)=(-u(n-1)+u(n))/dx;
      ux(n  )=(-u(n-1)+u(n))/dx;
      end
