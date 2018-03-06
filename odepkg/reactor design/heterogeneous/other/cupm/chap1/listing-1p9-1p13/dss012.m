      function [ux]=dss012(xl,xu,n,u,v)
%...
%...  FUNCTION DSS012 IS AN APPLICATION OF FIRST-ORDER DIRECTIONAL
%...  DIFFERENCING IN THE NUMERICAL METHOD OF LINES.  IT IS INTENDED
%...  SPECIFICALLY FOR THE ANALYSIS OF CONVECTIVE SYSTEMS MODELLED BY
%...  FIRST-ORDER HYPERBOLIC PARTIAL DIFFERENTIAL EQUATIONS WITH THE
%...  SIMPLEST FORM
%...
%...                            U  + V*U  = 0                        (1)
%...                             T      X
%...
%...  THE FIRST FOUR PARAMETERS, XL, XU, N AND U, ARE THE SAME AS
%...  FOR FUNCTION DSS002.  THE FIFTH PARAMETER, V, MUST BE PROVIDED
%...  TO DSS012 SO THAT THE DIRECTION OF FLOW IN EQUATION (1) CAN BE
%...  USED TO SELECT THE APPROPRIATE FINITE DIFFERENCE APPROXIMATION
%...  FOR THE FIRST-ORDER SPATIAL DERIVATIVE IN EQUATION (1), U .
%...  THE CONVENTION FOR THE SIGN OF V IS                      X
%...
%...     FLOW LEFT TO RIGHT             V GT 0
%...     (I.E., IN THE DIRECTION        (I.E., THE SIXTH ARGUMENT IS
%...     OF INCREASING X)               POSITIVE IN CALLING DSS012)
%...
%...     FLOW RIGHT TO LEFT             V LT 0
%...     (I.E., IN THE DIRECTION        (I.E., THE SIXTH ARGUMENT IS
%...     OF DECREASING X)               NEGATIVE IN CALLING DSS012)
%...
%...  COMPUTE THE SPATIAL INCREMENT, THEN SELECT THE FINITE DIFFERENCE
%...  APPROXIMATION DEPENDING ON THE SIGN OF V IN EQUATION (1).
      dx=(xu-xl)/(n-1);
      if v > 0
%...
%...     (1)  FINITE DIFFERENCE APPROXIMATION FOR POSITIVE V
              ux(1)=(u(2)-u(1))/dx;
              for i=2:n
                 ux(i)=(u(i)-u(i-1))/dx;
              end
      end
%...
%...     (2)  FINITE DIFFERENCE APPROXIMATION FOR NEGATIVE V
      if v < 0
              nm1=n-1;
              for i=1:nm1
                 ux(i)=(u(i+1)-u(i))/dx;
              end
              ux(n)=(u(n)-u(n-1))/dx;
      end
