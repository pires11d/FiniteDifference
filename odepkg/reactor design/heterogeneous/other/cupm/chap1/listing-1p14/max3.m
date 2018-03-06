%DECK max3
      function [xmax]=max3(x1,x2,x3)
         xmax=x1;
         if(x2 > x1)xmax=x2; end
         if(x3 > xmax)xmax=x3; end

