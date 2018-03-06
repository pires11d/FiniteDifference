%DECK min3
      function [xmin]=min3(x1,x2,x3)
         xmin=x1;
         if(x2 < x1)xmin=x2; end
         if(x3 < xmin)xmin=x3; end

