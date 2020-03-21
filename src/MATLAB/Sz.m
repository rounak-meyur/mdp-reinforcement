classdef Sz
   % size class makes comfortable to calculate row and columns of a matrix
   methods (Static)
      function r = r(matrix) % rows of a matrix
         D=size(matrix);
         r=D(1);
      end
      function c = c(matrix) % columns of a matrix
         D=size(matrix);
         c=D(2);
      end
      function z = z(matrix) % zero or no
         z = (Sz.r(matrix)==0);
      end
   end
end

