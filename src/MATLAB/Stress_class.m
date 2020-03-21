classdef Stress_class
   %UNTITLED Summary of this class goes here
   %   Detailed explanation goes here
   
   properties
   end
   
   methods(Static)
      function [rnc] = stress (rnc,coef)
         rnc.bus(:,3) = rnc.bus(:,3)*coef;
         rnc.gen(:,9) = rnc.gen(:,9)*coef;
      end
      function [suc_coef] = maximum_stress(rnc)
         coef = 1; % stress coefficient
         k=0; k_max = 10; % maximum number of iterations
         success = 1; rnc0 = rnc; fal_coef = 0;
         while k<k_max
            k=k+1; rnc = rnc0;
            if (success == 1)
               suc_coef = coef;
               if (fal_coef == 0)
                  coef = coef*2;
               else
                  coef = (suc_coef+fal_coef)/2;
               end
            else
               fal_coef = coef;
               coef = (coef+suc_coef)/2;
            end
            rnc = Stress_class.stress (rnc, coef);
            mpopt = mpoption('OUT_ALL',0,'VERBOSE',0);
            [rnc, success] = rundcopf(rnc,mpopt);
         end
         if (success == 1)
            suc_coef = coef;
         end
      end
   end
   
end

