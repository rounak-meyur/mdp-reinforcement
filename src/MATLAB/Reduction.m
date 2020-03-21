classdef Reduction
  % Reduction functions for a grid
   
   properties
   end
   
   methods(Static)
      function [new_rnc] = remove_parralel(rnc)
         x = rnc.branch(:,4); new_rnc = rnc;
         E = Reduction.check_parralel(rnc.branch(:,1:2));
         x(:,2)=ones(Sz.r(x),1);
         for i=1:Sz.r(E)
            if (E(i,3)~=i)
               ASSERT(E(i,3)<i);
               x(E(i,3),1)=x(E(i,3),1)+x(i,1);
               x(E(i,3),2)=x(E(i,3),2)+1;
               E(i,3) = 0; x(i,1)=0;
            end
         end
         x(:,1)=x(:,1)./x(:,2);
         new_rnc.branch(E(:,3)==0,:)=[]; 
      end
      function [structure] = check_parralel(E) % method checks if there are parralel lines and gives back appropriate massive
         n = Sz.r(E); % number of lines
         for i = 1:n
            similar1 = find((E(1:(i-1),2)==E(i,2))&(E(1:(i-1),1)==E(i,1)));
            similar2 = find((E(1:(i-1),2)==E(i,1))&(E(1:(i-1),1)==E(i,2)));
            if (~Sz.z(similar1))
               if (~Sz.z(similar2))
                  E(i,3)=min(similar1(1),similar2(1));
               else
                  E(i,3)=similar1(1);
               end
            else
               if (~Sz.z(similar2))
                  E(i,3)=similar2(1);
               else
                  E(i,3)=i;
               end
            end
         end
         structure = E; % old E + column of remapped lines 
      end
      function [new_E,new_C,new_P,new_g,new_x,new_lim,new_oldf,rnc] = remove_leafes(E,C,P,g,x,lim,oldf,rnc) % g - ground bus
         w = zeros(Sz.c(C),1); P(:,2)=P(:,1); E(:,3)=ones(Sz.r(E),1); E(:,4)=1:Sz.r(E);
         for i=1:Sz.c(C)
            w(i)=sum(abs(C(:,i))); % number of neighbours of each bus
         end
         while (sum(w(:,1)==1)>0) % While we have buses that have 1 neighbour
            for i=1:Sz.r(w)
               if (w(i)==1)
                  if (i==565)
                     q=1;
                  end
                  line = find((abs(C(:,i))==1)&(E(:,3)~=0));
                  E(line,3)=0;
                  if (E(line,1)==i)
                     if (i==g)  % change ground bus if it becomes reduced
                        g=E(line,2);
                     end
                     w(E(line,2))=w(E(line,2))-1;
                     P(E(line,2),2)=P(E(line,2),2)+P(i,2);
                  else
                     if (i==g)
                        g=E(line,1);
                     end
                     w(E(line,1))=w(E(line,1))-1;
                     P(E(line,1),2)=P(E(line,1),2)+P(i,2);
                  end
                  
                  w(i)=0;
               end
            end
         end
         w(:,2)=1:Sz.r(w);
         w(w(:,1)==0,:)=[];
         lines = find(E(:,3)~=0);
         E(E(:,3)==0,:)=[];
         for i=(1:Sz.r(E))
            E(i,1)=find(w(:,2)==E(i,1));
            E(i,2)=find(w(:,2)==E(i,2));
         end
         
         new_g=find(w(:,2)==g);
         new_oldf = oldf(lines(:));
         new_lim = lim(lines(:));
         new_E = E(:,1:2);
         new_P = P(w(:,2),2);
         % rnc.bus = rnc.bus(w(:,2),:);
         % rnc.branch = rnc.branch(lines(:)); rnc.branch(:,9) = new_lim;
         new_C = C(lines(:),w(:,2));
         new_x = x(lines(:));
      end
      function [area] = connectivity_analysis(E) % check if graph is connected
         % E should be such that there is at least one line for each bus
         area=1:max(E(:));
         for i=(1:Sz.r(E))
            if area(E(i,1))~=area(E(i,2))
               qM = max(area(E(i,1)),area(E(i,2)));
               qm = min(area(E(i,1)),area(E(i,2)));
               area((area(:)==qM))=qm;
            end
         end
      end
      function [rnc,map] = remap_grid(rnc) % remaps buses in the grid
         % because some cases have messed up numbering
         map(:,1) = rnc.bus(:,1);
         map(:,2) = 1:Sz.r(rnc.bus);
         for i=1:Sz.r(rnc.gen)
            rnc.gen(i,1) = map((map(:,1)==rnc.gen(i,1)),2);
         end
         for i=1:Sz.r(rnc.branch)
            rnc.branch(i,1) = map((map(:,1)==rnc.branch(i,1)),2);
            rnc.branch(i,2) = map((map(:,1)==rnc.branch(i,2)),2);
         end
         rnc.bus(:,1) = map(:,2);
      end
      
      %function [rnc] = aggregate_generation(rnc)
      %   j=1;
      %   while (j<Sz.r(rnc.gen))
      %      i=j+1;
      %      while (i<=Sz.r(rnc.gen)) 
      %         if (rnc.gen(j,1)==rnc.gen(i,1))
      %            rnc.gen(j,[2:5,9:10])=rnc.gen(j,[2:5,9:10])+rnc.gen(i,[2:5,9:10]);
      %            rnc.gen(i,:)=[];
      %            rnc.gencost(i,:)=[];
      %         else
      %            i=i+1;
      %         end
      %      end
      %      j=j+1;
      %   end
      %end

end
   
end

