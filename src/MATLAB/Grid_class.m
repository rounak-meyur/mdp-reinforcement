classdef Grid_class
   %   GRID_CLASS is class for grid object, which contains information 
   %   about the graph of the grid as well as some methods on the grid
   
   properties(GetAccess = 'public', SetAccess = 'public')
      N; % number of nodes
      M; % number of lines
      G; % number of generators
      rnc; % non-parsed case
      x; % array of conductivities
      f; % array of power flows
      E; % connectivity matrix ( look up real name )
      A; % adjacency matrix
      P; % vector of powers
      C; % c matrix (connection matrix)
      B; % susceptance matrix
      lim; % vector of limits
      L; % LODF matrix
      Mask;
      filtered_size; % size of the filtered set after fast algorithm
      C1_isl; % set of signle islanding outages
      C2_isl; % set of double islanding outages
      gr_bus; % number of ground bus
      oldf; % buffer variable to remember flows on previous step
      brute_cont; % Brute-force calculated contingencies
      t_fast; % fast algorithm completion time
      t_brute; % brute force completion time 
      brute_cont_fake; % Brute-force fake contingencies ( -.-.-.-.- )
      case_name; % name of the case
      map; % mapping of old (messed up numbering and new correct one)
      error='';
      err=0;
      path='';
   end
   
   methods  
      function path = setpath(obj)
         fullpath = fileparts(mfilename('fullpath')); path = fullpath;
%          cd(fullpath); cd('../'); fullpath=pwd;
%          path40 = strcat(fullpath,'/matpower4_0');path41 = strcat(fullpath,'/matpower4_1');
%          addpath(path41); addpath(path41);
%          if (version==0), rmpath(path41);
%          else rmpath(path40); end      
      end
      
      function obj = Grid_class(rc,case_name) % Class Constructor
         obj.case_name = case_name;
         obj.path = obj.setpath();
         [rc,obj.map] = Reduction.remap_grid(rc); % correct numbering of buses in the grid 
         rc.branch(rc.branch(:,11)==0,:)=[]; % delete "off" lines
         rc.branch(:,9)=1; rc.branch(:,10)=0; % delete information about phase shifter and transformators
         obj.rnc = Reduction.remove_parralel(rc); % remove all parallel lines
         obj = obj.dcopf();

         obj.E = obj.rnc.branch(:,1:2); % Connectivity matrix
         obj.P = obj.create_P(); % Power vector
         obj.gr_bus = find(obj.rnc.bus(:,2)==3); % Ground bus
         obj.x = obj.rnc.branch(:,4); % vector of conductivities
         obj.C = obj.create_C(); % C matrix (each row - one line that has 1 at 'from' bus and -1 at 'to' bus)
         obj.oldf = obj.rnc.branch(:,14); % Power flows calculated by MATPOWER
         obj.lim = obj.rnc.branch(:,6); % Limits
                  
         % Remove leafes
         %[obj.E,obj.C,obj.P,obj.gr_bus,obj.x,obj.lim,obj.oldf,obj.rnc] = Reduction.remove_leafes(obj.E,obj.C,obj.P,obj.gr_bus,obj.x,obj.lim,obj.oldf,obj.rnc);
         
         % Calculate essential matrices
         obj.B = obj.create_B(); % Susceptance matrix
         obj = obj.ground_bus(); % Ground a ground bus
         obj.f = diag(1./obj.x)*obj.C*(obj.B\obj.P); % calculate flow vector
      
      end
      function P = create_P(obj) % collects P vector
         P = zeros(Sz.r(obj.rnc.bus),1);
         for i=1:Sz.r(P)
            P(i)=-obj.rnc.bus(i,3);
         end
         for i=1:Sz.r(obj.rnc.gen)
            P(obj.rnc.gen(i,1))=P(obj.rnc.gen(i,1))+obj.rnc.gen(i,2);
         end
      end
      function C = create_C(obj) % C matrix
         C = zeros(Sz.r(obj.E),Sz.r(obj.rnc.bus));
         for i=(1:Sz.r(obj.E))
            C(i,obj.E(i,1))=1;  % from
            C(i,obj.E(i,2))=-1; % to
         end
      end
      function B = create_B(obj) % B matrix
         B = obj.C'*diag(1./obj.x)*obj.C;
      end
      function obj = ground_bus(obj) % grounding one bus to make susceptance matrix invertible
         obj.C(:,obj.gr_bus)=[];
         obj.B(:,obj.gr_bus)=[];
         obj.B(obj.gr_bus,:)=[];
         obj.P(obj.gr_bus)=[];
      end
      
      function obj = dcopf(obj) % runs Optimal Power Flow from MATPOWER package
         mpopt = mpoption('OUT_ALL',0,'VERBOSE',0);
         fname = strcat(obj.path,'/logs/log_last_opf');
         solvedcase = strcat(obj.path,'/logs/opf_solved_case');
         [obj.rnc, success] = rundcopf(obj.rnc,mpopt,fname,solvedcase);         
         if (success==1), fprintf('OPF was successfuly solved');
         else fprintf('Error calculating OPF\n');
              obj.error='OPF error'; obj.err=1;
         end
      end
      function [number_of_violations,margin_absolute,margin_relative,top]=N_0_analysis(obj)
         [number_of_violations] = sum(abs(obj.f(:))>obj.lim(:));
         [margin_absolute] = min(obj.lim(:)-abs(obj.f(:)));        
         [margin_relative] = (obj.lim(:)-abs(obj.f(:)))./obj.lim(:);
         [dangerous,ind] = sort((obj.lim(:)-abs(obj.f(:)))./obj.lim(:));
         [top]=[dangerous(1:10),ind(1:10)];
      end
      function obj = N_1_analysis_alternative(obj)
         Grid_class.beauty_print('   Start N-1 analysis  ');
         invB = inv(obj.B);
         reverseStr=''; k=0; lines = zeros(Sz.r(obj.E),1);
         margins = zeros(Sz.r(obj.E),1); % vector of maximum violations
         L = zeros(Sz.r(obj.E),Sz.r(obj.E));
         C1 = diag(1./obj.x)*obj.C; % calculate this matrix once for better performance
         C1_isl = zeros(Sz.r(obj.E),1); % vector of radial contingencies
         for i=1:Sz.r(obj.E)
            % Turn off one line, invert matrix
            % using Woodbury formula inv(A+U'CU) =
            % inv(A)-inv(A)*U'*inv(inv(C)+U*inv(A)*U')*U*inv(A)
            newInvB = Grid_class.woodbury_inverse(invB,obj.C(i,:),1/obj.x(i));
            if (newInvB ~= Inf)
               flows = C1*(newInvB*obj.P);flows(i)=0;
               if (obj.f(i)==0)
                  L(:,i)=0; L(i,i) = -1;
               else
                  L(:,i)=(flows-obj.f)/obj.f(i); % LODF
               end
               qq = abs(flows(:))./obj.lim(:);
               number_of_violations = sum(abs(flows)>obj.lim(:));
               margins = max(margins(:),qq(:)); % Update maximum violations
               lines(abs(flows)>obj.lim(:))=lines(abs(flows)>obj.lim(:))+1; % mark dangerous lines
               if (number_of_violations>0)
                  k=k+1; % number of 'dangerous' N-1 contingecies
               end
               % Beautiful printing
               msg = sprintf('\tProcessed %d/%d. Number of dangerous N-1 contingecies is %d : %d', i, Sz.r(obj.E),k);
               fprintf([reverseStr, msg]);
               reverseStr = repmat(sprintf('\b'), 1, length(msg));
            else
               C1_isl (i) = 1;
               L(i,i) = -1;
            end
         end
         
         obj.L = L;
         fprintf('\n The size of N-1 islanding set is %d',sum(C1_isl(:)));
         fprintf('\n\n\tN-1 analysis was performed,%d dangerous N-1 contigencies were found, %d lines are violated \n',k,sum(lines(:)>0));
         fullpath = fileparts(mfilename('fullpath'));
         cd(fullpath);
         pathfile = strcat(pwd,'/results/N_1_analysis_dangerous_lines_',obj.case_name,'.mat');
         limits = obj.lim;
         save(pathfile,'lines','margins','L','limits','C1_isl');
      end
      
      function obj = N_2_analysis(obj,approach)
         [exists,structure] = Grid_class.file_exists('/results/N_1_analysis_dangerous_lines',obj.case_name);
         if (exists==1)
            % fprintf('\tSuccessfuly loaded N-1 contingency analysis file\n');
            if (sum(structure.lines(:)>0)>0)
               fprintf('Grid is not N-1 secure. Automatically increasing limits through lines \n');
               obj = obj.N_1_protect(structure.lines,structure.margins);
               obj = obj.N_1_analysis(); % run N-1 analysis again
               obj = obj.N_2_analysis(approach);
            else
               obj.L = structure.L;
               obj.lim = structure.limits;
               obj.C1_isl = structure.C1_isl;
               if (strcmp(sprintf(approach),sprintf('bruteforce')))
                  Grid_class.beauty_print('Start bruteforce N-2 analysis');
                  obj = obj.run_N_2_bruteforce(ones(Sz.r(obj.E),Sz.r(obj.E))-eye(Sz.r(obj.E)));
                  fprintf('\n\tRunning time for brute force algorithm is %d sec \n',obj.t_brute);
                  [pathfile] = obj.create_full_path('/results/N_2_analysis_brute_force_algorithm',obj.case_name);
                  cont_brute_force_algorithm = obj.brute_cont;
                  save(pathfile,'cont_brute_force_algorithm');
               else
                  if (strcmp(sprintf(approach),sprintf('fast')))
                     obj = obj.run_N_2_fast();
                  else
                     fprintf('Undefined approach for N-2 contingecy analysis');
                  end
               end
            end
         else
            fprintf('Run N-1 contingency analysis at first\n');
            obj = obj.N_1_analysis();
            obj = obj.N_2_analysis(approach);
         end
      end
      
      function obj = run_N_2_fast(obj)
         % K.S.Turitsyn and P.A.Kaplunovich algorithm
         C1_isl = obj.C1_isl;
         C2_isl = zeros(Sz.r(obj.E),Sz.r(obj.E)); % double outage islanding contingencies
         Grid_class.beauty_print('Start fast N-2 analysis');
         tstart = tic;
         A0 = ones(Sz.r(obj.E),Sz.r(obj.E))-eye(Sz.r(obj.E)); % matrix that depics the A set (1 - pair is in set, 0 is not)
         B0 = ones(Sz.r(obj.E),Sz.r(obj.E)); % ---=--- B set
         A  = zeros(Sz.r(obj.E),Sz.r(obj.E)); % Axy matrix
         B  = zeros(Sz.r(obj.E),Sz.r(obj.E)); % Bxz matrix
         Denominator = ones(Sz.r(obj.E),Sz.r(obj.E));
         Numerator = ones(Sz.r(obj.E),Sz.r(obj.E));
         
         A0(C1_isl(:)==1,:) = 0;
         A0(:,C1_isl(:)==1) = 0;
         A0(abs(obj.f(:))<1e-8,:)=0; % if f(i)=0 then we make an element of matrix A equal to zero
         A0(:,abs(obj.f(:))<1e-8)=0;
         qq = obj.L.*obj.L';
         tr = 1e-8; % treshold
         A0(abs(qq-1)<=tr)=0; % if lines are consequtive we don't consider them (islanding case)
         C2_isl(abs(qq-1)<=tr) = 1;
         A0(abs(qq-1)<=tr)=0;
         C2_isl(abs(qq-1)<=tr)= 1;
         % to be 100% correct i need to filter out L'.*L = 1
         
         fprintf('\t Size of C2_isl is %d\n',(sum(C2_isl(:))-Sz.r(C2_isl))/2);
         
         % now for all potential pairs that are still there we compute
         % elements of matrix Axy
         Denominator = Denominator-obj.L.*obj.L';
         Numerator = Numerator+diag(1./obj.f)*obj.L*diag(obj.f);
         A(A0~=0) = Numerator(A0~=0)./Denominator(A0~=0);
         % elements of matrix Bxz
         Bp = diag(1./(obj.lim-obj.f))*obj.L*diag(obj.f);
         Bn = -diag(1./(obj.lim+obj.f))*obj.L*diag(obj.f);
         Bn = Bn - diag(diag(Bn));
         Bp = Bp - diag(diag(Bp)); B0 = B0 - diag(diag(B0));
         %B0(:,C1_isl(:)==1) = 0;
         
         k = 0;  changing=1; 
         number_islanding_contingencies = sum(obj.C1_isl(:))*Sz.r(obj.E)-sum(obj.C1_isl(:)) + sum(C2_isl(:))/2;
         kmax = 10; % maximum number of iterations before stop
         str = {};
         while (changing==1)&&(k<kmax)
            oldA = sum(A0(:)); oldB = sum(B0(:)); % remember number of pair to undertand when to stop
            fprintf('\t %d iteration: number of potential confingecies::%d, B::%d;   islanding contingencies: %d\n',k,oldA/2,oldB,number_islanding_contingencies);
            % PHASE I
            Wbuf1 = max((diag(max(Bp))*A),(diag(min(Bp))*A));
            Wbuf2 = max((diag(max(Bn))*A),(diag(min(Bn))*A));
            W = max(Wbuf1+Wbuf1',Wbuf2+Wbuf2'); 
            
            str{k+1,1} = A0; str{k+1,2} = B0;
            str{k+1,6} = W;  str{k+1,3} = A; str{k+1,4} = Bp; str{k+1,5} = Bn;
            
            A0(W<=1)=0; A(A0==0)=0;
            
            
            % PHASE II
            Wbuf1 = max(max(Bp,[],2)*max(A),min(Bp,[],2)*min(A));
            Wbuf2 = max(max(Bn,[],2)*max(A),min(Bn,[],2)*min(A));
            Wb1 = max(Bp*diag(max(A,[],2))+Wbuf1,Bp*diag(min(A,[],2))+Wbuf1);
            Wb2 = max(Bn*diag(max(A,[],2))+Wbuf2,Bn*diag(min(A,[],2))+Wbuf2); 
            W = max(Wb1,Wb2); % bounding matrix for the set B
            
            str{k+1,7} = W;
            
            B0(W<=1)=0; Bn(B0==0)=0;Bp(B0==0)=0; k=k+1;
            
            
            if (oldA == sum(A0(:)))&&(oldB == sum(B0(:)))
               changing = 0;
            end
         end
         save('str.mat','str','-v7.3');
         obj.Mask = A0;
         fprintf('Running time of pruning loop is %d',toc(tstart));
         obj = obj.run_N_2_bruteforce(A0);
         obj.filtered_size = sum(A0(:))/2;
         obj.t_fast = toc(tstart);
         obj.C2_isl = C2_isl;
         fprintf('\n\tRunning time for fast algorithm is %d sec \n',obj.t_fast);
         [pathfile] = obj.create_full_path('/results/N_2_analysis_fast_algorithm',obj.case_name);
         cont_fast_algorithm = obj.brute_cont; 
         save(pathfile,'cont_fast_algorithm','C1_isl','C2_isl');
      end
      
      function obj = run_N_2_bruteforce(obj,A0)
         k=0; reverseStr='';
         C2_isl = zeros(Sz.r(obj.E),Sz.r(obj.E));
         brute_cont = []; 
         str = sprintf('Bruteforce enumeration over %d pairs',sum(A0(:))/2);
         fprintf('\n\t%s \n',str);
         tstart = tic;
         for i=1:(Sz.r(obj.L)-1)
            if (sum(A0(i,:))>0)
               %m = find(A0(i,:)==1);
               %m = m(find(m(:)>i));
               for j=(i+1):Sz.r(obj.L)
                  if A0(i,j)~=0 % if we consider this pair
                     if abs(det(obj.L([i,j],[i,j])))<1e-10
                        C2_isl(i,j) = 1;
                        C2_isl(j,i) = 1;
                     else
                        xq = (obj.L([i,j],[i,j])\obj.f([i,j]));
                        f_new = obj.f(:) - obj.L(:,[i,j])*xq;
                        f_new(i)=0; f_new(j)=0;
                        if sum(obj.lim(:)-abs(f_new(:))<1e-10)>0
                           k=k+1;
                           brute_cont(k,1:3)=[i,j,sum(obj.lim<abs(f_new))];
                        end
                     end
                     % Beautiful printing
                     completed = (Sz.r(obj.L)*(i-1)+1-(i+1)*i/2+j-i)/((Sz.r(obj.L)-1)*Sz.r(obj.L)/2);
                     msg = sprintf('\tProcessed %0.0f percent. Number of contingencies %d; fake %d',100*completed,k,sum(C2_isl(:))/2);
                     fprintf([reverseStr, msg]);
                     reverseStr = repmat(sprintf('\b'), 1, length(msg));
                  end
               end
            end
         end
         msg = sprintf('\tProcessed %0.0f percent. Number of contingencies %d; fake %d',100,k,sum(C2_isl(:))/2);
         fprintf([reverseStr, msg]);
         obj.C2_isl = C2_isl;
         obj.brute_cont = brute_cont;
         obj.t_brute = toc(tstart);
      end
      
      function obj = N_1_protect(obj,lines,margins) % increase limits on lines that caused N-1 contingency
         Grid_class.beauty_print('N-1 protect safe');
         mm = max(margins(lines(:)==0)); % maximum margin of safe lines
         obj.lim(lines(:)>0)=margins(lines(:)>0).*obj.lim(lines(:)>0)/mm; %increase 
      end
      
      function obj = N_1_analysis(obj)
         Grid_class.beauty_print('   Start N-1 analysis  ');
         ptdf = makePTDF(obj.rnc);
         L = makeLODF(obj.rnc.branch,ptdf);
         reverseStr=''; 
         k=0; tr = 1e-10; 
         lines = zeros(Sz.r(obj.E),1);
         margins = zeros(Sz.r(obj.E),1); % vector of maximum violations
         C1_isl = zeros(Sz.r(obj.E),1); % vector of radial contingencies
         L(isnan(L)) = -1;
         L(isinf(L)) = -1;
         C1_isl(all(L==-1,1)) = 1;
         
         for i=1:Sz.r(obj.E)
             L(i,i) = -1;
             flows = obj.f + L(:,i) * obj.f(i);
             flows(i)=0;
             qq = abs(flows(:))./obj.lim(:);
             number_of_violations = sum(abs(flows(:))>obj.lim(:));
             margins = max(margins(:),qq(:)); % Update maximum violations
             
             lines(obj.lim(:)-abs(flows(:))<tr)=lines(obj.lim(:)-abs(flows(:))<tr)+1; % mark dangerous lines
             if (number_of_violations>0)
                 k=k+1; % number of 'dangerous' N-1 contingecies
             end
             % Beautiful printing
             msg = sprintf('\tProcessed %d/%d. Number of dangerous N-1 contingecies is %d : %d', i, Sz.r(obj.E),k);
             fprintf([reverseStr, msg]);
             reverseStr = repmat(sprintf('\b'), 1, length(msg));
         end
         
         obj.L = L;
         obj.C1_isl = C1_isl;
         
         fprintf('\n The size of N-1 islanding set is %d',sum(obj.C1_isl(:)));
         fprintf('\n\n\tN-1 analysis was performed,%d dangerous N-1 contigencies were found, %d lines are violated \n',k,sum(lines(:)>0));
         fullpath = fileparts(mfilename('fullpath'));
         cd(fullpath);
         pathfile = strcat(pwd,'/results/N_1_analysis_dangerous_lines_',obj.case_name,'.mat');
         limits = obj.lim;
         save(pathfile,'lines','margins','L','limits','C1_isl');
      end
   end
   
   methods(Static)
      function beauty_print(text)
         fprintf('\n***********************************************\n');
         fprintf('************%s************\n',text);
         fprintf('***********************************************\n');
      end
      function cm = cmp(A,B)
         RE = 1e-8;
         cm = abs(A-B) < RE*max(abs(A),abs(B));
      end
      function [inversed] = woodbury_inverse(invA,U,C)
         if ~Grid_class.cmp(det(inv(C)+U*invA*U'),0)
            inversed = invA-invA*U'*((inv(C)+U*invA*U')\U)*invA;
         else
            inversed = Inf;
         end
      end
      function [exists,structure] = file_exists(relative_name,case_name)
         pathfile = Grid_class.create_full_path(relative_name,case_name);
         if (exist(pathfile,'file')==2)
            exists=1;
            structure = load(pathfile);
         else
            exists=0;
            structure = 0;
         end
      end
      function [pathfile] = create_full_path(relative_name,case_name)
         fullpath = fileparts(mfilename('fullpath'));
         cd(fullpath);
         pathfile = strcat(pwd,relative_name,'_',case_name,'.mat');
      end
   end
end

