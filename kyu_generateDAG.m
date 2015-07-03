function dag_random = kyu_generateDAG(dim,density,dag_caus)

% randomly generates a DAG with a given dimension (dim) 
% and the density of nodes (density:0-1)

dag_random = zeros(dim);
isdag=false;
%count = 1;
while ~isdag
   %count = count+1;
   %disp(count)
   R = sprand(dim(1),dim(2),density);
   for i = 1:dim(1)
       for j = 1:dim(2)
         if i==j 
             R(i,j) = 0; 
         elseif R(i,j)>0.5
             R(i,j) = 1;
         else
             R(i,j) = 0;
         end
       end
   end
   % check cyclicity
   RR = reachability_graph(R);
   
   % check causality
   noncausal = 0;
   if ~isempty(dag_caus)
       for i = 1:dim(1)
            for j = 1:dim(2)
                if dag_caus(i,j) == 0 && R(i,j) == 1
                    noncausal = 1;
                end
            end
       end
       isdag = ~any(diag(RR)==1) && ~noncausal;
   end
end
dag_random = full(R);  
