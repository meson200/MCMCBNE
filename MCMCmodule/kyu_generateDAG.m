function dag_random = kyu_generateDAG(dag_caus,density,maxfi)

% randomly generates a DAG from a given causality constraint (dag_caus) 
% the DAG has to pass two tests:
% 1) no loop
% 2) fan in to every node <= maxfi
% an input variable 'density' is the number of links you desire in the
% randomly sampled graph
<<<<<<< HEAD

function is_dag = isdag(R,dim,maxfi)
=======

dim = size(dag_caus);
[ix,iy] = find(dag_caus);
isdag = false;

n_possible_links = numel(find(dag_caus(:)));
nlinks = round(density*n_possible_links)+1;

while ~isdag

    R = zeros(dim); 
    indx_s = randsample(length(ix),nlinks,false);
    %R([ix(indx_s) iy(indx_s)]) = 1;
    for i = 1:nlinks
       R(ix(indx_s(i)),iy(indx_s(i))) = 1;
    end
    % check cyclicity
>>>>>>> 2b47fa88270021c22d16c9ae0e57f768d24431c7
    RR = reachability_graph(R);
    % check maximum fan in
    tmp = sum(R,1) - maxfi*ones(1,dim(1));
    ind_c = find(max(0,tmp)); % detect a child which would exceed max num of parents
<<<<<<< HEAD
    is_dag = ~any(diag(RR)==1) && isempty(ind_c);
end

dim = size(dag_caus);
R = zeros(dim);
[ix,iy] = find(dag_caus);
n_possible_links = numel(find(dag_caus(:)));
nlinks = round(density*n_possible_links)+1;
indx_s = randsample(length(ix),nlinks,false);
% randomly generate 1 DAG
for i = 1:nlinks
       R(ix(indx_s(i)),iy(indx_s(i))) = 1;
end

if isdag(R,dim,maxfi)
   dag_random = R; 
else
   while ~isdag(R,dim,maxfi)
        % remove edges until it becomes DAG
        [non_x,non_y] = find(R);
        which = randsample(length(non_x),1,false);
        R(non_x(which),non_y(which)) = 0;
   end
   dag_random = R;
end

end
=======
    isdag = ~any(diag(RR)==1) && isempty(ind_c);
   
end
dag_random = R;
>>>>>>> 2b47fa88270021c22d16c9ae0e57f768d24431c7
