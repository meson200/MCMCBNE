function dag_random = kyu_generateDAG(dag_caus,density,maxfi)

% randomly generates a DAG from a given causality constraint (dag_caus) 
% the DAG has to pass two tests:
% 1) no loop
% 2) fan in to every node <= maxfi
% an input variable 'density' is the number of links you desire in the
% randomly sampled graph

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
    RR = reachability_graph(R);
    % check maximum fan in
    tmp = sum(R,1) - maxfi*ones(1,dim(1));
    ind_c = find(max(0,tmp)); % detect a child which would exceed max num of parents
    isdag = ~any(diag(RR)==1) && isempty(ind_c);
   
end
dag_random = R;
