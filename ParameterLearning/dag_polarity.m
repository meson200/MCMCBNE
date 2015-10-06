function [pos,neg] = dag_polarity(bnet,engine)

% returns the directed matrix that denotes the 
% polarity of interaction
% M(i,j) = 1: higher value of i increases the chances of higher j
% M(i,j) = -1: higher value of i decrases the chances of higher j

dag = bnet.dag;
nvar = size(dag,1);
pos = zeros(nvar,nvar);
neg = zeros(nvar,nvar);
for i = 1:nvar
    for j = 1:nvar
        if dag(i,j) > 0
            % i: parent node, j: child node
            % construct P(j|i) = P(i,j)/P(i)
            % we are interested in P(j=2|i=1) and P(j=2|i=2)
            ev2 = cell(1,nvar);
            ev2{i} = 2;
            ev1 = cell(1,nvar);
            ev1{i} = 1;     
            [eng1, ~] = enter_evidence(engine, ev1);
            marg = marginal_nodes(eng1, j);
            p_low = marg.T(2);
            [eng2, ~] = enter_evidence(engine, ev2);
            marg = marginal_nodes(eng2, j);
            p_hi = marg.T(2);
            if p_hi > p_low
                pos(i,j) = 1;
            else
                neg(i,j) = 1;
            end
          
        end
    end
end
    