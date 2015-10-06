function [Gs, op, nodes, A] = mk_nbrs_of_digraph_kyu(G0,A,maxfi,dag_caus)
% MY_MK_NBRS_OF_DIGRAPH Make all digraphs that differ from G0 by a single edge deletion, addition or reversal, subject to acyclicity
% [Gs, op, nodes, A] = my_mk_nbrs_of_digraph(G0,<A>)
%
% G0 is an adj matrix s.t. G0(i,j)=1 iff i->j in graph
% A is the ancestor matrix for G0  (opt, creates if necessary)
%
% Gs(:,:,i) is the i'th neighbor
% op{i} = 'add', 'del', or 'rev' is the operation used to create the i'th neighbor. 
% nodes(i,1:2) are the head and tail of the operated-on arc.
% Modified from mk_nbrs_of_digraph by Sonia Leach
%
% Modified by Sonia Leach Feb 02

if isempty(A), A = reachability_graph(G0'); end

n = length(G0);
[I,J] = find(G0); % I(k), J(k) is the k'th edge
E = length(I);    % num edges present in G0


% SINGLE EDGE DELETIONS
% all deletions are valid wrt acyclity

Grep = repmat(G0(:), 1, E); % each column is a copy of G0
% edge_ndx(k) is the scalar location of the k'th edge 
edge_ndx = find(G0);

% edge_ndx = subv2ind([n n], [I J]); % equivalent
% We set (ndx(k), k) to 0 for k=1:E in Grep
ndx = subv2ind(size(Grep), [edge_ndx(:) (1:E)']);
G1 = Grep;
G1(ndx) = 0;
Gdel = reshape(G1, [n n E]);


% SINGLE EDGE REVERSALS

% SML: previously Kevin had that legal structure was if
% A(P,i)=1 for any P = { p | p in parents(j), p~=i}
% specifically he said 
%  "if any(A(ps,i)) then there is a path i -> parent of j -> j
%   so reversing i->j would create a cycle"
% Thus put in another way:
%    for each i,j if sum(G0(:,j)' * A(:,i)) > 0, reversing i->j
% is not legal.
%
% Ex. Suppose we want to check if 2->4 can be reversed in the 
% following graph: 
% G0 =                               A =
%     0     0     1     0               0     0     0     0
%     0     0     1     1               0     0     0     0
%     0     0     0     1               1     1     0     0
%     0     0     0     0               1     1     1     0
% 
% Then parents(4) = G0(:,4) = [0 1 1 0]'
% and A(:,2) = [0 0 1 1]. Thus G0(:,4)'*A(:,2) = 1 b/c 3 is
% an ancestor of 4 and a child of 2. Note that this works b/c
% matrix multiplication has the effect of ANDing the two vectors 
% and summing up the result (equiv. to the any(A(ps,i)) in kevin's code)
%
% So, we vectorize and check for all i,j pairs by looking for
% 1's in L = (G0'*A)' which has L(i,j)=1 if rev(i,j) not legal
% Note that this will give 1's where there are none in the G0
% so we do a L=max(0, G0-L) to cancel out only the existing edges that 
% aren't legal (subtracting where both are 1 and setting where
% G0=0 and A=1 back to 0).

L = max(0, G0-(G0'*A)');
[IL, JL] = find(L);  % I(k), J(k) is the k'th legal edge to rev.
[IL,JL] = filter_maxfanin(G0,IL,JL,n,maxfi); % Maximum fan-in check
EL = length(IL);

if EL>0

    % SML: First we have to DELETE THE EDGES WE ARE REVERSING
    % We can't use G1 w/ reversed edges already deleted (as
    % Kevin did) b/c the space of possible deletions are different 
    % now (some reverses aren't legal)

    Grep = repmat(G0(:), 1, EL); % each column is a copy of G0
    % edge_ndx(k) is the scalar location of the k'th edge 
    edge_ndx = subv2ind([n n], [IL JL]); 
    % We set (ndx(k), k) to 0 for k=1:E in Grep
    ndx = subv2ind(size(Grep), [edge_ndx(:) (1:EL)']);
    G1 = Grep;
    G1(ndx) = 0;

    % SML: Now we add in our REVERSED EDGES
    % rev_edge_ndx(k) is the scalar location of the k'th legal reversed edge
    rev_edge_ndx = subv2ind([n n], [JL IL]);

    % We set (rev_edge_ndx(k), k) to 1 for k=1:EL in G1
    % We have already deleted i->j in the previous step
    ndx = subv2ind(size(Grep), [rev_edge_ndx(:) (1:EL)']);
    G1(ndx) = 1;
    Grev = reshape(G1, [n n EL]);
else
    Grev = [];
end
% SINGLE EDGE ADDITIONS

% SML: previously Kevin had that any addition was legal if A(i,j)=0
% however, you can not add i->j  if j is a descendent of i.
% Thus, we create all possible additions in Gbar and then
% subtract the descendants of each edge as possible parents
% This means the potential parents of i (i.e. Gbar(:,i))
% can not also be descendants if i i.e. (A(:,i)) which is accomplished
% by subtracting (Gbar-A == 1 iff Gbar=1 & A=0)

Gbar = ~G0;  % Gbar(i,j)=1 iff there is no i->j edge in G0
Gbar = setdiag(Gbar, 0); % turn off self loops

GbarL = Gbar-A;
[IbarL,JbarL] = find(GbarL);  % I(k), J(k) is the k'th legal edge to add
[IbarL,JbarL] = filter_maxfanin(G0,IbarL,JbarL,n,maxfi); % Maximum fan-in check

EbarL = length(IbarL);

if EbarL>0

    %bar_edge_ndx = find(GbarL);
    bar_edge_ndx = subv2ind([n n],[IbarL, JbarL]);
    %isequal(bar_edge_ndx,bar_edge_ndx2)

    Grep = repmat(G0(:), 1, EbarL); % each column is a copy of G0
    ndx = subv2ind(size(Grep), [bar_edge_ndx(:) (1:EbarL)']);
    Grep(ndx) = 1;
    Gadd = reshape(Grep, [n n EbarL]);

else
    Gadd = [];
end

Gs = cat(3, Gdel, Grev, Gadd);

nodes = [I J;
     IL JL;
   IbarL JbarL];

op = cell(1, E+EL+EbarL);
op(1:E) = {'del'};
op(E+(1:EL)) = {'rev'};
op((E+EL+1):end) = {'add'}; 

% check if the proposed move is not causal

[Gs,op,nodes] = filter_causality(Gs,op,nodes,dag_caus);

end


% max fan-in restriction 
% added by Kyu on Aug 18 2012

function [ILnew JLnew] = filter_maxfanin(G0,IL,JL,n,maxfi,ops)

ILnew = []; JLnew = [];    
for i = 1:length(IL)    
    Gmoved = G0;
    Gmoved(IL(i),JL(i)) = 1;
    tmp = sum(Gmoved,1) - maxfi*ones(1,n);
    ind_c = find(max(0,tmp)); % detect a child which would exceed max num of parents
    if isempty(ind_c) % all edges satisfy the limit
        ILnew = [ILnew IL(i)];
        JLnew = [JLnew JL(i)];
    end
    
end
ILnew = ILnew';
JLnew = JLnew';

% deleted_edge = length(IL) - length(ILnew)

end

function [Gs_new,op_new,nodes_new] = filter_causality(Gs_old,op_old,nodes_old,dag_caus)
Gs_new = Gs_old;
op_new = op_old;
nodes_new = nodes_old;

if ~isempty(dag_caus)
    sz = size(dag_caus);
    ndag = length(op_old);
    toremove = zeros(ndag,1);
    for u = 1:length(op_old)
        for i = 1:sz(1)
            for j = 1:sz(2)
                if dag_caus(i,j) == 0 && Gs_old(i,j,u) == 1
                    toremove(u) = 1;
                end
            end
        end
    end
    Gs_new(:,:,find(toremove)) = [];
    op_new(find(toremove)) = [];
    nodes_new(find(toremove),:) = [];
end
    
end



