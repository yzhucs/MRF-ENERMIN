function [t_weights e_weights] = ...
    construct_st_graph_vec(values, measured_data, neighbors, lik_func, ...
                           penalty_func, local_operator)
% FUNCTION CONSTRUCT_ST_GRAPH_VEC constructs the s-t graph on which we solve
% the s-t min-cut problem for binary minimization. For a given s-t cut, if
% a node x is on the souce side, it is interpreted as x=0, and otherwise
% x=1. When embed a binary minimization inside a discrete optimization, x=0
% indicates the values stay the same, and x=1 indcates the values changed
% to the value by applying the local operator. It's a vector version
% of construct_st_graph where we try to use MATLAB vector operations as
% much as possible, in replace of slow loops.
% values, the current value assigned to the target variables.
% measured_data, the measurements.
% neighbors, is an edge_num*n matrix (where n>=2), where edge_num is the
% number of edges. its 1st column stores the index of the starting node,
% 2nd column stores the index of the ending node. The extra columns may
% store neighboring information specific to an application, .e.g., whether
% an edge is a wrap around edge.
% lik_func, the function handle to implement a likelihood function specific
% to an application. This handle could be empty if we don't care the unary
% terms in defining the energy function.
% penalty_func, the function handle to implement a pairwise penalty
% function specific to an application.
% local_operator, the function handle to implement a local search
% operator specific to an application. 
% return - t_weights, the array storing the weights on terminal edges.
% return - e_weights, the array storing the weights on non-terminal edges.
%
% Author: Yao Zhu (yzhucs@gmail.com)
%
%   Reference: [1] V. Kolmogorov, R. Zabih, What Energy Functions Can be
%              Minimized via Graph Cuts?

% construct part of the graph according to neighboring system and penalty
% function.

% here start_idx and end_idx are vectors.
start_idx = neighbors(:,1);
end_idx = neighbors(:,2);
% +1 for MATLAB indexing.
start_idx = start_idx+1;
end_idx = end_idx+1;
% construction according to section 4 of reference [1].
a = penalty_func(values(start_idx),values(end_idx))';
b = penalty_func(values(start_idx),local_operator(values(end_idx)))';
c = penalty_func(local_operator(values(start_idx)),values(end_idx))';
d = penalty_func(local_operator(values(start_idx)), local_operator(values(end_idx)))';

e_weights = (c + b - a - d);

% get the total number of variables in MRF.
num_nodes = numel(values);
t_weights = zeros(1,2*num_nodes);

for neighbor_idx = 1:size(neighbors,1)
    % here start_idx and end_idx are scalars.
    start_idx = neighbors(neighbor_idx,1);
    end_idx = neighbors(neighbor_idx,2);  
    % +1 for MATLAB indexing.
    start_idx = start_idx+1;
    end_idx = end_idx+1;
    if (c(neighbor_idx) > a(neighbor_idx))
        t_weights(2*start_idx-1) = t_weights(2*start_idx-1) + c(neighbor_idx) - a(neighbor_idx);
    else
        t_weights(2*start_idx) = t_weights(2*start_idx) + a(neighbor_idx) - c(neighbor_idx);
    end
    if (d(neighbor_idx) > c(neighbor_idx))
        t_weights(2*end_idx-1) = t_weights(2*end_idx-1) + d(neighbor_idx) - c(neighbor_idx);
    else
        t_weights(2*end_idx) = t_weights(2*end_idx) + c(neighbor_idx) - d(neighbor_idx);
    end 
    if (b(neighbor_idx) > a(neighbor_idx))
        t_weights(2*end_idx-1) = t_weights(2*end_idx-1) + b(neighbor_idx) - a(neighbor_idx);
    else
        t_weights(2*end_idx) = t_weights(2*end_idx) + a(neighbor_idx) - b(neighbor_idx);
    end 
    if (d(neighbor_idx) > b(neighbor_idx))
        t_weights(2*start_idx-1) = t_weights(2*start_idx-1) + d(neighbor_idx) - b(neighbor_idx);
    else
        t_weights(2*start_idx) = t_weights(2*start_idx) + b(neighbor_idx) - d(neighbor_idx);
    end
end

% construct part of the graph according to measurements and likelihood function.
if ~isempty(lik_func)
    sink_lik = lik_func(values, measured_data);
    src_lik = lik_func(local_operator(values), measured_data);
    index = find(src_lik > sink_lik);
    t_weights(2*index-1) = t_weights(2*index-1) + src_lik(index) - sink_lik(index);
    index = find(src_lik < sink_lik);
    t_weights(2*index) = t_weights(2*index) + sink_lik(index) - src_lik(index);
end