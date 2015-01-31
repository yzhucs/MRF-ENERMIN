function obj_value = energy_vec(values, measured_data, neighbors, lik_func, penalty_func)
% FUNCTION ENERGY_VEC implements the energy function of the discrete
% minimization problem. It's a vector version of energy where we try to use
% MATLAB vector operations as much as possible, in replace of slow loops. 
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
%
% Author: Yao Zhu (yzhucs@gmail.com)

if ~isempty(lik_func)
    % compute the unary terms from likelihood function.
    lik = lik_func(values, measured_data);
    obj_value = sum(lik);
else
    obj_value = 0;
end

% here start_idx and end_idx are vectors.
start_idx = neighbors(:,1);
end_idx = neighbors(:,2);
% +1 for MATLAB indexing.
start_idx = start_idx+1;
end_idx = end_idx+1;
obj_value = obj_value + sum(penalty_func(values(start_idx),values(end_idx)));