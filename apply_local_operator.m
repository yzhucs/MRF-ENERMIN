function values = apply_local_operator(values, cut, local_operator)
% FUNCTION APPLY_LOCAL_OPERATOR applies the local search operator to the
% target variable. This function will be called in solving a discrete
% optimization problem through solving a sequence of s-t min-cut problems.
% values, the current value of the target variables.
% cut, a {-1,+1} indicator variable as returned by solving an s-t min-cut
% problem. cut(i)=1 means node i is on the source side, otherwise it's on
% the sink side.
% local_operator, a function handle implementing the admissible local
% search operations.
% return - values, the updated value of target variables after applying the
% local search operation.
%
% Author: Yao Zhu (yzhucs@gmail.com)

% values and cut must be of the same size.
if sum(size(values) == size(cut)) ~= length(size(values))
    error('values and cut dimension not match in apply_local_operator.');
end

% the value assigned to node i changes only if it is on the sink side.
index = (cut <= 0);
values(index) = local_operator(values(index));