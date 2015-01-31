% This is a templated script to conduct Markove random field (MRF) energy
% minimization using the GraphCut approach. This script is templated in the
% sense that you just need to implement and provide the following variables and
% functions and then plug them into this template code.
% neighbors - the array defining the graph structure of MRF. It is an
%             edge_num*n matrix (where n>=2), where edge_num is the
%             number of edges. its 1st column stores the index of the
%             starting node, 2nd column stores the index of the ending node.
%             The extra columns may store neighboring information specific
%             to an application, .e.g., whether an edge is a wrap around edge.
% lik_func - the function handle to implement a likelihood function specific
%            to an application. This handle could be empty if we don't care
%            the unary terms in defining the energy function.
% penalty_func - the function handle to implement a pairwise penalty function
%                specific to an application.
% local_operator - the function handle to implement a local search operator
%                  specific to an application. 
%
% Author: Yao Zhu (yzhucs@gmail.com)

% load your measurement data, let the variable be measured_data.
% If you don't have measurement data, set measured_data = [];

% construct the graph structure of your MRF, let the variable be neighbors,
% which follow the format describe above.

% Implement and define your likelihood function.
% If you don't care the likelihood function, set lik_func=[];
% If you care the likelihood function, set as the following
lik_func = @(values, measured_data)your_likelihood(values, measured_data, list_of_application_parameters);

% Implement and define your pairwise penalty function.
penalty_func = @(val_start, val_end)your_likelihood(val_start, val_end, list_of_application_parameters);

% Implement and define your local search operator.
local_operator = @(values)your_local_operator(values, list_of_application_parameters);

%%%%%% solve a sequence of s-t min-cut problems %%%%%%

% maximum number of binary minimizations.
max_iter = 10;
% the vector saving the evaluation of energy function on the values of the
% discrete variables.
energy_value = zeros(max_iter+1,1);
% evaluate the energy of the initial value.
energy_value(1) = energy_vec(values, measured_data, neighbors, lik_func, penalty_func);
disp('initial energy value is: ');
disp(energy_value(1));

disp('start the sequence of binary minimization...');
for i = 1:max_iter
    disp(['solving the ' num2str(i) '-th s-t min-cut problem...']);
    % construct s-t graph.
    disp('      construct the s-t graph');
    [t_weights e_weights] = construct_st_graph_vec(values, measured_data, neighbors, ...
                                                   lik_func, penalty_func, local_operator);
    % form the adjacency matrix representing the undirected s-t graph.
    ei = neighbors(:,1);
    ej = neighbors(:,2);
    % let the node index starts from 1.
    ei = ei+1;
    ej = ej+1;
    % form an undirected graph.
    A = sparse([ei;ej],[ej;ei],[e_weights'; e_weights']);
    % form the whole graph. we form GA as a block symmetric matrix.
    GA = [0 t_weights(1:2:end) 0;
         t_weights(1:2:end)' A t_weights(2:2:end)'
         0 t_weights(2:2:end) 0];
    s = 1;
    t = size(GA,1);

    % solve the s-t min-cut problem represented by (GA, s, t). The result
    % is a cut indicator vector. Let the variable be cut. For this purpose,
    % you need a s-t min-cut solver.
    
    % You may use the max-flow solver max_flow provide by MatlabBGL
    % https://www.cs.purdue.edu/homes/dgleich/packages/matlab_bgl/
    % OR you may use the PIRMCut solver that can handle floating-point valued
    % problem instances. See http://arxiv.org/abs/1501.03105
    % CAUTION: If you use a min-cut solver that can only handle interger
    % values, you need to scale and rounding GA properly to get meaningful
    % results.
    
    % apply the local operator according to the cut indicator vector.
    values = apply_local_operator(values, cut, local_operator);
    % evalute the energy on the updated values.
    energy_value(i+1) = energy_vec(values, measured_data, neighbors, lik_func, penalty_func);
    disp('current energy value is: ');
    disp(energy_value(i+1));
    % check whether we have seen a local optimum.
    if energy_value(i+1) >= energy_value(i)
        iter = i;
        break;
    end
end
disp('Done with the sequence of binary minimization.');
energy_value(iter+2:end)=[];