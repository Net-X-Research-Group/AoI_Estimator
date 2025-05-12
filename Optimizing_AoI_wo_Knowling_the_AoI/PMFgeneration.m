% This function generates a vector of arrival probabilities at each slot without any condition in an iterative mannerfunction [f] = PMFgeneration

function [f] = PMFgeneration

global Period prob_index Uniform

% Initialize PMF based on the type of arrival process
if prob_index == 1
    f = 0; % Bernoulli arrival (placeholder, no implementation provided)
    
elseif prob_index == 2
    % Periodic arrival process
    T = Period(:, 1); 
    num_bins = 1000; % Number of bins for the PMF
    x = 1:num_bins;
    f = zeros(length(T), num_bins);
    
    % Generate PMF for each node
    for node_index = 1:length(T)
        for n = 1:num_bins
            if mod(n, T(node_index)) == 0
                f(node_index, n) = 1;
            end
        end
    end
    
elseif prob_index == 3
    % Uniform arrival process
    startvalue = Uniform(:, 1); 
    endvalue = Uniform(:, 2);  
    num_bins = 1000; % Number of bins for the PMF
    x = 1:num_bins;
    pdf_values = zeros(1, num_bins);
    f = zeros(length(endvalue), num_bins);
    
    % Generate PMF for each node
    for node_index = 1:length(endvalue)
        pdf_values(startvalue(node_index):endvalue(node_index)) = ...
            ones(endvalue(node_index) - startvalue(node_index) + 1, 1) / ...
            (endvalue(node_index) - startvalue(node_index) + 1);
        f(node_index, :) = pdf_values / sum(pdf_values);
    end
end

end