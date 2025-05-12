% This function simulates the Stale Max-Weight (MW) scheduling policy in a network
% with M streams and N maximum scheduled sources per slot. It calculates the
% actual throughput and optimal cost of the scheduling policy over a number of iterations.
%
% Parameters:
% K - Total time slots
% T - Frame length
% M - Number of streams in the network
% h1 - Initial age
% z1 - Initial system time
% Delay - Transmission delay
% num_iterations - Number of iterations for averaging results
% N - Number of maximum scheduled sources per slot

function [actual_throughput, optimal_cost] = Stale_MW_simulation(K, T, M, h1, z1, Delay, num_iterations, set_Policy, N)

global probabilityS probabilityD cost alpha arrival Period prob_index Uniform 

%% Simulation Setup

% Initialize probability, alpha, cost, and arrival arrays
p = zeros(M, 1);
p = [p + probabilityS; 0];

pD = zeros(M, 1);
pD = [pD + probabilityD; 0]; % Ensure "probability" has the right size

A = zeros(M, 1);
A = [A + alpha; 0]; % Ensure "alpha" has the right size

C = zeros(M, 1);
C = [C + cost; 0]; % Ensure "cost" has the right size

Arr = zeros(M, 1);
Arr = [Arr + arrival; 0]; % Ensure "arrival" has the right size

Per = Period;
optimal_cost = 0;

%% Simulation

real_q = zeros(M, num_iterations);
AoI = zeros(M, num_iterations);
node_index = 0;

for iteration = 1:num_iterations

    % Initialize variables for the simulation
    delivered = zeros(M + 1, 1);    
    h_history = h1 * ones(M, K);
    h_BS = h1 * ones(M, 1);
    z = z1 * ones(M, 1);
    vector_successful_tx = ones(M, 1);
    arrival_factor = zeros(M, 1);
    interarrivalfactor = ones(M, 1);
    servenode = zeros(N, K);

    for slot_index = 1:K
     
        if slot_index < 2 * Delay + 1
            h_BS(:) = h_BS(:) + 1;
        else 
            h_BS(:) = h_history(:, slot_index - 2 * Delay);
        end

        for node_index = 1:M
            if prob_index == 1
                if rand < Arr(node_index)
                    z(node_index) = 0;
                    vector_successful_tx(node_index) = 0; % new arrival
                else
                    z(node_index) = z(node_index) + 1;
                end
            else
                if arrival_factor(node_index) + interarrivalfactor(node_index) == slot_index
                    z(node_index) = 0;
                    arrival_factor(node_index) = slot_index;
                    if prob_index == 2
                        interarrivalfactor(node_index) = Per(node_index);
                    elseif prob_index == 3
                        interarrivalfactor(node_index) = randi([Uniform(node_index, 1), Uniform(node_index, 2)]);
                    end
                    vector_successful_tx(node_index) = 0; % new arrival
                else
                    z(node_index) = z(node_index) + 1;
                end
            end
        end

        z_hat = 1 ./ Arr(1:M, :);
        if set_Policy == 1 % The LINEAR Max Weight Policy is used                
            [~, sorted_indices] = sort(sqrt(A(1:M) .* p(1:M) .* pD(1:M)) .* (h_BS - z_hat), 'descend');
            servenode(:, slot_index) = sorted_indices(1:N);
        end         
        
        AoI(:, iteration) = AoI(:, iteration) + A(1:M) .* h_history(:, slot_index);       
        h_history(:, slot_index + 1) = h_history(:, slot_index) + 1;                
        
        for i = 1:N
            node_index = servenode(i, slot_index);
            if (node_index ~= (M + 1)) && (rand < p(node_index) * pD(node_index)) && (vector_successful_tx(node_index) == 0)
                delivered(node_index) = delivered(node_index) + 1;
                vector_successful_tx(node_index) = 1;                
                h_history(node_index, slot_index + 1) = z(node_index) + Delay;
            end  
        end
 
    end
    real_q(:, iteration) = delivered(1:M) / K;
end

% Calculate average AoI and optimal cost
AoI_average_GRE = mean(AoI, 2); % Averaging over iterations
optimal_cost_clients = T * AoI_average_GRE(:) / K;
optimal_cost = sum(optimal_cost_clients);
actual_throughput = mean(real_q, 2);

end
