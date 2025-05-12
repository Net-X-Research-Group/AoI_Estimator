% This function calculates the arrival probabilities for a given source based on different probability distributions (Bernoulli, Periodic, Uniform).
% and Algorithm 3 in the paper


function arrival_probabilities = calculate_arrival_probabilities(source_index, K)
    global Period arrival Uniform prob_index 
    Arr = arrival;
    per = Period;
   
    % Determine the arrival rate based on the probability index
    if prob_index == 1
        arr = Arr(source_index); % Bernoulli arrival
    elseif prob_index == 2 
        arr = per(source_index); % Periodic arrival
    elseif prob_index == 3 
        arr = Uniform(source_index, :); % Uniform arrival
    end

    slot = 1000;
    timelength = K - 1;
    % Generate PMF based on the probability index and arrival rate
    pmf = generate_pmf(prob_index, slot, arr);

    arrival_probabilities = zeros(1, timelength);

    % Get the length of the PMF
    pmf_length = length(pmf);

    % Calculate the arrival probability for each time slot
    for t = 1:timelength
        % Initialize the probability for the current time slot
        probability = 0;
        
        % Iterate to calculate cumulative arrival probability
        for k = 1:pmf_length
            if t - k > 0
                probability = probability + arrival_probabilities(t - k) * pmf(k);
            end
        end

        % Add the base probability for the current time slot
        if t <= pmf_length
            probability = probability + pmf(t);
        end

        arrival_probabilities(t) = probability;
    end
end

% This function generates the Probability Mass Function (PMF) based on the
% probability index and the arrival rate for different distributions.
function pmf = generate_pmf(prob_index, slot, arr)
    % Generate PMF based on the type of arrival process
    switch prob_index
        case 1 % Bernoulli arrival, Erlang distribution
            p = arr; % Arrival probability
            pmf = p * (1 - p) .^ (0:slot-1);
            pmf = pmf / sum(pmf); % Normalize
        case 2 % Periodic arrival, constant distribution
            pmf = zeros(1, slot);
            interval = arr; % Period
            if interval <= slot
                pmf(interval) = 1; % Constant arrival time
            end
        case 3 % Uniform arrival
            pmf = zeros(1, slot);
            pmf(arr(1):arr(2)) = 1 / (arr(2) - arr(1) + 1); % Uniformly distributed between arr(1) and arr(2)
        otherwise
            error('Invalid prob_index. Must be 1, 2, or 3.');
    end
end
