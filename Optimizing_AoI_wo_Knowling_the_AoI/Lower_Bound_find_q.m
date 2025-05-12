% This function calculates the optimal Q values (optQ) for the lower bound 
% of the scheduling problem in a network with M streams and N maximum scheduled sources per slot.

function Lower_Bound_find_q(M, N)

global probabilityS probabilityD alpha arrival optQ

iterations = 200; % Number of iterations for the binary search

% Initialize probability, alpha, and arrival arrays
p = zeros(M, 1);
probability = p + probabilityS .* probabilityD; % Ensure "probability" has the correct size

A = zeros(M, 1);
A = A + alpha; % Ensure "alpha" has the correct size

Arr = zeros(M, 1);
Arr = Arr + arrival; % Ensure "arrival" has the correct size
Arr = min(Arr, probability); % Limit arrival to the probability

% Initialize optQ
optQ = zeros(M, 1);

% If the sum of Arr/probability is less than or equal to N, assign directly
if sum(Arr ./ probability) <= N
    for node = 1:M
        optQ(node) = Arr(node);
    end
else
    % Calculate Gamma_star and gamma according to the new algorithm description
    Gamma_star = (sum(sqrt(A ./ probability)))^2 / (2 * M * N^2);
    gamma = sqrt(A .* probability ./ (2 * M * Gamma_star));
    
    % Calculate GammaI and initialize max and min Gamma
    GammaI = A .* probability ./ (2 * M * (Arr.^2));
    maxGamma = max(GammaI);
    minGamma = 0;

    % Use binary search to iteratively update Gamma_aux
    for count = 1:iterations
        Gamma_aux = (maxGamma + minGamma) / 2;
        Q_aux = min(Arr, sqrt(A .* probability ./ (2 * M * Gamma_aux)));
        S = sum(Q_aux ./ probability);
        if S <= N
            maxGamma = Gamma_aux;
        else
            minGamma = Gamma_aux;
        end
    end

    % Finally, calculate optQ based on the optimized Gamma
    optQ = min(Arr, sqrt(A .* probability * N^2 ./ (2 * M * maxGamma)));
end
end
