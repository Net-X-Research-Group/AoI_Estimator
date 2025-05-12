% This script corresponds to Figures 5, 6 in the paper.
% The main function is to compare the AoI of different schedulers varying channel reliability.
% Period: July 2024

clear
clc

global probabilityS probabilityD cost alpha arrival optQ optMU Period prob_index Uniform 

tic;

%% Model Setup
M = 8; % Number of streams in the network
N = 2; % Number of maximum scheduled sources in each slot

averagePeriod_3 = 1; % Average period/3
vector_Channelreliablity = 0.2:0.05:1; % Range of channel reliabilities
ChannelreliablityS = 0.5; % Source to BS channel reliability
ChannelreliablityD = 0.8; % BS to destination channel reliability
Delay = 5; % Transmission delay
num_setups = length(vector_Channelreliablity); % Number of different reliability setups
num_iterations = 5; % Number of iterations 
T = 1; % Frame length
h1 = 1; % Initial age
z1 = 0; % Initial system time
prob_index = 3; % Type of arrival process: 1. Bernoulli 2. Periodic 3. Uniform interarrival time

%% Parameters Setup
optimal_cost_Rand = zeros(1, num_setups);
optimal_cost_Max_Weight = zeros(1, num_setups);
optimal_cost_LowerBound_stoch = zeros(1, num_setups);
optimal_cost_Max_Weight_estimation = zeros(1, num_setups);
optimal_cost_Max_Weight_Ideal = zeros(1, num_setups);
optimal_cost_Max_Weight_estimation_nofeedback = zeros(1, num_setups);
%% Simulation
for count = 1:num_setups 
    disp(count) % Display current setup count
    K = 100000; % Some constant parameter for the simulation
    
    % Period vector for the current setup
    Period = averagePeriod_3*3 * ones(M, 1); 
    cost = zeros(M, 1);
    alpha = [4; 3; 2; 1; 5; 4; 1; 2]; % Priorities for each stream
    arrival = ones(M, 1); % Arrival rates
    probabilityD = ones(M, 1) * vector_Channelreliablity(count); % BS to destination reliability
    probabilityS = ones(M, 1) * ChannelreliablityD; % Source to BS reliability
    Uniform = repmat([2*averagePeriod_3 4*averagePeriod_3], M, 1);% Gaussian arrival times

    % Ensuring all parameters have the correct dimensions
    p = probabilityS; 
    pD = probabilityD; 
    A = alpha; 
    Arr = arrival; 
    
    % Finding optimal Q for Lower Bound
    Lower_Bound_find_q(M, N);             
    optimal_cost_LowerBound_stoch(count) = sum(alpha .* (1 + 1 ./ optQ + 2 * Delay)) / (2 * M);
    disp('Lower Bound for Stochastic Arrivals')         
    
    % Finding optimal MU for Randomized Policy
    randomized_find_optMU(M, N); 
        
    % Calculating optimal cost for Randomized Policy based on arrival process type
    if prob_index == 1
        optimal_cost_Rand(count) = sum(A .* (1 ./ Arr - 1 + 1 ./ optMU ./ p ./ pD + Delay)) / M; 
    elseif prob_index == 2
        optimal_cost_Rand(count) = sum(A .* (1 / 2 / Arr - 1 + 1 ./ optMU ./ p ./ pD + Delay)) / M; 
    elseif prob_index == 3
        optimal_cost_Rand(count) = sum(A .* (29 * averagePeriod_3 / 18 - 1 + 1 ./ optMU ./ p ./ pD + Delay)) / M;
    end

    % Simulating Max-Weight Policy with Stale Information
    set_Policy = 1; % LINEAR MaxWeight    
    [aux, optimal_cost_Max_Weight(count)] = Stale_MW_simulation(K, T, M, h1, z1, Delay, num_iterations, set_Policy, N);    
    optimal_cost_Max_Weight(count) = optimal_cost_Max_Weight(count) / (M * T);
    disp('MW with Stale Side Information')

    % Simulating Ideal Max-Weight Policy
    set_Policy = 1; % LINEAR MaxWeight with Estimation   
    [aux, optimal_cost_Max_Weight_Ideal(count)] = MW_F_simulation(K, T, M, h1, z1, Delay, num_iterations, set_Policy, N);    
    optimal_cost_Max_Weight_Ideal(count) = optimal_cost_Max_Weight_Ideal(count) / (M * T);
    disp('Ideal MaxWeight')

    % Simulating Max-Weight Policy with Estimation HoLFT Delay
    set_Policy = 1; % LINEAR MaxWeight with Estimation   
    [aux, optimal_cost_Max_Weight_estimation(count)] = MW_E_simulation(K, T, M, h1, z1, Delay, num_iterations, set_Policy, N);    
    optimal_cost_Max_Weight_estimation(count) = optimal_cost_Max_Weight_estimation(count) / (M * T);
    disp('Linear Max-Weight policy with Optimal Buffer with Estimation HoLFT')

    % Simulating Max-Weight Policy with Estimation No Feedback
    set_Policy = 1; % LINEAR MaxWeight with Estimation   
    [aux, optimal_cost_Max_Weight_estimation_nofeedback(count)] = Nofeedback_MW_E_simulation(K, T, M, h1, z1, Delay, num_iterations, set_Policy, N);    
    optimal_cost_Max_Weight_estimation_nofeedback(count) = optimal_cost_Max_Weight_estimation_nofeedback(count) / (M * T);
    disp('Linear Max-Weight policy with Optimal Buffer with Estimation HoLFT')

    toc
end

% Plotting results
figure(3)
hold on
plot(vector_Channelreliablity, optimal_cost_Max_Weight, 'mx--', 'LineWidth', 3, 'MarkerSize', 10);
plot(vector_Channelreliablity, optimal_cost_Max_Weight_Ideal, 'r*-', 'LineWidth', 3, 'MarkerSize', 10);
plot(vector_Channelreliablity, optimal_cost_Max_Weight_estimation, 'cv-', 'LineWidth', 3, 'MarkerSize', 10);
plot(vector_Channelreliablity, optimal_cost_Max_Weight_estimation_nofeedback, 'gdiamond--', 'LineWidth', 3, 'MarkerSize', 10);
plot(vector_Channelreliablity, optimal_cost_Rand, 'bo-', 'LineWidth', 3, 'MarkerSize', max(1, 10));
plot(vector_Channelreliablity, optimal_cost_LowerBound_stoch, 'k-', 'LineWidth', 3, 'MarkerSize', 2);
legend({'Randomized', 'Lower Bound', 'MW-E', 'MW-E without Feedback', 'MW-F', 'MW-S'}, 'Location', 'Northeast');
ylabel('Expected Weighted Sum AoI')
xlabel('Channel Reliability p^D_i')
xlim([0.2 1])
ylim([0 100])
hold off
toc;
