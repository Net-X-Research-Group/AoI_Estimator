% This script corresponds to Figures 7 in the paper.
% The main function is to compare the AoI of different schedulers varying Delay.
% Period: July 2024

clear
clc

global probabilityS probabilityD cost alpha arrival optQ optMU Period prob_index Uniform 

tic;

%% Model Setup
M = 8; % Number of streams in the network
N = 2; % Number of maximum scheduled sources in each slot

vector_Delay = 0:1:20; % Range of delays
ChannelreliablityS = 0.5; % Source to BS channel reliability
ChannelreliablityD = 0.8; % BS to destination channel reliability
Period = ones(M, 1) * 3; % Period of arrivals

num_setups = length(vector_Delay); % Number of different delay setups
num_iterations = 10; % Number of iterations for averaging results
T = 1; % Frame length
set_simplified = 1; % Simplification flag
h1 = 1; % Initial age
z1 = 0; % Initial system time, must be lower than h1
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
    disp(count)
    K = 100000; % Time Horizon
    cost = zeros(M, 1);
    alpha = [4; 3; 2; 1; 5; 4; 1; 2]; % Priorities for each stream
    arrival = ones(M, 1) ./ Period; % Arrival rates
    probabilityS = ones(M, 1) * ChannelreliablityS; % Source to BS reliability
    probabilityD = ones(M, 1) * ChannelreliablityD; % BS to destination reliability
    Uniform = repmat([2 4], M, 1); % Uniform arrival times

    % Ensuring all parameters have the correct dimensions
    p = zeros(M, 1) + probabilityS;
    pD = zeros(M, 1) + probabilityD;
    A = zeros(M, 1) + alpha;
    Arr = zeros(M, 1) + arrival;

    Delay = vector_Delay(count); % Current delay setup
    
    Lower_Bound_find_q(M, N); % Find the optimal Q for Lower Bound
    optimal_cost_LowerBound_stoch(count) = sum(alpha .* (1 + 1 ./ optQ+ 2*Delay) ) / ( 2* M);
    disp('Lower Bound for Stochastic Arrivals')         
    
    randomized_find_optMU(M, N); % Find the optimal MU for Randomized Policy

    % Calculating optimal cost for Randomized Policy based on arrival process type
    if prob_index == 1
        optimal_cost_Rand(count) = sum(A .* (1 ./ Arr - 1 + 1 ./ optMU ./ p ./ pD + Delay)) / M; 
    elseif prob_index == 2
        optimal_cost_Rand(count) = sum(A .* (1 / 2 ./ Arr - 1 + 1 ./ optMU ./ p ./ pD + Delay)) / M; 
    elseif prob_index == 3
        optimal_cost_Rand(count) = sum(A .* (29 * 3 / 18 - 1 + 1 ./ optMU ./ p ./ pD + Delay)) / M;
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
[Linear_Random_plot] = plot(vector_Delay, optimal_cost_Rand, 'bo-', 'LineWidth', 3, 'MarkerSize', 10);
[Linear_Max_Weight_plot] = plot(vector_Delay, optimal_cost_Max_Weight, 'mx-', 'LineWidth', 3, 'MarkerSize', 10);
[Linear_Max_Weight_Ideal_plot] = plot(vector_Delay, optimal_cost_Max_Weight_Ideal, 'r*-', 'LineWidth', 3, 'MarkerSize', 10);
[Linear_Max_Weight_Estimation_HoLFT_delay_plot_nofeedback] = plot(vector_Delay, optimal_cost_Max_Weight_estimation_nofeedback, 'cv-', 'LineWidth', 3, 'MarkerSize', 10);
[Linear_Max_Weight_Estimation_HoLFT_delay_plot] = plot(vector_Delay, optimal_cost_Max_Weight_estimation, 'gdiamond--', 'LineWidth', 3, 'MarkerSize', 10);
[Lower_Bound] = plot(vector_Delay, optimal_cost_LowerBound_stoch, 'k-', 'LineWidth', 3, 'MarkerSize', 10);
legend([ Lower_Bound, Linear_Max_Weight_Estimation_HoLFT_delay_plot, Linear_Max_Weight_Estimation_HoLFT_delay_plot_nofeedback, Linear_Max_Weight_Ideal_plot], ...
     'Lower Bound', 'MW with Estimation', 'MW with No-Feedback Estimation', 'MW with Full Knowledge', 'Location', 'Northeast');
ylabel('Expected Weighted Sum AoI')
xlabel('Delay')
xlim([0 20])
ylim([0 100])

% zoom in
% xZoom = [1, 3];
% yZoom = [6, 25];
% 
% rectangle('Position', [xZoom(1), min(yZoom), diff(xZoom), range(yZoom)], 'EdgeColor', 'k', 'LineWidth', 1.5);
% 
% ax = axes('Position', [0.15, 0.6, 0.25, 0.25]);
% hold(ax, 'on');
% 
% plot(ax, vector_Delay, optimal_cost_Max_Weight, 'mx-', 'LineWidth', 3, 'MarkerSize', 10);
% plot(ax, vector_Delay, optimal_cost_Max_Weight_estimation_Ideal, 'r*-', 'LineWidth', 3, 'MarkerSize', 10);
% plot(ax, vector_Delay, optimal_cost_Max_Weight_estimation, 'g-', 'LineWidth', 3, 'MarkerSize', 10);
% plot(ax, vector_Delay, optimal_cost_Max_Weight_estimation_nofeedback, 'gx--', 'LineWidth', 3, 'MarkerSize', 10);
% plot(ax, vector_Delay, optimal_cost_Rand, 'bo-', 'LineWidth', 3, 'MarkerSize', 10);
% plot(ax, vector_Delay, optimal_cost_LowerBound_stoch, 'k-', 'LineWidth', 3, 'MarkerSize', 10);
% 
% xlim(ax, xZoom);
% ylim(ax, yZoom);
% 
% hold off;
