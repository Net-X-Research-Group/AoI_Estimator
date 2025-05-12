% This function simulates the Max-Weight (MW) scheduling policy with estimation without feedback in a network
% with M streams and N maximum scheduled sources per slot. It calculates the
% actual throughput and optimal cost of the scheduling policy over a number of iterations.

function [actual_throughput,optimal_cost] = Nofeedback_MW_E_simulation(K,T,M,h1,z1,Delay,num_iterations,set_Policy,N)

global probabilityS probabilityD cost alpha Period prob_index Uniform 

%% Simulation Setup
p=zeros(M,1);
p=[p+probabilityS;0];

pD=zeros(M,1);
pD=[pD+probabilityD;0]; %Making sure that "probability" has the right size

A =zeros(M,1);
A =alpha; %Making sure that "alpha" has the right size

C =zeros(M,1);
C =[C+cost;0]; %Making sure that "cost" has the right size

Arr =zeros(M,1);
Arr = 1./Period;  %Making sure that "arrival" has the right size

Per = Period;

f = PMFgeneration;

%% Simulation
real_q=zeros(M,num_iterations);
AoI=zeros(M,num_iterations);

for iteration=1:1:num_iterations

    delivered=zeros(M+1,1);  
    c=zeros(M+1,1);  
    h_history=h1*ones(M,K);
    tauD_history = zeros(M,K);
    D_history = zeros(M,K);
    D_history(:,1)=ones(M,1);
    z=z1*ones(M,K);
    vector_successful_tx=zeros(M,1);
    t_star = zeros(M, 1);
    t_star_1 = zeros(M, 1);    
    bar_t_star = zeros(M, 1);
    bar_t_star_1 = zeros(M, 1);    
    interarrivalfactor = ones(M,1);
    g_phi1_phi = cell(M,K);
    servenode = zeros(N,K);
    for i = 1:1:M
        g_phi1_phi{i,1} =  1;
    end
    lambda=Arr.*ones(M,K);
    lambda(:,1) = ones(M,1);
    for i = 1:1:M
        lambda(i,2:K)=calculate_arrival_probabilities(i,K);
    end
    lambdaj=lambda;
    arrival_factor = zeros(M,1);
    tauD_hat_lastslot  = zeros(M,1);
    realvector_successful_tx = zeros(M,1);
    d = zeros(M,K);
    transmissiontime = zeros(M,K);
    tauD_hat_history=zeros(M,K);
    t_star_D =zeros(M,K); 
    success_factor = zeros(N,K);
    lastdeliveredpackets = zeros(M,1);
    for slot_index=1:1:K-1      
        % Estimation of real time h
        tauD_hat = zeros(M,1);
        tauA_hat = zeros(M,1);
        for node_index=1:1:M
            % Packets Arrival
                if prob_index == 1
                    if rand<Arr(node_index)
                        z(node_index,slot_index)=0;
                        vector_successful_tx(node_index)=0;% new arrival
                        realvector_successful_tx(node_index)=0;
                    elseif slot_index>1
                        z(node_index,slot_index)=z(node_index,slot_index-1)+1;
                    end
                else
                    if arrival_factor(node_index) + interarrivalfactor(node_index) == slot_index
                        z(node_index,slot_index)=0;
                        arrival_factor(node_index) = slot_index;
                        if prob_index == 2
                            interarrivalfactor(node_index) = Per(node_index);
                        elseif prob_index == 3
                            interarrivalfactor(node_index) = randi([Uniform(node_index,1), Uniform(node_index,2)]);
                        end
                        vector_successful_tx(node_index)=0;% new arrival
                        realvector_successful_tx(node_index)=0;
                    elseif slot_index>1
                        z(node_index,slot_index)=z(node_index,slot_index-1)+1;
                    end
                end

        end
%% MMSE Estimator
        for i =1:1:M
            for tau = t_star_1(i)+1:1:t_star(i)
                if tau>1
                    g_phi1_phi{i, tau} = lambdaj(i, tau);
                    pfactor=1-lambdaj(i, tau);
                    if tau > 1
                        for tauj = tau-1:-1:bar_t_star_1(i)+1 
                            g_phi1_phi{i, tau}(end+1) = pfactor*lambdaj(i, tauj);
                            pfactor = pfactor*(1-lambdaj(i, tauj));
                        end
                    end
                    g_phi1_phi{i, tau} = g_phi1_phi{i, tau}/sum(g_phi1_phi{i, tau});
                end
            end
            for tau = t_star(i)+1:1:slot_index
                if tau>1
                    if prob_index ~= 1
                        if bar_t_star(i)>1
                            lambda(i,tau) = 0;
                            for j = t_star_1(i):1:tau-1
                                if tau-j<=length(g_phi1_phi{i,tau-1})
                                    lambda(i,tau) = lambda(i,tau) + g_phi1_phi{i,tau-1}(tau-j)*f(i,tau-j); 
                                end
                            end
                        end
                    end
                    lambda(i,t_star(i)+1:bar_t_star(i))=0;
                    lambdaj(i,t_star(i)+1:bar_t_star(i))=0;
                    lambda(i,t_star_1(i)+1:bar_t_star_1(i))=0;
                    lambdaj(i,t_star_1(i)+1:bar_t_star_1(i))=0;
                    if t_star(i)>bar_t_star_1(i)
                        lambdaj(i,bar_t_star_1(i)+1:t_star(i))=lambda(i,bar_t_star_1(i)+1:t_star(i))/(1-prod(1-lambda(i,bar_t_star_1(i)+1:t_star(i))));
                    end
                    g_phi1_phi{i, tau} = lambdaj(i, tau);
                    pfactor=1-lambdaj(i, tau);
                    if tau > 1
                        for tauj = tau-1:-1:bar_t_star_1(i)+1 
                            g_phi1_phi{i, tau}(end+1) = pfactor*lambdaj(i, tauj);
                            pfactor = pfactor*(1-lambdaj(i, tauj));
                        end
                    end
                    g_phi1_phi{i, tau} = g_phi1_phi{i, tau}/sum(g_phi1_phi{i, tau});
                end
            end
        end
        % Estimate tauA
        for i=1:1:M
            for tau = slot_index:-1:bar_t_star_1(i)+1
                tauA_hat(i) = tauA_hat(i) + tau*g_phi1_phi{i, slot_index}(slot_index+1-tau);
            end
        end
        tauD_history(:,slot_index) = slot_index-h_history(:,slot_index);
        tauD_hat_new = zeros(M,1);
        
        % Estimate tauD[Di(t)]
        for i=1:1:M
                for tau = t_star(i):-1:bar_t_star_1(i)+1
                    tauD_hat_new(i) = tauD_hat_new(i) + tau*g_phi1_phi{i, t_star(i)}(t_star(i)+1-tau);
                    tauD_hat_history(i,d(i,slot_index)) = tauD_hat_new(i);
                end
        end
        % Estimate tauD(t+theta_i)
        for i=1:1:M
                if c(i)==1
                    tauD_hat(i)  = tauD_hat_lastslot(i) *(1-pD(i))+tauD_hat_new(i) *pD(i);
                    c(i)=0;
                else
                    tauD_hat(i)  = tauD_hat_lastslot(i);
                end
        end

        tauD_hat_lastslot = tauD_hat;
        h_hat = slot_index - tauD_hat;
        z_hat = slot_index - tauA_hat;
%% Max_Weight policy
        if set_Policy==1 % Max Weight                
            [~, sorted_indices] = sort(sqrt(A(1:M).*p(1:M).*pD(1:M)).*(h_hat-z_hat),'descend');
            servenode(:,slot_index) = sorted_indices(1:N);
        end    
        AoI(:,iteration)=AoI(:,iteration)+A(1:M).*h_history(:,slot_index);       
        h_history(:,slot_index+1)=h_history(:,slot_index)+1;
        d(:,slot_index+1)=d(:,slot_index);
        
        for i=1:1:N
            % receive process at the destination
                success_factor(i,slot_index) = rand;
                if slot_index > Delay
                    node_index = servenode(i,slot_index-Delay);
                    if (node_index~=(M+1)) && (success_factor(i,slot_index-Delay)<p(node_index)*pD(node_index)) && (d(node_index,slot_index-Delay+1)>lastdeliveredpackets(node_index))
                        delivered(node_index)=delivered(node_index)+1;     
                        h_history(node_index,slot_index+1)=z(node_index,slot_index-Delay)+Delay;
                        D_history(node_index,slot_index+1)=1;
                        lastdeliveredpackets(node_index)=d(node_index,slot_index-Delay+1);
                    end
                end
%% observation process at the BS
                node_index = servenode(i,slot_index);

                if (node_index~=(M+1)) && (success_factor(i,slot_index)<p(node_index)) && (vector_successful_tx(node_index)==0)
                    vector_successful_tx(node_index)=1;
                    d(node_index,slot_index+1)=d(node_index,slot_index)+1;
                    transmissiontime(node_index,d(node_index,slot_index+1))=1;
                    c(node_index)= 1;
                    t_star_1(node_index) = t_star(node_index); 
                    t_star(node_index) = slot_index;
                    t_star_D(node_index,d(node_index,slot_index+1))=slot_index;
                    bar_t_star_1(node_index) = bar_t_star(node_index);
                    bar_t_star(node_index)=t_star(node_index);
    
                elseif (node_index~=(M+1)) && (success_factor(slot_index)<p(node_index))&& (vector_successful_tx(node_index)==1)
                    bar_t_star(node_index) = slot_index;
                    transmissiontime(node_index,d(node_index,slot_index+1))=transmissiontime(node_index,d(node_index,slot_index+1))+1;
                    c(node_index)= 1;
                end

        end
    end
    real_q(:,iteration)=delivered(1:M)/K;

end
AoI_average_GRE=mean(AoI,2); 
optimal_cost_clients=T*AoI_average_GRE(:)/K;
optimal_cost=sum(optimal_cost_clients);
actual_throughput=mean(real_q,2);
end
