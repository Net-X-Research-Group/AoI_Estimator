function [actual_throughput,optimal_cost]= MW_F_simulation(K,T,M,h1,z1,Delay,num_iterations,set_Policy,N)

global probabilityS probabilityD cost alpha arrival optMU Period prob_index Uniform Gaussin

%% Simulation Setup

p=zeros(M,1);
p=[p+probabilityS;0];

pD=zeros(M,1);
pD=[pD+probabilityD;0]; %Making sure that "probability" has the right size

A=zeros(M,1);
A =[A+alpha;0]; %Making sure that "alpha" has the right size

C=zeros(M,1);
C =[C+cost;0]; %Making sure that "cost" has the right size

Arr=zeros(M,1);
Arr =[Arr+arrival;0]; %Making sure that "arrival" has the right size

Per = Period;

optimal_cost=0;

%% Simulation

real_q=zeros(M,num_iterations);
AoI=zeros(M,num_iterations);

for iteration=1:1:num_iterations

    delivered=zeros(M+1,1);    
    h_history=h1*ones(M,K);
    h_BS=h1*ones(M,1);
    z=z1*ones(M,K);
    vector_successful_tx=zeros(M,1);
    arrival_factor = zeros(M,1);
    interarrivalfactor = ones(M,1);
    servenode = zeros(N,K);
    for slot_index=1:1:K       

        for node_index=1:1:M

                if prob_index == 1
                    if rand<Arr(node_index)
                        z(node_index,slot_index)=0;
                        vector_successful_tx(node_index)=0;
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
                        elseif prob_index == 4
                            interarrivalfactor(node_index) = round(normrnd(Gaussin(node_index,1), Gaussin(node_index,2)));
                        end
                        vector_successful_tx(node_index)=0;% new arrival
                    elseif slot_index>1
                        z(node_index,slot_index)=z(node_index,slot_index-1)+1;
                    end
                end

        end

        % AoI knowledge at the BS
        h_BS(:) = h_history(:,slot_index);

        if set_Policy==1 %Max Weight                 
                [~, sorted_indices] = sort(sqrt(A(1:M).*p(1:M).*pD(1:M)).*( h_BS-z(:,slot_index)),'descend');
                servenode(:,slot_index) = sorted_indices(1:N);
        end       

        AoI(:,iteration)=AoI(:,iteration)+A(1:M).*h_history(:,slot_index);         
        h_history(:,slot_index+1)=h_history(:,slot_index)+1; 
        % receive process at the destination
        for i=1:1:N
            if slot_index > Delay
                node_index = servenode(i,slot_index-Delay);
                if (node_index~=(M+1)) && (rand<p(node_index)*pD(node_index)) && (vector_successful_tx(node_index)==0)
                    delivered(node_index)=delivered(node_index)+1;
                    vector_successful_tx(node_index)=1;                
                    h_history(node_index,slot_index+1)=z(node_index,slot_index-Delay)+Delay;
                end  
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

