function randomized_find_optMU(M, N)

global probabilityS probabilityD alpha optMU

iterations = 200;

% 计算合并的概率值
probability = probabilityS .* probabilityD; % 确保"probability"大小正确

% 确保alpha和weight有正确的大小
A = alpha;


% 初始化optMu
optMU = zeros(M, 1);

% 计算初始值
theta = sum(sqrt(probability  ./ A)) / (N^2);
vi = (A ./ (probability));
theta_i = max(vi);

% 计算mu_i
mu = min(1, sqrt(vi / theta_i));
S = sum(mu);

% 使用二分法优化theta
while S > N && theta_i > 0
    theta_i = theta_i * 1.01; % 减少theta_i
    mu = min(1, sqrt(vi / theta_i));
    S = sum(mu);
end

% 最终结果
optMU = mu;
end



% function randomized_find_optMU(M,N)
% 
% global probabilityS probabilityD alpha arrival optMU
% 
% iterations=200;
% 
% p=zeros(M,1);
% p=p+probabilityS.*probabilityD;%Making sure that "probability" has the right size
% 
% A=zeros(M,1);
% A =A+alpha; %Making sure that "alpha" has the right size
% 
% Arr=zeros(M,1);
% Arr=Arr+arrival; %Making sure that "arrival" has the right size
% 
% 
% optMU=zeros(M,1);
% 
% % if sum(Arr./p)<=N
% %     for node=1:M
% %         optMU(node)=Arr(node);
% % %     end
% % else    
%     Gamma_star=(sum(sqrt(p./A)))^2/N/N;
%     QI=A./p;
%     if sum((QI-Arr)>0)==0
%         for node=1:M
%             optMU(node)=QI(node);
%         end
%     else
%         GammaI=A./p;
%         maxGamma=max(GammaI);
%         minGamma=0;
%         %Q_initial=sqrt(A./(M*maxGamma));
%         %S=sum(mu_initial);
%         for count=1:iterations
%             Gamma_aux=(maxGamma+minGamma)/2;
%             mu_aux=zeros(M,1);
%             for node=1:M
%                 mu_aux(node)=min([1,GammaI(node)./maxGamma]);
%                 % Q_aux(node)=sqrt(A(node).*p(node)./(M*Gamma_aux));
%             end    
%             S=sum(mu_aux);
%             if S<=N
%                 maxGamma=Gamma_aux;
%             else
%                 minGamma=Gamma_aux;
%             end
%         end
%         optMU=zeros(M,1);
%         for node=1:M
%             % optQ(node)=sqrt(A(node).*p(node)./(M*Gamma_aux));
%             optMU(node)=min([1,GammaI(node)./maxGamma]);
%         end    
%     end    
% % end
% end