clc;
clear;
time = datetime;
version = [num2str(time.Year) num2str(time.Month) num2str(time.Day) num2str(time.Hour)];
%% Parameters Setting
% parameters of network
M = 5; % AP个数
N = 5; % 用户个数
epsilon_c = 1.5;
epsilon_p = 1;

angle=rand(1,M)*2*pi; %(0,2*pi)之间均匀分布数据点 
r=200*sqrt(rand(1,M)); %（0,1）之间r^2均匀分布数据点
x_ap=r.*cos(angle); %AP的x坐标
y_ap=r.*sin(angle); %AP的y坐标

angle_1=rand(1,N)*2*pi; %(0,2*pi)之间均匀分布数据点 
r_1=200*sqrt(rand(1,N)); %（0,1）之间r^2均匀分布数据点
x_user=r_1.*cos(angle_1); % User的x坐标
y_user=r_1.*sin(angle_1); % USer的y坐标


% parameters of APs
An=2; % 每个AP的天线个数
B_m = randi([10,50],1,M);   % AP m的带宽资源(MHz)
C_m = randi([5,20],1,M);  % AP m的计算资源(CPU频率-GHz)
d = zeros(M,M);  % 存放AP间距离的矩阵
AP_SINR = zeros(M,M);
P_ap = 1;  % AP通信时的发射功率
P_I = 0.2; % AP通信时其他AP的干扰功率
H = sqrt(1/2) * (randn(M,M) + 1i * randn(M,M));  % AP间的信道矩阵
z_2 = 10^(-3.5);  %z_2噪声功率
I_AP = zeros(M,M);
for i =1:M
    for j=1:M
        d(i,j) = sqrt((x_ap(1,i) - x_ap(1,j))^2+(y_ap(1,i) - y_ap(1,j))^2);  
        I_AP(i,j) = P_I*abs(H(i,j))^2 / d(i,j);
        if i == j
            I_AP(i,j) = 0;
        end
    end
end
% 计算AP_SINR
for i = 1:M
    for j=1:M
        if i==j
            AP_SINR(i,j) = 0;
            continue
        end
        AP_SINR(i,j) = (P_ap*abs(H(i,j))^2/d(i,j)) / (sum(I_AP(i,:))-I_AP(i,i)-I_AP(i,j) + (z_2)^2 );
    end
end



% parameters of Tasks
P=1; %用户发送功率
Task_L = rand(1,N)*200; %User的Task数据量大小(MB)
Task_rho = 4000; % Task的计算密度(cycles*MB)
Task_D_th = rand(1,N)*10; % Task的容忍时延(s)



%parametes of ADMM
theta = 0.01; %迭代停止条件
q = 100; % 拉格朗日惩罚因子初始值
t = 1; %迭代次数index
lambda = zeros(M,N); %拉格朗日因子初始值
t_max = 2000;
B_min = 0.0001;
C_min = 0.0001;
zeta = 10;


%parameters of blockchain/RAFT
S_b = 500;  % 区块数据大小（B）
S_com = 10; % 通信数据大小（B）
D_g = 5;  % 区块间隔
R_max = max(C_m)/S_b + max(B_m)/(S_com+S_b);

%% Initialization
% Initialize A and B
A = zeros(M,N);  
for i=1:N  
    x = rand(1,M);
    y = sum(x);
    r = x/y;
    A(:,i) = r;
end

cnt = zeros(1,M);

B= zeros(M,N);
for i=1:M
    for j=1:N
        if A(i,j) ~= 0
            cnt(1,i) = cnt(1,i) +1;
        end
    end
end

for i=1:M
    for j=1:N
        if A(i,j) ~= 0
            B(i,j) = B_m(1,i) / cnt(1,i) / 2;
        end
    end
end

C = zeros(M,N);
iter_index = 0;
iter_max_GBD = 10;

for iter_index =1:iter_max_GBD
    clear Kappa
    % Calculate the SINR
    User_SINR = User_SINR_Caculate(M,N,A,An);

   %  Solving the problem when fixing A & B
    cvx_begin quiet
    variable C_re(M,N)
    % 目标函数


    for i = 1:M
        for j=1:N
            Kappa(i,j) = (C_m(1,i) - sum(C_re(i,:))) /  (B_m(1,i) - sum(B(i,:)));
        end
    end

    E_user_up = sum(epsilon_c * sum((ones(M,1) * Task_L) ./ (ones(M,1)*(log2(ones(1,N) + User_SINR))).*A./B));
    E_ap_edge = sum(epsilon_p * sum((ones(M,1) * (Task_L.*Task_rho)).*A .* inv_pos(C_re)));
    E_ap_leader_block = sum(sum(epsilon_p / R_max + epsilon_p*S_b / ((S_b+S_com)*R_max)*(inv_pos(Kappa))));
    
    AP_SINR_sum = log2(1+AP_SINR);
    AP_SINR_sum = sum(AP_SINR_sum(i,:));
    AP_SINR_average = AP_SINR_sum/(M-1);

    E_ap_leader_com1_1 = sum(sum(3*epsilon_c / (S_b*R_max) * 3*(M-1)*(S_com+S_b)/ AP_SINR_sum * (inv_pos(Kappa))));
    E_ap_leader_com2_1 = sum(sum((M-1)*(S_com+S_b)/ AP_SINR_sum * (inv_pos(Kappa))));
    E_ap_leader_com1_2 = 3*epsilon_c/ (S_b+S_com)*R_max * 3*(M-1)*(S_com+S_b) / AP_SINR_sum;
    E_ap_leader_com2_2 = epsilon_c/ R_max * (M-1) / AP_SINR_sum;
    E_ap_leader = E_ap_leader_block + E_ap_leader_com1_1 + E_ap_leader_com2_1 + E_ap_leader_com1_2 + E_ap_leader_com2_2;
    E_ap_follower_com_1 = sum(S_com / log2(1+AP_SINR_average) * inv_pos((B_m(1,i)-sum(inv_pos(B)))));
    E_ap_follower_com_2 = sum(sum(epsilon_c*S_com/(S_b*R_max*log2(1+AP_SINR_average)) * (inv_pos(Kappa))));
    E_ap_follower_com_3 = epsilon_c*S_com/((S_com+S_b)*R_max*log2(1+AP_SINR_average));
    E_ap_follower = E_ap_follower_com_1 + E_ap_follower_com_2 + E_ap_follower_com_3;
    
    

    E_total = E_user_up + E_ap_edge + E_ap_leader + E_ap_follower ;
    
    D_off = sum((ones(M,1) * Task_L) ./ (ones(M,1) * log2(1+User_SINR)) .* (A./B) + (ones(M,1) * Task_L).* (ones(M,1) * Task_rho) .* A .* inv_pos(C_re));
    minimize(E_total);
    subject to      
      %D_off   <= Task_D_th; 
      sum(C_re,2) <= C_m';
      0 <= C_re; 
        
        
       

        
cvx_end

C = C_re;


end





   %  Solving the problem when fixing B & C
    cvx_begin quiet
    variable A_re(M,N)
    % 目标函数


    for i = 1:M
        for j=1:N
            Kappa(i,j) = (C_m(1,i) - sum(C_re(i,:))) /  (B_m(1,i) - sum(B(i,:)));
        end
    end

    E_user_up = sum(epsilon_c * sum((ones(M,1) * Task_L) ./ (ones(M,1)*(log2(ones(1,N) + User_SINR))).*A./B));
    E_ap_edge = sum(epsilon_p * sum((ones(M,1) * (Task_L.*Task_rho)).*A .* inv_pos(C_re)));
    E_ap_leader_block = sum(sum(epsilon_p / R_max + epsilon_p*S_b / ((S_b+S_com)*R_max)*(inv_pos(Kappa))));
    
    AP_SINR_sum = log2(1+AP_SINR);
    AP_SINR_sum = sum(AP_SINR_sum(i,:));
    AP_SINR_average = AP_SINR_sum/(M-1);

    E_ap_leader_com1_1 = sum(sum(3*epsilon_c / (S_b*R_max) * 3*(M-1)*(S_com+S_b)/ AP_SINR_sum * (inv_pos(Kappa))));
    E_ap_leader_com2_1 = sum(sum((M-1)*(S_com+S_b)/ AP_SINR_sum * (inv_pos(Kappa))));
    E_ap_leader_com1_2 = 3*epsilon_c/ (S_b+S_com)*R_max * 3*(M-1)*(S_com+S_b) / AP_SINR_sum;
    E_ap_leader_com2_2 = epsilon_c/ R_max * (M-1) / AP_SINR_sum;
    E_ap_leader = E_ap_leader_block + E_ap_leader_com1_1 + E_ap_leader_com2_1 + E_ap_leader_com1_2 + E_ap_leader_com2_2;
    E_ap_follower_com_1 = sum(S_com / log2(1+AP_SINR_average) * inv_pos((B_m(1,i)-sum(inv_pos(B)))));
    E_ap_follower_com_2 = sum(sum(epsilon_c*S_com/(S_b*R_max*log2(1+AP_SINR_average)) * (inv_pos(Kappa))));
    E_ap_follower_com_3 = epsilon_c*S_com/((S_com+S_b)*R_max*log2(1+AP_SINR_average));
    E_ap_follower = E_ap_follower_com_1 + E_ap_follower_com_2 + E_ap_follower_com_3;
    
    

    E_total = E_user_up + E_ap_edge + E_ap_leader + E_ap_follower ;
    
    D_off = sum((ones(M,1) * Task_L) ./ (ones(M,1) * log2(1+User_SINR)) .* (A./B) + (ones(M,1) * Task_L).* (ones(M,1) * Task_rho) .* A .* inv_pos(C_re));
    minimize(E_total);
    subject to      
      %D_off   <= Task_D_th; 
      sum(A_re) == ones(1,N);
      0 <= A_re <= 1; 
              
cvx_end

A = A_re;


clear Kappa
   %  Solving the problem when fixing A & C
        cvx_begin quiet
    variable B_re(M,N)
    % 目标函数


    for i = 1:M
        for j=1:N
            Kappa(i,j) = (C_m(1,i) - sum(C_re(i,:))) /  (B_m(1,i) - sum(B(i,:)));
        end
    end

    E_user_up = sum(epsilon_c * sum((ones(M,1) * Task_L) ./ (ones(M,1)*(log2(ones(1,N) + User_SINR))).*A./B));
    E_ap_edge = sum(epsilon_p * sum((ones(M,1) * (Task_L.*Task_rho)).*A .* inv_pos(C_re)));
    E_ap_leader_block = sum(sum(epsilon_p / R_max + epsilon_p*S_b / ((S_b+S_com)*R_max)*(inv_pos(Kappa))));
    
    AP_SINR_sum = log2(1+AP_SINR);
    AP_SINR_sum = sum(AP_SINR_sum(i,:));
    AP_SINR_average = AP_SINR_sum/(M-1);

    E_ap_leader_com1_1 = sum(sum(3*epsilon_c / (S_b*R_max) * 3*(M-1)*(S_com+S_b)/ AP_SINR_sum * (inv_pos(Kappa))));
    E_ap_leader_com2_1 = sum(sum((M-1)*(S_com+S_b)/ AP_SINR_sum * (inv_pos(Kappa))));
    E_ap_leader_com1_2 = 3*epsilon_c/ (S_b+S_com)*R_max * 3*(M-1)*(S_com+S_b) / AP_SINR_sum;
    E_ap_leader_com2_2 = epsilon_c/ R_max * (M-1) / AP_SINR_sum;
    E_ap_leader = E_ap_leader_block + E_ap_leader_com1_1 + E_ap_leader_com2_1 + E_ap_leader_com1_2 + E_ap_leader_com2_2;
    E_ap_follower_com_1 = sum(S_com / log2(1+AP_SINR_average) * inv_pos((B_m(1,i)-sum(inv_pos(B)))));
    E_ap_follower_com_2 = sum(sum(epsilon_c*S_com/(S_b*R_max*log2(1+AP_SINR_average)) * (inv_pos(Kappa))));
    E_ap_follower_com_3 = epsilon_c*S_com/((S_com+S_b)*R_max*log2(1+AP_SINR_average));
    E_ap_follower = E_ap_follower_com_1 + E_ap_follower_com_2 + E_ap_follower_com_3;
    
    

    E_total = E_user_up + E_ap_edge + E_ap_leader + E_ap_follower ;
    
    D_off = sum((ones(M,1) * Task_L) ./ (ones(M,1) * log2(1+User_SINR)) .* (A./B) + (ones(M,1) * Task_L).* (ones(M,1) * Task_rho) .* A .* inv_pos(C_re));
    minimize(E_total);
    subject to      
      %D_off   <= Task_D_th; 
      sum(B_re,2) <= B_m';
      0 <= B_re; 
        
        
     B = B_re;
     % E_offloading = (E_user_up + E_ap_edge) / 10000
     % E_consensus = (E_ap_leader + E_ap_follower) / 10000
     E_total / 10000
     clear Kappa
        
cvx_end
% consensus delay 



