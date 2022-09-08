clc
clear
% 不同user数目下四种算法的对比 offloading_energy
MC_OO_list = zeros(1,10);
MC_NO_list = zeros(1,10);
MC_RO_list = zeros(1,10);
MC_proposed_list = zeros(1,10);
% consensus_energy
MC_OO_list_con = zeros(1,10);
MC_NO_list_con = zeros(1,10);
MC_RO_list_con = zeros(1,10);
MC_proposed_list_con = zeros(1,10);
epsilon_value = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5];

% 保证每次迭代参数一致
M = 5; % AP个数
C_m = randi([5,20],1,M);  % AP m的计算资源(CPU频率-GHz)
B_m = randi([10,50],1,M);   % AP m的带宽资源(MHz)
H = sqrt(1/2) * (randn(M,M) + 1i * randn(M,M));  % AP间的信道矩阵
Task_L_0 = randi([100,200],1,1); %User的Task数据量大小(MB)
Task_rho_0 = rand(1,1)*10; % Task的计算密度(cycles*MB)
Task_D_th_0 = rand(1,1)*10; % Task的容忍时延(s)

for index = 1:1:10
    
% Initilazation
%% Parameters Setting
% parameters of network

N = 5; % 用户个数

epsilon_p = 1;
epsilon_c = epsilon_p * epsilon_value(index);
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

d = zeros(M,M);  % 存放AP间距离的矩阵
AP_SINR = zeros(M,M);
P_ap = 1;  % AP通信时的发射功率
P_I = 0.2; % AP通信时其他AP的干扰功率

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
Task_L = ones(1,N)*Task_L_0; %User的Task数据量大小(MB)
Task_rho = ones(1,N)*Task_rho_0; % Task的计算密度(cycles*MB)
Task_D_th = ones(1,N)*Task_D_th_0; % Task的容忍时延(s)
%Task_O = randi([75,100],1,N);


%parametes of ADMM
theta = 0.01; %迭代停止条件
q = 100; % 拉格朗日惩罚因子初始值
t = 1; %迭代次数index
lambda = zeros(M,N); %拉格朗日因子初始值
t_max = 20;
B_min = 0.0000001;
C_min = 0.0000001;

%parameters of blockchain/RAFT
S_b = 500;  % 区块数据大小（KB）
S_com = 10; % 通信数据大小（KB）
D_g = 1;  % 区块间隔(S)
R_max = max(C_m)/S_b + max(B_m)/(S_com+S_b);
P_Leader = rand(1,1);
P_Follower = 1 - P_Leader;



%% Initialization
% 首先生成一个初始的A
A_global = zeros(M,N);  % 分簇矩阵
A_local = zeros(M,N); % 分簇矩阵的本地副本
for i=1:N  % 每一列生成和为1的M个随机数
    x = rand(1,M);
    y = sum(x);
    r = x/y;
    A_local(:,i) = r;
end
C_re_global = zeros(M,N);  % 计算资源决策变量的倒数
B_re_global = zeros(M,N);  % 带宽资源决策变量的倒数
Chi_global = zeros(M,N); % 辅助变量\chi =  A * B_re
Psi_global = zeros(M,N);  % 辅助变量\psi = A * C_re

%% Optocal Offloading
while  t<=t_max  
%% Global Variables Update
norm_final = norm(A_global - A_local);



% 计算用户SINR
User_SINR = User_SINR_Caculate(M,N,A_local,An);
V_result = zeros(1,M);
V_result_con = zeros(1,M);
for i=1:M  %每个AP都求解一个优化问题
% CVX 求解全局变量优化问题
cvx_begin quiet
    % 优化变量
    variable A(1,N)
    variable B_re(1,N)
    variable C_re(1,N)
    variable Chi(1,N)
    variable Psi(1,N)

    % 目标函数
    E_user_up = epsilon_c * sum(Task_L / log2(ones(1,N) + User_SINR)*Chi(1,:));
    E_ap_edge = epsilon_p * sum(Task_L.*Task_rho.*Psi(1,:));
    
    
    V = E_user_up + E_ap_edge;
    L = sum(lambda(i,:) .* (A - A_local(i,:))) + q/2 * sum((A-A_local(i,:)).^2);
    %O =  sum(A .* Task_O);
    minimize(V+L);
    subject to
        Chi >= A / B_m(i);  % (33)
        Chi <= B_re + A/B_m(i) + 1/B_m(i);
        Chi <= A / B_min;
        Chi >= A/B_min - 1/B_min + B_re;
        
        Psi >= A/C_m(i);  % (34)
        Psi <= C_re + A/C_m(i) - 1/C_m(i);
        Psi <= A / C_min;
        Psi >= A/C_min - 1/C_min + C_re;

        
        Task_L ./ log2(1+User_SINR) .* Chi + Task_L .* Task_rho .* Psi <= Task_D_th; %(37)
        sum(inv_pos(B_re)) <= B_m(i);  %(25c-25e)
        sum(inv_pos(C_re)) <= C_m(i);
        for j = 1:N
        1/B_m(i) <= B_re(1,j) <= 1/B_min;
        1/C_m(i) <= C_re(1,j) <= 1/C_min;
        end
        
        0 <= A <= 1; 
cvx_end


% 共识阶段能量消耗
Kappa = (C_m(i) - sum(1./C_re)) / (B_m(i) - sum(1./C_re));
E_ap_leader_block = epsilon_p / R_max + epsilon_p*S_b / ((S_b+S_com)*R_max)*(inv_pos(Kappa));
    
AP_SINR_sum = log2(1+AP_SINR* 10^8);
AP_SINR_sum = sum(AP_SINR_sum(i,:));
AP_SINR_average = AP_SINR_sum/(M-1);

E_ap_leader_com1_1 = 3*epsilon_c / (S_b*R_max) * 3*(M-1)*(S_com+S_b)/ AP_SINR_sum * (inv_pos(Kappa));
E_ap_leader_com2_1 = (M-1)*(S_com+S_b)/ AP_SINR_sum * (inv_pos(Kappa));
E_ap_leader_com1_2 = 3*epsilon_c/ (S_b+S_com)*R_max * 3*(M-1)*(S_com+S_b) / AP_SINR_sum;
E_ap_leader_com2_2 = epsilon_c/ R_max * (M-1) / AP_SINR_sum;
E_ap_leader = E_ap_leader_block + E_ap_leader_com1_1 + E_ap_leader_com2_1 + E_ap_leader_com1_2 + E_ap_leader_com2_2;
E_ap_follower_com_1 = S_com / log2(1+AP_SINR_average) * inv_pos((B_m(1,i)-sum(inv_pos(B_re))));
E_ap_follower_com_2 = epsilon_c*S_com/(S_b*R_max*log2(1+AP_SINR_average)) * (inv_pos(Kappa));
E_ap_follower_com_3 = epsilon_c*S_com/((S_com+S_b)*R_max*log2(1+AP_SINR_average));
E_ap_follower = E_ap_follower_com_1 + E_ap_follower_com_2 + E_ap_follower_com_3;

V_result(i) = V;  
V_result_con(i) = E_ap_leader + E_ap_follower;

A_global(i,:) = A;
B_re_global(i,:) = B_re;
C_re_global(i,:) = C_re;
Chi_global(i,:) = Chi;
Psi_global(i,:) = Psi;
E_edge  = E_user_up + E_ap_edge;
E_com = E_ap_leader + E_ap_follower;
end

V_sum = sum(V_result);
V_sum_con = sum(V_result_con);
MC_OO_list(1,index) = V_sum
MC_OO_list_con(1,index) = V_sum_con
%% Local Variables Update
cvx_begin quiet
variable A_local_copy(M,N)
J_1 = lambda .* (A_global-A_local_copy);
J_2 = q/2 .* ((A_global - A_local_copy).^2);
J = sum(J_1(:)) + sum(J_2(:));
minimize(J);
subject to
sum(A_local_copy) == 1;
0 <= A_local_copy <= 1; 

cvx_end
A_local = A_local_copy;

%% ﻿Dual Variable update
lambda = lambda + q * (A_global - A_local);
t = t + 1;




end

t=1;
%% Single_AP
%caculate cluster matrix
A = zeros(M,N);  %"分簇"矩阵
d_user_ap = zeros(M,N);
for i = 1:N
    for j = 1:M
        d_user_ap(j,i) = sqrt((x_user(1,i) - x_ap(1,j))^2+(y_user(1,i) - y_ap(1,j))^2);
    end
end

for i=1:N
    [max_d,max_index] = max(d_user_ap);
    i_index = max_index(1,i);
    A(i_index,i) = 1;
end

Task_L_value = zeros(M,N); %User的Task数据量大小(KB)
Task_rho_value = zeros(M,N); % Task的计算密度(10^6cycles*KB)
for i=1:M
    Task_rho_value(i,:) = Task_rho;
    Task_L_value(i,:) = Task_L;
end
V_com_total = 0;
% 直接用CVX求解B,C
cvx_begin quiet
    % 优化变量
    variable B(M,N)
    variable C(M,N)
    % 目标函数
    
    E_user_up = epsilon_c * Task_L_value .* A ./ log2(1 + User_SINR) .* inv_pos(B) ; % MxN
    E_ap_edge = epsilon_p * Task_L_value .* Task_rho_value.* A .* inv_pos(C);  % MxN
    for i=1:M
    E_ap_leader_block = P_Leader *  S_b * inv_pos(C_m(i) - sum(C(i,:)));
    
    AP_SINR_sum = log2(1+AP_SINR* 10^8);
    AP_SINR_sum = sum(AP_SINR_sum(i,:));
    AP_SINR_average = AP_SINR_sum/(M-1);

    E_ap_leader_com1 = 3*epsilon_c  * 3*(M-1)*(S_com+S_b)/ AP_SINR_sum * inv_pos(B_m(i) - sum(B(i,:)));
    E_ap_leader_com2 = (M-1)*(S_com+S_b)/ AP_SINR_sum * inv_pos(B_m(i) - sum(B(i,:)));
    E_ap_leader = E_ap_leader_block + E_ap_leader_com1 + E_ap_leader_com2;
    
    E_ap_follower_com = S_com / log2(1+AP_SINR_average) * inv_pos(B_m(i) - sum(B(i,:)));
    E_ap_follower = E_ap_follower_com;
    V_com = E_ap_leader + E_ap_follower;
    V_com_total = V_com_total + V_com;
    % L = sum(lambda(i,:) .* (A - A_local(i,:))) + q/2 * sum((A-A_local(i,:)).^2);
    % O =  sum(A .* Task_O);
    end
    V_edge = sum(E_user_up(:)) + sum(E_ap_edge(:));
    minimize(V_edge + V_com_total);
    subject to
%for j=1:N
%       max(Task_L(j) / log2(1+User_SINR(j)) .* inv_pos(B(:,j)) .* A(:,j) + Task_L_value(j) * Task_rho(j) .* inv_pos(C(:,j)) .* A(:,j)) <= Task_D_th(j); %(37)
%end
for i=1:M
        sum(B(i,:)) <= B_m(i);  %(25c-25e)
        sum(C(i,:)) <= C_m(i);
end
        
cvx_end

MC_NO_list(1,index) = V_edge
MC_NO_list_con(1,index) = V_com_total


t=1;
%% RO
while t<=t_max
%% Global Variables Update
norm_final = norm(A_global - A_local);



% 计算用户SINR
User_SINR = User_SINR_Caculate(M,N,A_local,An);
V_result = zeros(1,M);
for i=1:M  %每个AP都求解一个优化问题
% CVX 求解全局变量优化问题
cvx_begin quiet
    % 优化变量
    variable A(1,N)
    variable B_re(1,N)
    variable C_re(1,N)
    variable Chi(1,N)
    variable Psi(1,N)

    %% 目标函数
    % 卸载阶段能量消耗
    E_user_up = epsilon_c * sum(Task_L / log2(ones(1,N) + User_SINR)*Chi(1,:));
    E_ap_edge = epsilon_p * sum(Task_L.*Task_rho.*Psi(1,:));
    
    % 共识阶段能量消耗
    E_ap_leader_block = P_Leader *  S_b * inv_pos(C_m(i) - sum(inv_pos(C_re)));
    
    AP_SINR_sum = log2(1+AP_SINR* 10^8);
    AP_SINR_sum = sum(AP_SINR_sum(i,:));
    AP_SINR_average = AP_SINR_sum/(M-1);
    
    E_ap_leader_com1 = 3*epsilon_c  * 3*(M-1)*(S_com+S_b)/ AP_SINR_sum * inv_pos(B_m(i) - sum(inv_pos(B_re)));
    E_ap_leader_com2 = (M-1)*(S_com+S_b)/ AP_SINR_sum * inv_pos(B_m(i) - sum(inv_pos(B_re)));
    E_ap_leader = E_ap_leader_block + E_ap_leader_com1 + E_ap_leader_com2;
    
    E_ap_follower_com = S_com / log2(1+AP_SINR_average) * inv_pos(B_m(i) - sum(inv_pos(B_re)));
    E_ap_follower = E_ap_follower_com;

    E_edge  = E_user_up + E_ap_edge;
    E_com = E_ap_leader + E_ap_follower;
    V = E_edge + E_com;
    L = sum(lambda(i,:) .* (A - A_local(i,:))) + q/2 * sum((A-A_local(i,:)).^2);
    %O =  sum(A .* Task_O);
    minimize(V+L);
    subject to
        Chi >= A / B_m(i);  % (33)
        Chi <= B_re + A/B_m(i) + 1/B_m(i);
        Chi <= A / B_min;
        Chi >= A/B_min - 1/B_min + B_re;
        
        Psi >= A/C_m(i);  % (34)
        Psi <= C_re + A/C_m(i) - 1/C_m(i);
        Psi <= A / C_min;
        Psi >= A/C_min - 1/C_min + C_re;

        
        Task_L ./ log2(1+User_SINR) .* Chi + Task_L .* Task_rho .* Psi <= Task_D_th; %(37)
        sum(inv_pos(B_re)) <= B_m(i);  %(25c-25e)
        sum(inv_pos(C_re)) <= C_m(i);
        for j = 1:N
        1/B_m(i) <= B_re(1,j) <= 1/B_min;
        1/C_m(i) <= C_re(1,j) <= 1/C_min;
        end
        
        0 <= A <= 1; 
cvx_end
V_result(i) = E_edge;
V_result_con(i) = E_com;
A_global(i,:) = A;
B_re_global(i,:) = B_re;
C_re_global(i,:) = C_re;
Chi_global(i,:) = Chi;
Psi_global(i,:) = Psi;

end
V_sum = sum(V_result);
V_sum_con = sum(V_result_con);
MC_RO_list(1,index) = V_sum
MC_RO_list_con(1,index) = V_sum_con
L_result = L;
%% Local Variables Update
cvx_begin quiet
variable A_local_copy(M,N)
J_1 = lambda .* (A_global-A_local_copy);
J_2 = q/2 .* ((A_global - A_local_copy).^2);
J = sum(J_1(:)) + sum(J_2(:));
minimize(J);
subject to
sum(A_local_copy) == 1;
0 <= A_local_copy <= 1; 

cvx_end
A_local = A_local_copy;

%% ﻿Dual Variable update
lambda = lambda + q * (A_global - A_local);
t = t + 1;


end

t=1;
%% proposed

C_global = zeros(M,N);  % 计算资源决策变量
B_global = zeros(M,N);  % 带宽资源决策变量
Chi_global = zeros(M,N); % 辅助变量\chi
Psi_global = zeros(M,N);  % 辅助变量\psi
Kappa_global = zeros(M,1); %辅助变量\kappa

zeta = 10;
R_max = max(C_m)/S_b + max(B_m)/(S_com+S_b);
while t<=t_max
%% Global Variables Update
norm_final = norm(A_global - A_local);


% 计算用户SINR
User_SINR = User_SINR_Caculate(M,N,A_local,An);
for i=1:M  %每个AP都求解一个优化问题

% D.C programming initialization B_re(1,N),C_re(1,N)
B_re_0 = ones(1,N) * (2*N/B_m(i));
C_re_0 = ones(1,N) * (2*N/C_m(i));
t_max_DC = 1; 

for t_DC = 1:t_max_DC
gradient_C_re = sum(-1/B_m(i) * (C_re_0).^2);
gradient_B_re = -C_m(i)-inv_pos(B_re_0)./ (B_m(i)- sum(inv_pos(B_re_0)))^2;
% CVX 求解D.C全局变量优化问题
cvx_begin quiet
    variable A(1,N)
    variable B_re(1,N)
    variable C_re(1,N)
    variable Chi(1,N)
    variable Psi(1,N)
    variable Kappa
    %variable E_user_up
    %variable E_ap_edge
    %variable E_ap_con
    % 目标函数
    E_user_up = epsilon_c * sum(Task_L / log2(ones(1,N) + User_SINR)*Chi(1,:));
    E_ap_edge = epsilon_p * sum(Task_L.*Task_rho.*Psi(1,:));
    E_ap_leader_block = epsilon_p / R_max + epsilon_p*S_b / ((S_b+S_com)*R_max)*(inv_pos(Kappa));
    
    AP_SINR_sum = log2(1+AP_SINR*10^8);
    AP_SINR_sum = sum(AP_SINR_sum(i,:));
    AP_SINR_average = AP_SINR_sum/(M-1);

    E_ap_leader_com1_1 = 3*epsilon_c / (S_b*R_max) * 3*(M-1)*(S_com+S_b)/ AP_SINR_sum * (inv_pos(Kappa));
    E_ap_leader_com2_1 = (M-1)*(S_com+S_b)/ AP_SINR_sum * (inv_pos(Kappa));
    E_ap_leader_com1_2 = 3*epsilon_c/ (S_b+S_com)*R_max * 3*(M-1)*(S_com+S_b) / AP_SINR_sum;
    E_ap_leader_com2_2 = epsilon_c/ R_max * (M-1) / AP_SINR_sum;
    E_ap_leader = E_ap_leader_block + E_ap_leader_com1_1 + E_ap_leader_com2_1 + E_ap_leader_com1_2 + E_ap_leader_com2_2;
    E_ap_follower_com_1 = S_com / log2(1+AP_SINR_average) * inv_pos((B_m(1,i)-sum(inv_pos(B_re))));
    E_ap_follower_com_2 = epsilon_c*S_com/(S_b*R_max*log2(1+AP_SINR_average)) * (inv_pos(Kappa));
    E_ap_follower_com_3 = epsilon_c*S_com/((S_com+S_b)*R_max*log2(1+AP_SINR_average));
    E_ap_follower = E_ap_follower_com_1 + E_ap_follower_com_2 + E_ap_follower_com_3;
    
    Penalty = zeta * ((C_m(i) - sum(inv_pos(C_re_0)))/B_m(i) - C_m(i)*inv_pos((B_m(i) - sum(inv_pos(B_re_0)))) + C_m(i)/B_m(i));
    Penalty = -Penalty;

    V = E_user_up + E_ap_edge + E_ap_leader + E_ap_follower;
    L = sum(lambda(i,:) .* (A - A_local(i,:))) + q/2 * sum((A-A_local(i,:)).^2);

    DC_penalty_B_re = zeta * sum(gradient_B_re .* (B_re - B_re_0));
    DC_penalty_C_re = zeta * sum(gradient_C_re * (C_re - C_re_0));
    minimize(V+L-Penalty-DC_penalty_C_re - DC_penalty_B_re);
    subject to
        Chi >= A / B_m(i);  % (33)
        Chi  <= B_re + A/B_m(i) - 1/B_m(i);
        Chi <= A / B_min;
        Chi >= A/B_min - 1/B_min + B_re;
        
        Psi >= A/C_m(i);  % (34)
        Psi <= C_re + A/C_m(i) - 1/C_m(i);
        Psi <= A / C_min;
        Psi >= A/C_min - 1/C_min + C_re;

        %Kappa >= (C_m(i) - sum(inv_pos(C_re)))/B_m(i);  % (35)
        %Kappa <= C_m(i)*inv_pos((B_m - sum(inv_pos(B_re)))) - C_m(i)/B_m(i);
        %Kappa <= (C_m(i)-sum(inv_pos(C_re))) / B_min;
        %Kappa >= (C_m(i)-sum(inv_pos(C_re))) / B_min + C_m(i) * inv_pos((B_m(i) - sum(inv_pos(B_re)))) - C_m/B_min;
        
        Task_L / log2(1+User_SINR) .* Chi + Task_L .* Task_rho .* Psi <= Task_D_th; %(37)
        sum(inv_pos(B_re)) <= B_m(i);  %(25c-25e)
        sum(inv_pos(C_re)) <= C_m(i);
        0 <= A <= 1; 
        
       

        
cvx_end

C_re_0 = C_re;
B_re_0 = B_re;
A_0 =A;
Chi_0 = Chi;
Psi_0 = Psi;
Kappa_0 = Kappa;
end


A_global(i,:) = A;
B_global(i,:) = B_re;
C_global(i,:) = C_re;
Chi_global(i,:) = Chi;
Psi_global(i,:) = Psi;
Kappa_global(i,1) = Kappa;
V_result(i) =  E_user_up + E_ap_edge;
V_result_con(i) = E_ap_leader + E_ap_follower;
end
V_sum = sum(V_result);
V_sum_con = sum(V_result_con);
MC_proposed_list(1,index) = V_sum
MC_proposed_list_con(1,index) = V_sum_con
%% Local Variables Update
cvx_begin quiet
variable A_local_copy(M,N)
J_1 = lambda .* (A_global-A_local_copy);
J_2 = q/2 .* ((A_global - A_local_copy).^2);
J = sum(J_1(:)) + sum(J_2(:));
minimize(J);
subject to
sum(A_local_copy) == 1;
0 <= A_local_copy <= 1; 

cvx_end
A_local = A_local_copy;

%% ﻿Dual Variable update
lambda = lambda + q * (A_global - A_local);
t = t + 1;


end

end

IP_C_energy_list = zeros(4,10); % OO,NO,RO,proposed

IP_C_energy_list(1,:) =  MC_OO_list_con + MC_OO_list;
IP_C_energy_list(2,:) =  MC_NO_list_con + MC_NO_list;
IP_C_energy_list(3,:) =  MC_RO_list_con + MC_RO_list;
IP_C_energy_list(4,:) =  MC_proposed_list_con + MC_proposed_list;

IP_epsilon_energy_list = IP_C_energy_list;
