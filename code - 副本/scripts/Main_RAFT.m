% clc;
clear;
time = datetime;
version = [num2str(time.Year) num2str(time.Month) num2str(time.Day) num2str(time.Hour)];
%% Parameters Setting
% parameters of network
M = 11; % AP个数
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
Task_L = randi([30,50],1,N); %User的Task数据量大小(KB)
Task_rho = rand(1,N)*0.2; % Task的计算密度(cycles*KB)
Task_D_th = randi([3,5],1,N); % Task的容忍时延(s)
Task_O = randi([75,100],1,N);


%parametes of ADMM
theta = 0.01; %迭代停止条件
q = 0.5; % 拉格朗日惩罚因子初始值
t = 1; %迭代次数index
lambda = zeros(M,N); %拉格朗日因子初始值
t_max = 5;
B_min = 0.0000001;
C_min = 0.0000001;

%parameters of blockchain/RAFT
S_b = 50;  % 区块数据大小（KB）
S_com = 10; % 通信数据大小（KB）
D_g = 0.5;  % 区块间隔(S)
%R_max = max(C_m)/S_b + max(B_m)/(S_com+S_b);
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
    
    E_ap_leader_com1 = P_Leader*3*epsilon_c  * 3*(M-1)*(S_com+S_b)/ AP_SINR_sum * inv_pos(B_m(i) - sum(inv_pos(B_re)));
    E_ap_leader_com2 = P_Leader*(M-1)*(S_com+S_b)/ AP_SINR_sum * inv_pos(B_m(i) - sum(inv_pos(B_re)));
    E_ap_leader = E_ap_leader_block + E_ap_leader_com1 + E_ap_leader_com2;
    
    E_ap_follower_com = S_com / log2(1+AP_SINR_average) * inv_pos(B_m(i) - sum(inv_pos(B_re)));
    E_ap_follower = P_Follower*E_ap_follower_com;

    E_edge  = E_user_up + E_ap_edge;
    E_com = E_ap_leader + E_ap_follower;
    V = E_edge + E_com;
    L = sum(lambda(i,:) .* (A - A_local(i,:))) + q/2 * sum((A-A_local(i,:)).^2);
    O =  sum(A .* Task_O);
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

        
        %Task_L ./ log2(1+User_SINR) .* Chi + Task_L .* Task_rho .* Psi <= Task_D_th; %(37)
        sum(inv_pos(B_re)) <= B_m(i);  %(25c-25e)
        sum(inv_pos(C_re)) <= C_m(i);
        for j = 1:N
        1/B_m(i) <= B_re(1,j) <= 1/B_min;
        1/C_m(i) <= C_re(1,j) <= 1/C_min;
        end
        
        0 <= A <= 1; 
cvx_end


Task_L ./ log2(1+User_SINR) .* Chi + Task_L .* Task_rho .* Psi;
V_result(i) = V;



A_global(i,:) = A;
B_re_global(i,:) = B_re;
C_re_global(i,:) = C_re;
Chi_global(i,:) = Chi;
Psi_global(i,:) = Psi;

end
V_sum = sum(V_result);
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
E_edge + E_com