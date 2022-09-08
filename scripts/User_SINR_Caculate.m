function User_sinr = User_SINR_Caculate(M,N,A,An)
%Caculate User SINR(1,N)
User_sinr = zeros(1,N);
for i=1:N  % 计算一个特定用户n的SINR
 if A(:,i) == zeros(M,1)
        continue
 end
Phi = find(A(:,i));   % 计算服务该用户的AP簇
Omega = zeros(1,N);  % 用户的簇间用户集合表示矩阵，值为1代表对应序号的用户为该用户的簇内干扰用户
 for j=1:N
     if j==i
         continue
     end
     Phi_n = find(A(:,j));
     if isequal(Phi,Phi_n)
         continue
     else
         Omega(1,j) = 1;  %该用户是目标用户的簇间干扰用户
     end
 end

% 计算干扰矩阵I
clus_size = length(Phi);   %AP簇中的AP个数
I = eye(An*clus_size, An*clus_size);

% 计算信道系数矩阵G(|Phi|*An,|Omega|-1),不包括当前用户
ch = sqrt(1/2) * (randn(M*An,N) + 1i * randn(M*An, N));
G_u = ch;
 for bs=M:-1:1
    if A(bs,i) ==0
        G_u((bs-1)*An+1:bs*An,:) = [];
    end
 end
 for user=N:-1:1
     if Omega(1,user)==0 || user==i
         G_u(:,user) = [];
     end
 end
% 计算g_n(|Phi|*A,U)
g_u = ch;
for bs=M:-1:1
    if A(bs,i) == 0
        g_u((bs-1)*An+1:bs*An,:) = [];
    end
end
% 计算波束形成矢量
temp = (I-G_u*pinv(G_u))*g_u(:,i);
beamform_u = temp./norm(temp);
% test
test = beamform_u' * G_u; %结果应该十分接近0，表示intra干扰消除   
% 求用户u的信干噪比
P = 1; % 用户发射功率
signal = P*abs(beamform_u'*g_u(:,i))^2;
interf = 0;
for user=1:N
    if sum(A(:,user).*A(:,i))==0
        interf = interf + P*abs(beamform_u'*g_u(:,user))^2;
    end
end
noise = 10^(-3.5);
User_sinr(1,i) = signal/(interf+noise^2);



end

end

