
close all; clc;clear;

%%%%%%%%%%%%%% 确定位置 %%%%%%%%%%%
x_BS = 0;  % BS的x坐标
y_BS = 0;  % BS的y坐标
x_RIS = 30;  % RIS的x坐标
y_RIS = 5;   % RIS的y坐标
x_t = 60;   % 目标的x坐标
y_t = 0;    % 目标的y坐标
% 设定用户的位置，用户的坐标在圆形区域内均匀随机生成
radius = 10;  % 半径为10米
center_x = 10; % 圆心x坐标
center_y = 30; % 圆心y坐标
K = 4; % 用户数量

% 随机生成用户位置
kappa = 2 * pi * rand(1, K);  % 用户的角度
r = radius * sqrt(rand(1, K));  % 用户到圆心的距离，均匀分布

x_USER = center_x + r .* cos(kappa);  % 用户x坐标
y_USER = center_y + r .* sin(kappa);  % 用户y坐标

%%%%%  距离、相位路径损耗、功率、rician因子参数 %%%%%%%%%
dBR = sqrt((x_BS-x_RIS)^2+(y_BS-y_RIS)^2);%BS-RIS
dRk = sqrt((x_RIS-x_USER).^2+(y_RIS-y_USER).^2);%RIS-User
dBk = sqrt((x_BS-x_USER).^2+(y_BS-y_USER).^2);%BS-User
dBt = sqrt((x_t - x_BS)^2 + (y_t - y_BS)^2);  % BS-target
dRt = sqrt((x_t - x_RIS)^2 + (y_t - y_RIS)^2);  % RIS-target
N_r = 6;%Bs ant number
N_t = 10;%Bs ant number
L = 8;  % 预设不同 RIS 反射面数量
L2 = 16;
L3 = 32;
K=4;

delta_0 = 10^(-3);%PL(d) 
d_0 = 1;%baseline
alpha_Rk = 2.2;%path loss exponents
alpha_BR = 2;
alpha_Bk = 3.2;
alpha_Bt = 2.2;
alpha_Rt = 2;
Rician_factor_dB_1 = 2.2;%Rician factors for BR Rk
Rician_factor_dB_2 = 0;%Rician factors for  Bk 
sigma_k2_db = -114; % dBm
sigma_s2_db = -150; % dBm
sigma_k2 = 10^(sigma_k2_db / 10) * 1e-3; % 转换为瓦特
sigma_s2 = 10^(sigma_s2_db / 10) * 1e-3;
User_ant = 1;%user ant number
C = 3; % 干扰信号的数量
alpha_s = 1; % 目标的复振幅 alpha_s 
alpha_c = 0.8;% clutter的复振幅 alpha_c 

%%%%%   生成信道%%%%%%%%

[PL_BR,PL_Rk,PL_Bk,PL_Bt,PL_Rt] = generate_channel_loss(dBk,dBR,dRk,dBt,dRt, delta_0, alpha_Rk,alpha_BR,alpha_Bk,alpha_Bt,alpha_Rt);
[hR, hB , G] = generate_channel(PL_BR,PL_Rk, PL_Bk, Rician_factor_dB_1,Rician_factor_dB_2, L, N_t,K);
%[hR2, hB2 , G2] = generate_channel2(PL_BR,PL_Rk, PL_Bk, Rician_factor_dB_1,Rician_factor_dB_2, L2, N_t,K);
%[hR3, hB3 , G3] = generate_channel3(PL_BR,PL_Rk, PL_Bk, Rician_factor_dB_1,Rician_factor_dB_2, L3, N_t,K);
%%%%%%%%%% steering_vector %%%%%%%%%%%%
Phi_s = -pi/4;% sensing signal入射角
Phi_c = [60, 70, 80] * pi / 180; %clutter 3个入射角
f = 2.2; %GHz
a_s = steering_vector(Phi_s , N_t );
a_c = zeros(N_t,C);
for k = 1:C
a_c(:,k) = steering_vector(Phi_c(k), N_t );%a_s is target reflection,a_c is clutter reflection
end
lambda_c = 0.5;
lambda_s = 0.5;%weight of s&c 

 %%%%%%%%%%构建beam pattern gain图%%%%%%%%%%%%%%
R = solveR(N_t,Phi_s); %%%%协方差矩阵
                                
load('v_storage.mat'); 
v = v_storage(:, 1); 
Theta = diag(v); % 构造对角矩阵 Θ = diag(exp(v));
%%%%%%%%%%%%%% L=32 %%%%%%%%%%%%%
%v2 = repmat(v, 2, 1);  % 将 v 扩展为 16×1
%Theta2 = diag(v2); % 构造 16×16 矩阵 Theta2
%%%%%%%%%%%%%% L=32 %%%%%%%%%%%%%
%v3 = repmat(v, 4, 1);  % 将 v 扩展为 32×1
%Theta3 = diag(v3); % 构造 32×32 矩阵 Theta3
%%%%%%% 计算 hk %%%%%%%%%%%%%%
h_k = zeros(N_t, K); % 初始化 h_k 
for k=1:K
    h_k(:,k) = (hB(:,k)' + (hR(:,k)' * Theta * G))';    %% theta is random
end
%%%%%%%%%%%%%% L=16 %%%%%%%%%%%%%
%h_k2 = zeros(N_t, K); % 初始化 h_k 
%for k=1:K
   % h_k2(:,k) = (hB2(:,k)' + (hR2(:,k)' * Theta2 * G2))';    %% theta is random
%end
%%%%%%%%%%%%%% L=32 %%%%%%%%%%%%%
%h_k3 = zeros(N_t, K); % 初始化 h_k 
%for k=1:K
 %   h_k3(:,k) = (hB3(:,k)' + (hR3(:,k)' * Theta3 * G3))';    %% theta is random
%end
%%%%%%% 计算 ZF 预编码矩阵 %%%%%%%%%%%%%
w = zeros(N_t,K);
H = h_k';
W = H' * pinv(H * H');
for k = 1:K
    W(:, k) = W(:, k) / norm(W(:, k));
    [U, S, ~] = svd(W);
    w(:, k) = U(:, 1) * sqrt(S(1,1)); 
    w(:, k) = w(:, k) / norm(w(:, k));
end
%%%%%%%%%%%%%% L=16 %%%%%%%%%%%%%
%w2 = zeros(N_t,K);
%H2 = h_k2';
%W2 = H2' * pinv(H2 * H2');
%for k = 1:K
  %  W2(:, k) = W2(:, k) / norm(W2(:, k));
  %  [U, S, ~] = svd(W2);
 %   w2(:, k) = U(:, 1) * sqrt(S(1,1)); 
 %   w2(:, k) = w2(:, k) / norm(w2(:, k));
%end
%%%%%%%%%%%%%% L=32 %%%%%%%%%%%%%
%w3 = zeros(N_t,K);
%H3 = h_k3';
%W3 = H3' * pinv(H3 * H3');
%for k = 1:K
   % W3(:, k) = W3(:, k) / norm(W3(:, k));
   % [U, S, ~] = svd(W3);
  %  w3(:, k) = U(:, 1) * sqrt(S(1,1)); 
   % w3(:, k) = w3(:, k) / norm(w3(:, k));
%end
%%%%%%%%%%% no ris %%%%%%%%%
   %h_k_noris = zeros(N_t, K); % 初始化 h_k 
%for k=1:K
   % h_k_noris(:,k) = (hB(:,k)');
%end
%H_noris = h_k_noris'; 
%W_noris = H_noris' * pinv(H_noris * H_noris'); 
%for k = 1:K
   % W_noris(:, k) = W_noris(:, k) / norm(W_noris(:, k));
%end
%[U, S, V] = svd(W_noris, 'econ');  
%w_noRIS = U * sqrt(S); 

%%%%%%%%%%% 计算辅助变量f_k %%%%%%%%%%%%%
f_k = zeros(K,1); % 初始化 f_k
for k = 1:K
    numerator = h_k(:,k)' * w(:, k); % 计算分子
    denominator = 0; % 初始化分母，包含噪声项
    for i = 1:K
        if i ~= k
            denominator = denominator + abs(h_k(:,k)' * w(:, i))^2; % 计算干扰项
        end
    end
    f_k(k) = numerator / denominator+sigma_k2; % 计算 f_k
end
%%%%%%%%%%%%%% L=16 %%%%%%%%%%%%%
%f_k2 = zeros(K,1); % 初始化 f_k
%for k = 1:K
 %  numerator6 = h_k2(:,k)' * w2(:, k); % 计算分子
  %  denominator6 = 0; % 初始化分母，包含噪声项
   % for i = 1:K
       % if i ~= k
       %     denominator6 = denominator6 + abs(h_k2(:,k)' * w2(:, i))^2; % 计算干扰项
     %   end
 %   end
 %   f_k2(k) = numerator6 / denominator6+sigma_k2; % 计算 f_k
%end
%%%%%%%%%%%%%% L=32 %%%%%%%%%%%%%
%f_k3 = zeros(K,1); % 初始化 f_k
%for k = 1:K
   % numerator7 = h_k3(:,k)' * w3(:, k); % 计算分子
   % denominator7 = 0; % 初始化分母，包含噪声项
  % for i = 1:K
      %  if i ~= k
          %  denominator7 = denominator7 + abs(h_k3(:,k)' * w3(:, i))^2; % 计算干扰项
     %   end
   % end
  %  f_k3(k) = numerator7 / denominator7+sigma_k2; % 计算 f_k
%end

%%%%%%%%%%% no ris f_k %%%%%%%%%%%%%%%
%f_k_noris = zeros(K,1); % 初始化f_k_noris
%for k = 1:K
   % numerator = h_k_noris(:,k)' * w_noRIS(:, k); % 计算分子
  %  denominator = 0; % 初始化分母，包含噪声项
   % for i = 1:K
       % if i ~= k
      %     denominator = denominator + abs(h_k_noris(:,k)' * w_noRIS(:, i))^2; % 计算干扰项
       % end
   % end
   % f_k_noris(k) = numerator / denominator+sigma_k2; % 计算 f_k_noris
%end
%%%%%%%%%%% 计算辅助变量f_K+1 %%%%%%%%%%%%%%%%%%%
denominator3 = 0;
for k = 1:K
    numerator3 = alpha_s * (a_s)' * w(:,k);
    for c = 1:C
    denominator3 = denominator3 + abs(alpha_c * (a_c(:,c))' * w(:,k))^2;
    end
    f_K1 = numerator3 / (denominator3 + sigma_s2); 
end
%%%%%%%%%%%%%% L=16 %%%%%%%%%%%%%
%denominator8 = 0;
%for k = 1:K
  %  numerator8 = alpha_s * (a_s)' * w2(:,k);
  %  for c = 1:C
  %  denominator8 = denominator8 + abs(alpha_c * (a_c(:,c))' * w2(:,k))^2;
   % end
   % f_K1_2 = numerator8 / (denominator8 + sigma_s2); 
%end
%%%%%%%%%%%%%% L=32 %%%%%%%%%%%%%
%denominator9 = 0;
%for k = 1:K
 %   numerator9 = alpha_s * (a_s)' * w3(:,k);
 %   for c = 1:C
 %  denominator9 = denominator9 + abs(alpha_c * (a_c(:,c))' * w3(:,k))^2;
 %   end
 %   f_K1_3 = numerator9 / (denominator9 + sigma_s2); 
%end
%%%%%%%%%%% 计算辅助变量f_K+1_noris %%%%%%%%%%%%%%%%%%%
%denominator2 = 0;
%for k = 1:K
   % numerator2 = alpha_s * (a_s)' * w_noRIS(:,k);
  %  for c = 1:C
  % denominator2 = denominator2 + abs(alpha_c * (a_c(:,c))' * w_noRIS(:,k))^2;
   % end
  %  f_K1_noris = numerator2 / (denominator2 + sigma_s2); 
%end

%%%%%%%%%%% AO 算法 %%%%%%%%%%%
iter = 1;
iter_max = 15;
tau_3 = 0.01; % 终止阈值
P_0_db = 30; % dBm (总功率)  
P_0 = 10^(P_0_db / 10) * 1e-3;
WSR_values = zeros(iter_max, 1); % 记录 WSR 变化
while(iter<=iter_max)
epsilon = 0.01;
 obj_prev = 0; 
%%%%%%%%% 求解W %%%%%%%%%%%%%
W = solveW(hB,hR,G,P_0,R,N_t,C,K,L,lambda_c,lambda_s,f_k,f_K1,v,sigma_k2,sigma_s2,a_c,a_s,alpha_c,alpha_s,epsilon);
%WnoRIS = solveWnoRIS(R,hB,P_0,N_t,C,K,lambda_c,lambda_s,f_k_noris,f_K1_noris,sigma_k2,sigma_s2,a_c,a_s,alpha_c,alpha_s,epsilon);
%W2 = solveW2(h_k2,P_0,R,N_t,C,K,lambda_c,lambda_s,f_k2,f_K1_2,sigma_k2,sigma_s2,a_c,a_s,alpha_c,alpha_s,epsilon);
%W3 = solveW3(h_k3,P_0,R,N_t,C,K,lambda_c,lambda_s,f_k3,f_K1_3,sigma_k2,sigma_s2,a_c,a_s,alpha_c,alpha_s,epsilon);
%%%%%%%% W 变 w%%%%%%%%%
w = zeros(N_t, K);  % 初始化 w
sigma_w = zeros(N_t, N_t);  % 预分配 sigma_w

for k = 1:K
    slice = W(:, :, k);  % 取 W 的第 k 个 31x31 矩阵
    [V, D] = eig(slice);  % 计算特征值和特征向量
    [max_val, idx] = max(diag(D));  % 提取最大特征值
    %w_title = sqrt(max_val) * V(:, idx);  % 取最大特征值对应的特征向量
   w(:, k) = sqrt(max_val) * V(1:N_t, idx);
% w(:, k) = w_title(1:N_t) / w_title(N_t + 1);  % 归一化，使其最后一维为 1
    sigma_w =  w(:, k) * w(:, k)';  % 计算 sigma_w
end 
rank_approx = max_val / sum(diag(D));
disp(['秩 1 近似程度: ', num2str(rank_approx)]);
%%%%%%%%%%%%%% L=16 %%%%%%%%%%%%%
%w2 = zeros(N_t, K);  % 初始化 w
%for k = 1:K
   % slice = W2(:, :, k);  % 取 W 的第 k 个 31x31 矩阵
  %  [V, D] = eig(slice);  % 计算特征值和特征向量
  %  [max_val, idx] = max(diag(D));  % 提取最大特征值
  %  w_title = sqrt(max_val) * V(:, idx);  % 取最大特征值对应的特征向量
  % w2(:, k) = sqrt(max_val) * V(1:N_t, idx);
% w(:, k) = w_title(1:N_t) / w_title(N_t + 1);  % 归一化，使其最后一维为 1
%end 
%%%%%%%%%%%%%% L=32 %%%%%%%%%%%%%
%w3 = zeros(N_t, K);  % 初始化 w
%for k = 1:K
%    slice = W3(:, :, k);  % 取 W 的第 k 个 31x31 矩阵
 %   [V, D] = eig(slice);  % 计算特征值和特征向量
 %   [max_val, idx] = max(diag(D));  % 提取最大特征值
 %   w_title = sqrt(max_val) * V(:, idx);  % 取最大特征值对应的特征向量
 %  w3(:, k) = sqrt(max_val) * V(1:N_t, idx);
% w(:, k) = w_title(1:N_t) / w_title(N_t + 1);  % 归一化，使其最后一维为 1
%end 

%%%%%%%%%%%%%%%%%% no ris %%%%%%%%%%%%%%%
%w_noRIS = zeros(N_t, K);  
%sigma_wnoRIS = zeros(N_t, N_t);  
%for k = 1:K
 %  slice_noRIS = WnoRIS(:, :, k);  
 %  [VnoRIS, DnoRIS] = eig(slice_noRIS);  
  %  [max_val, idx] = max(diag(DnoRIS));  
  %  w_noRIS(:, k) = sqrt(max_val) * V(1:N_t, idx);  
  %  sigma_wnoRIS =  w_noRIS(:, k) * w_noRIS(:, k)';  
%end
%sigma_w2 = w2 * w2';
%sigma_w3 = w3 * w3';
%%%%%%%%%%% 画图 %%%%%%%%%%%%%%%%
%plot_beampattern(N_t, w,w_noRIS,R)
%plot_beampattern_L(N_t, w,w3,w2)
%plot_CDF(N_t, K, hB, hR, G, w, sigma_k2,lambda_c,lambda_s,a_s,a_c,alpha_s,alpha_c,C,sigma_s2);
%%%%%%% 测算beampattern_error,Tx_power,R_k,R_s %%%%%%%%%%%
result(sigma_w ,R,h_k,w,sigma_k2,alpha_s,a_s,alpha_c,a_c,sigma_s2,K,C);
disp(norm(sigma_w  - sum(W(1:N_t, 1:N_t), 3), 'fro'));
%%%%%%%%% 求解v %%%%%%%%%%%%%
v = solvev(N_t,L,G,hR,hB,w,K,f_k,sigma_k2,lambda_c);

%%%%%%%%% 求解f %%%%%%%%%%%%%

f = solvef(G,w,hB,hR,v,sigma_k2,K);

%%%% 计算 obj_curr（目标函数）
%%%%%%%%%% 计算 R_k 
R_k = zeros(K, 1);
SINR = zeros(K, 1);
h_k(:,k) = hB(:,k)' + v.' * diag(hR(:,k)) * G;

for k = 1:K
    signal_power = abs(h_k(:,k)' * w(:,k))^2;
    interference_power = 0;
    for i = 1:K
        if i ~= k
            interference_power = interference_power + abs(h_k(:,k)' * w(:,i))^2;
        end
    end
    SINR(k) = signal_power / (interference_power + sigma_k2); % 信干噪比
    R_k(k) = log2(1 + SINR(k));
end
%%%%%%%%%% 计算 R_s 
denominator4 = 0;
numerator4 = 0;

for k = 1:K
    numerator4 = numerator4 + abs(alpha_s * (a_s)' * w(:,k))^2; % 计算分子
    for c = 1:C
        denominator4 = denominator4 + abs(alpha_c * a_c(:,c)' * w(:,k))^2;
    end
end
SCNR = numerator4 / (denominator4 + sigma_s2); % 计算信干噪比
R_s = log2(1 + SCNR);
    
    % 计算目标函数
    obj_curr = lambda_c * sum(R_k) + lambda_s * R_s;

% 记录当前迭代的 WSR 值
    WSR_values(iter) = obj_curr; % 之后的迭代正常计算

        % 打印当前迭代信息
    fprintf('Iteration %d: WSR = %.4f\n', iter, obj_curr);

    % 终止条件检查
    if abs(obj_curr - obj_prev) <= tau_3 
        fprintf('Converged at iteration %d\n', iter);
        break;
    end

    % 更新 obj_prev
    obj_prev = obj_curr;

    % 迭代次数 +1
    iter = iter + 1;
end
% 迭代次数数组
%iterations = 1:(length(WSR_values));
% 绘制加权和速率（WSR）随迭代次数的变化
%figure;
%plot(iterations, WSR_values(1:length(iterations)), '-o', 'LineWidth', 2, 'MarkerSize', 6, 'Color', 'b');
%grid on;
%xlabel('iterations');
%ylabel('Weighted Sum Rate（bps/Hz）');
%set(gca, 'FontSize', 12);
%legend('FP-SDP-SOCP', 'Location', 'southeast');
%xlim([1 length(iterations)]);
%ylim([min(WSR_values)-0.1, max(WSR_values)+0.1]);



 








