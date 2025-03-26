function [R] = solveR(N_t, Phi_s)

f = 1; % 归一化频率 (单位: 波长)
P0 = 1; % 功率预算

% 定义副瓣范围
Omega = [-pi/2, -15 * pi/180, 15 * pi/180, pi/2]; 

% 角度点数
M = 200; 

% 只取副瓣范围内的角度
phi_m = [linspace(Omega(1), Omega(2), M/2), linspace(Omega(3), Omega(4), M/2)];%%%sidelobe

% 计算 steering vector
 alpha_0= steering_vector(Phi_s, N_t);
 alpha_1= steering_vector(Phi_s+pi*15/180, N_t);
 alpha_2= steering_vector(Phi_s-pi*15/180, N_t);
 alpha_m = zeros(N_t,M);
 for k = 1:M
 alpha_m(:,k) = steering_vector(phi_m(k), N_t);
 end
% CVX optimization
cvx_begin sdp quiet
    variable R(N_t, N_t) hermitian semidefinite
    variable t

    % 目标函数
    minimize(-t)

    % 约束 (10a) 
for k = 1:M
        alpha_0' * R * alpha_0 - alpha_m(:,k)' * R * alpha_m(:,k) >= t;
end

% 约束 (10b) 

alpha_1' * R * alpha_1 <= alpha_0' * R * alpha_0 / 2;

% 约束 (10c) 
alpha_2' * R * alpha_2 <= alpha_0' * R * alpha_0 / 2;

% 约束 (10d) R 的特性
R >= 0;
    % 约束 (10e) 
diag(R) == P0*ones(N_t,1) / N_t;
cvx_end
