function [beampattern_error,Tx_power,R_k,SINR,R_s,SCNR]  = result(sigma_w ,R,h_k,w,sigma_k2,alpha_s,a_s,alpha_c,a_c,sigma_s2,K,C)

beampattern_error = norm( sigma_w  - R, 'fro')^2;
Tx_power = trace(sigma_w );

%%%%%%%%%% 计算 R_k 

R_k = zeros(K, 1);
SINR = zeros(K, 1);
SINR_dB = zeros(K, 1);
for k = 1:K
    signal_power = abs(h_k(:,k)' * w(:,k))^2; % 期望信号功率
    interference_power = 0; % 初始化干扰功率
    for i = 1:K
        if i ~= k
            interference_power = interference_power + abs(h_k(:,k)' * w(:,i))^2; % 累加干扰
        end
    end
    SINR(k) = signal_power / (interference_power + sigma_k2); % 计算 SINR
    SINR_dB(k) = 10 * log10(SINR(k));
    R_k(k) = log2(1 + SINR(k));
end
R_k(k) = log2(1 + SINR(k));
%%%%%%%%%% 计算 R_s 
denominator4 = 0;
for k = 1:K
numerator4 = abs(alpha_s * (a_s)' * w(:,k))^2;% 计算分子部分
end
for c = 1:C
denominator4 = denominator4 + abs(alpha_c * a_c(:,c)' * w(:,k))^2;
end
SCNR = numerator4 / (denominator4 + sigma_s2); % 加上噪声功率
R_s = log2(1 + SCNR); % 计算数据速率 R_s

disp(['beampattern_error is ',num2str(beampattern_error)]);
disp(['Tx_power is ',num2str(Tx_power)]);
disp(['R_k is ',num2str(sum(R_k))]);
disp(['R_s is ',num2str(R_s)]);
for i = 1:length(SINR)
    disp(['SINR', num2str(i), ' is ', num2str(SINR_dB(i))]);
end
