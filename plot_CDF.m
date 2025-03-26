function plot_CDF(N_t, K, hB, hR, G, w, sigma_k2,lambda_c,lambda_s,a_s,a_c,alpha_s,alpha_c,C,sigma_s2)
    load('v_storage.mat','v_storage');

    num_samples = size(v_storage,2);



% 初始化存储1000组 v 值的速率总和
    Rate_sum_all = zeros(num_samples, 1);

    % 遍历 13 组 v 值
    for idx = 1:num_samples
        v = v_storage(:,idx);

        % 计算信道系数 h_k
        h_k = zeros(N_t, K);
        for k = 1:K
            h_k(:,k) = hB(:,k).' + v.' * diag(hR(:,k)) * G;
        end

        % 计算 SINR 和 Rate
        SINR_RIS1 = zeros(K, 1);
        Rate_RIS1 = zeros(K, 1);
        interference_power = 0;
        for k = 1:K
            signal_power = abs(h_k(:,k)' * w(:,k))^2;
            for i = 1:K
             if i ~= k
            interference_power = interference_power + abs(h_k(:,k)' * w(:,i)).^2 ;
             end
            end
            SINR_RIS1(k) = signal_power / (interference_power + sigma_k2);
            Rate_RIS1(k) = log2(1 + SINR_RIS1(k));
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
        % 记录该组 v 值的速率总和
        Rate_sum_all(idx) = sum(Rate_RIS1)+R_s;  
    end


    % **绘制仿真 CDF**
    figure;
    cdfplot(Rate_sum_all);
    hold on;

% **调整理论 CDF 计算范围**
x = linspace(min(Rate_sum_all), max(Rate_sum_all), 100);  % 确保 x 从 Rate_sum_all 的最小值开始
mu = mean(Rate_sum_all);
sigma = std(Rate_sum_all);

% **使用正态分布 CDF 拟合**
y_theory = normcdf(x, mu, sigma);  % 计算理论 CDF

    % **绘制理论 CDF**
    plot(x, y_theory, 'r', 'LineWidth', 0.5);

    grid on;
    xlabel('Weighted Sum Rate (bits/s/Hz)');
    ylabel('CDF');
    legend(' FP-SDP-SOCP CDF', 'Theoretical CDF', 'Location', 'best');
    hold off;
end
