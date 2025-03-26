function plot_iterations_WSR(N_t, K, hB, hR, G, w, sigma_k2, alpha_s, alpha_c, a_c, a_s, C, sigma_s2, lambda_c, lambda_s)
    load('v_storage.mat', 'v_storage'); % 加载 v_storage

    num_samples = size(v_storage, 2); % 获取总样本数
    num_iterations = 30; % 选取 30 组 v 进行计算

    % 确保 v_storage 至少有 30 组数据
    if num_samples < num_iterations
        error('v_storage 数据不足 30 组，请检查 v_storage.mat');
    end

    % 随机选取 30 组 v（如果想要等间隔选取，也可以修改）
    selected_indices = linspace(1, num_samples, num_iterations);
    selected_indices = round(selected_indices); % 取整保证是整数索引

    % 初始化存储 WSR
    WSR_all = zeros(num_iterations, 1);

    % 遍历选定的 30 组 v 计算 WSR
    for idx = 1:num_iterations
        v = v_storage(:, selected_indices(idx));

        % 计算信道系数 h_k
        h_k = zeros(N_t, K);
        for k = 1:K
            h_k(:, k) = hB(:, k).' + v.' * diag(hR(:, k)) * G;
        end

        % 计算 SINR 和 Sum Rate
        SINR = zeros(K, 1);
         R_k = zeros(K, 1);
        for k = 1:K
            signal_power = abs(h_k(:, k)' * w(:, k))^2;
            interference_power = 0;
            for i = 1:K
                if i ~= k
                    interference_power = interference_power + abs(h_k(:, k)' * w(:, i))^2;
                end
            end
            SINR(k) = signal_power / (interference_power + sigma_k2);
            R_k(k) = log2(1 + SINR(k)); % 用户 k 的速率 R_k
        end

        % 计算 R_s (感知速率)
        numerator4 = 0;
        denominator4 = 0;
        for k = 1:K
            numerator4 = numerator4 + abs(alpha_s * a_s' * w(:, k))^2; % 计算分子部分
        end
        for c = 1:C
            denominator4 = denominator4 + abs(alpha_c * a_c(:, c)' * w(:, k))^2; % 计算干扰部分
        end
        SCNR = numerator4 / (denominator4 + sigma_s2); % 计算 SCNR
        R_s = log2(1 + SCNR); % 计算数据速率 R_s

        % 计算加权和速率 (WSR)
        WSR_all(idx) = lambda_c * sum( R_k) + lambda_s * R_s;
    end

    % 绘制 WSR 曲线（30 次迭代）
    figure;
    plot(1:num_iterations, WSR_all, '-o', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Iteration Number', 'FontSize', 12);
    ylabel('WSR (Weighted Sum Rate)', 'FontSize', 12);
    title('WSR vs. Iteration Number (30 Iterations)', 'FontSize', 14);
    grid on;
    xlim([1, num_iterations]);
    ylim([min(WSR_all)-1, max(WSR_all)+1]);
    legend('WSR Curve', 'Location', 'best');
end
