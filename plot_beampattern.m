function plot_beampattern(N_t, w,w_noRIS, R)
    % 构建角度范围
    phi_range = linspace(-90, 90, 181) * pi / 180; % 角度转换为弧度
    G_phi_w = zeros(size(phi_range)); % 用 w*w' 计算
    G_phi_R = zeros(size(phi_range)); % 用 R 计算
   G_phi_noRIS = zeros(size(phi_range)); % 无 RIS 情况




    % 计算 Beampattern Gain
    for i = 1:length(phi_range)
        a_phi = exp(1j * 2 * pi * (0:N_t-1)' * (1/2) * sin(phi_range(i))); % 导向矢量
        G_phi_w(i) = real(a_phi' * (w * w') * a_phi); % 通信波束
        G_phi_R(i) = real(a_phi' * R * a_phi); % 雷达波束
        G_phi_noRIS(i) = real(a_phi' * (w_noRIS * w_noRIS') * a_phi); % 无 RIS 波束
    end


% 设置 Beampattern Gain dB 范围
G_phi_noRIS_dB = 10 * log10(max(G_phi_noRIS, 1e-10)); % 避免 log(0)
G_phi_w_dB = 10 * log10(max(G_phi_w, 1e-10));  
G_phi_R_dB = 10 * log10(max(G_phi_R, 1e-10));

% 绘制 Beampattern
figure;
hold on;
plot(-90:90, G_phi_w_dB, 'b', 'LineWidth', 2);  % 绿色曲线
plot(-90:90, G_phi_R_dB, 'r--', 'LineWidth', 2); % 红色虚线
plot(-90:90, G_phi_noRIS_dB, 'k-.', 'LineWidth', 2); % 无 RIS 绿色点虚线
hold off;

% 设置图例
legend('RIS Communications Beamforming', 'Radar-Only Beamforming', 'No RIS Beamforming', 'Location', 'northwest');
%legend('RIS Communications Beamforming', 'Radar-Only Beamforming', 'Location', 'northwest');
% 轴标签
xlabel('Angle (degree)');
ylabel('Beampattern Gain (dB)');

% 调整坐标范围
xlim([-90, 90]);
ylim([-50, 30]); % 限制 dB 范围，避免极端值影响可视化

grid on;