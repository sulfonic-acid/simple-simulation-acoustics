% 两端固定弦的非对称振动 - 振型叠加动画
clear; clc; close all;

% 参数设置
L = 1;          % 弦长
c = 1;          % 波速
N_modes = 20;   % 考虑的模态数
duration = 5;   % 动画时长
fps = 25;       % 帧率

% 非对称激励位置（不在中点）
excitation_pos = 0.3 * L;

% 空间离散化
x = linspace(0, L, 200);

% 时间离散化
t = linspace(0, duration, duration * fps);

% 初始化图形窗口
figure('Position', [100, 100, 1200, 800]);

% 计算各阶模态的权重（非对称激励）
weights = zeros(1, N_modes);
for n = 1:N_modes
    % 模态形状函数（正弦函数）
    phi_n = sin(n * pi * excitation_pos / L);
    % 权重与激励位置的模态值成正比
    weights(n) = phi_n / n^2;  % n^2 使得高阶模态振幅衰减
end

% 归一化权重
weights = weights / max(abs(weights));

% 创建GIF文件
filename = 'string_vibration.gif';

for frame = 1:length(t)
    clf;
    
    % 初始化总位移和各阶模态位移
    total_displacement = zeros(size(x));
    mode_displacements = zeros(N_modes, length(x));
    
    % 计算各阶模态的响应
    for n = 1:N_modes
        % 第n阶模态的频率
        omega_n = n * pi * c / L;
        
        % 第n阶模态的形状函数
        phi_n = sin(n * pi * x / L);
        
        % 第n阶模态的时间函数（简谐振动）
        time_factor = cos(omega_n * t(frame));
        
        % 第n阶模态的位移
        mode_displacement = weights(n) * phi_n * time_factor;
        mode_displacements(n, :) = mode_displacement;
        
        % 叠加到总位移
        total_displacement = total_displacement + mode_displacement;
    end
    
    % 绘制总振动
    subplot(2, 2, [1, 2]);
    plot(x, total_displacement, 'b-', 'LineWidth', 3);
    hold on;
    plot([0, L], [0, 0], 'k--', 'LineWidth', 1);
    plot(excitation_pos, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    xlabel('位置 x', 'FontSize', 12);
    ylabel('位移 u(x,t)', 'FontSize', 12);
    title('两端固定弦的非对称振动 - 总响应', 'FontSize', 14);
    grid on;
    axis([0, L, -1.2, 1.2]);
    legend('总位移', '平衡位置', '激励位置', 'Location', 'northeast');
    
    % 绘制前几阶主要模态
    subplot(2, 2, 3);
    colors = {'r', 'g', 'b', 'm', 'c'};
    for n = 1:min(5, N_modes)
        plot(x, mode_displacements(n, :), colors{n}, 'LineWidth', 1.5);
        hold on;
    end
    xlabel('位置 x', 'FontSize', 12);
    ylabel('位移', 'FontSize', 12);
    title('主要模态分量', 'FontSize', 12);
    grid on;
    axis([0, L, -0.8, 0.8]);
    legend('1阶', '2阶', '3阶', '4阶', '5阶', 'Location', 'northeast');
    
    % 绘制模态权重
    subplot(2, 2, 4);
    stem(1:N_modes, abs(weights), 'filled', 'LineWidth', 1.5);
    xlabel('模态阶数 n', 'FontSize', 12);
    ylabel('权重 |a_n|', 'FontSize', 12);
    title('各阶模态的权重分布', 'FontSize', 12);
    grid on;
    xlim([0.5, N_modes + 0.5]);
    
    % 添加时间信息
    sgtitle(sprintf('非对称振动 (激励位置: %.1fL, 时间: %.2f s)', ...
           excitation_pos/L, t(frame)), 'FontSize', 16);
    
    % 捕获帧并写入GIF
    drawnow;
    frame_img = getframe(gcf);
    im = frame2im(frame_img);
    [imind, cm] = rgb2ind(im, 256);
    
    if frame == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 1/fps);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
    end
end

% 显示模态信息
fprintf('非对称振动参数:\n');
fprintf('弦长: L = %.1f\n', L);
fprintf('激励位置: x = %.2fL\n', excitation_pos/L);
fprintf('考虑的模态数: %d\n', N_modes);
fprintf('前5阶模态权重:\n');
for n = 1:min(5, N_modes)
    fprintf('  第%d阶: %.4f\n', n, weights(n));
end

fprintf('\nGIF动画已保存为: %s\n', filename);

% 额外分析：显示非对称性特征
fprintf('\n非对称性分析:\n');
even_modes = weights(2:2:end);
odd_modes = weights(1:2:end);
fprintf('偶阶模态平均权重: %.4f\n', mean(abs(even_modes)));
fprintf('奇阶模态平均权重: %.4f\n', mean(abs(odd_modes)));
fprintf('非对称系数: %.4f\n', mean(abs(even_modes))/mean(abs(odd_modes)));