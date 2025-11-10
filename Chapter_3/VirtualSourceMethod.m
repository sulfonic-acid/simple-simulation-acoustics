%% 虚源方法演示：一维波在硬边界(x=0)的反射
% 课程用图 - 展示虚源法的物理意义
clear; clc; close all;

%% 参数设置
c = 1;          % 波速
L = 10;         % 显示区域长度
T_total = 8;    % 总时间
dt = 0.1;       % 时间步长
x = linspace(-L, L, 1000);  % 空间坐标

% 初始高斯波包（向右传播）
x0 = 3;         % 初始位置
sigma = 0.5;    % 波包宽度
k = 2;          % 波数

%% 创建图形
figure('Position', [100, 100, 800, 600]);

% 预计算所有帧
t_values = 0:dt:T_total;
num_frames = length(t_values);

% 预分配存储
u_incident = zeros(length(x), num_frames);    % 入射波
u_image = zeros(length(x), num_frames);       % 虚源波  
u_total = zeros(length(x), num_frames);       % 总波场

%% 计算波场演化
for n = 1:num_frames
    t = t_values(n);
    
    % 第一行：入射波（从右向左传播）
    u_incident(:, n) = exp(-((x - x0 + c*t)/sigma).^2) .* cos(k*(x - x0 + c*t));
    
    % 第二行：虚源波（从左向右传播）
    u_image(:, n) = -exp(-((x + x0 - c*t)/sigma).^2) .* cos(k*(x + x0 - c*t));
    
    % 第三行：总波场（x>0区域）
    u_total(:, n) = zeros(size(x));
    pos_indices = x >= 0;  % 物理区域
    u_total(pos_indices, n) = u_incident(pos_indices, n) + u_image(pos_indices, n);
end

%% 创建GIF动画
filename = 'virtual_source_method.gif';

for n = 1:num_frames
    clf;
    
    % 第一行：入射波传播
    subplot(3,1,1);
    plot(x, u_incident(:, n), 'b-', 'LineWidth', 2);
    hold on;
    % 标记边界
    plot([0, 0], [-1, 1], 'k--', 'LineWidth', 2);
    % 标记传播方向
    if n < num_frames/2
        annotation('textarrow', [0.7, 0.6], [0.85, 0.85], 'String', '传播方向');
    end
    ylim([-1.2, 1.2]);
    ylabel('振幅');
    title('第一行: 从右向左传播的入射波');
    grid on;
    
    % 第二行：虚源波传播  
    subplot(3,1,2);
    plot(x, u_image(:, n), 'r-', 'LineWidth', 2);
    hold on;
    plot([0, 0], [-1, 1], 'k--', 'LineWidth', 2);
    % 标记传播方向
    if n > num_frames/2
        annotation('textarrow', [0.3, 0.4], [0.55, 0.55], 'String', '传播方向');
    end
    ylim([-1.2, 1.2]);
    ylabel('振幅');
    title('第二行: 从左向右传播的虚源波');
    grid on;
    
    % 第三行：总波场（物理区域）
    subplot(3,1,3);
    plot(x, u_total(:, n), 'g-', 'LineWidth', 3);
    hold on;
    plot([0, 0], [-1, 1], 'k--', 'LineWidth', 2);
    % 突出显示物理区域
    fill([0, L, L, 0], [-1.2, -1.2, -1.1, -1.1], 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    text(L/2, -1.15, '物理区域 (x>0)', 'HorizontalAlignment', 'center');
    
    ylim([-1.2, 1.2]);
    xlabel('位置 x');
    ylabel('振幅');
    title('第三行: 实际观测到的总波场 (入射波 + 虚源波)');
    grid on;
    
    % 添加时间信息
    sgtitle(sprintf('虚源方法演示 - 时间 t = %.1f (反射系数 r = -1)', t_values(n)), ...
           'FontSize', 14, 'FontWeight', 'bold');
    
    % 捕获帧并写入GIF
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    if n == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end

%% 额外：绘制关键时刻的对比图
key_times = [1, 3, 5, 7]; % 关键时间点
figure('Position', [200, 200, 1000, 800]);

for i = 1:length(key_times)
    [~, idx] = min(abs(t_values - key_times(i)));
    t = t_values(idx);
    
    subplot(2, 2, i);
    
    % 绘制三个波场
    plot(x, u_incident(:, idx), 'b-', 'LineWidth', 1.5, 'DisplayName', '入射波');
    hold on;
    plot(x, u_image(:, idx), 'r-', 'LineWidth', 1.5, 'DisplayName', '虚源波');
    plot(x, u_total(:, idx), 'g-', 'LineWidth', 2.5, 'DisplayName', '总波场');
    
    % 标记边界和区域
    plot([0, 0], [-1, 1], 'k--', 'LineWidth', 2, 'DisplayName', '硬边界');
    fill([0, L, L, 0], [-1.2, -1.2, -1.1, -1.1], 'y', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    
    ylim([-1.2, 1.2]);
    xlabel('位置 x');
    ylabel('振幅');
    title(sprintf('时间 t = %.1f', t));
    legend('Location', 'best');
    grid on;
end

sgtitle('虚源方法关键时间点对比', 'FontSize', 16, 'FontWeight', 'bold');

fprintf('GIF动画已保存为: %s\n', filename);
fprintf('关键时间点对比图已生成\n');

%% 物理原理说明
fprintf('\n=== 虚源方法物理原理说明 ===\n');
fprintf('1. 硬边界条件 (r = -1): 位移在边界处必须为零\n');
fprintf('2. 虚源方法: 在边界另一侧放置一个镜像源，其相位相反\n');
fprintf('3. 物理区域 (x>0) 的总波场 = 实际入射波 + 虚源产生的波\n');
fprintf('4. 自动满足边界条件: u(0,t) = u_incident(0,t) + u_image(0,t) = 0\n');
