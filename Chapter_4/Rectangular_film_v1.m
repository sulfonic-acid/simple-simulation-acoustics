%% 矩形薄膜振动模态分析主程序
clear; close all; clc;

%% 参数设置
a = 1.0;        % x方向长度(m)
b = 1.0;        % y方向长度(m)
T = 1000;       % 张力(N/m)
rho = 1.0;      % 面密度(kg/m²)
Nx = 50;        % x方向网格点数
Ny = 50;        % y方向网格点数

% 外力参数
F0 = 100;       % 外力幅值(N)
omega_f = 50;   % 外力频率(rad/s)
x_force = 0.3;  % 外力作用点x坐标
y_force = 0.7;  % 外力作用点y坐标

% 时间参数
t_end = 2;      % 总时间(s)
dt = 0.01;      % 时间步长(s)

%% 创建网格
x = linspace(0, a, Nx);
y = linspace(0, b, Ny);
[X, Y] = meshgrid(x, y);
dx = a/(Nx-1);
dy = b/(Ny-1);

%% 1. 自由振动模态分析
fprintf('=== 自由振动模态分析 ===\n');

% 计算固有频率和模态形状
modes_to_show = [1,1; 1,2; 2,1; 2,2; 1,3; 3,1];
num_modes = size(modes_to_show, 1);

figure('Position', [100, 100, 1200, 800]);
for i = 1:num_modes
    m = modes_to_show(i,1);
    n = modes_to_show(i,2);
    
    % 理论固有频率
    omega_mn = pi * sqrt(T/rho) * sqrt((m/a)^2 + (n/b)^2);
    freq_mn = omega_mn / (2*pi);
    
    % 模态形状
    phi_mn = sin(m*pi*X/a) .* sin(n*pi*Y/b);
    
    % 绘制模态
    subplot(2, 3, i);
    surf(X, Y, phi_mn, 'EdgeColor', 'none');
    title(sprintf('模态(%d,%d)\nf=%.2f Hz', m, n, freq_mn));
    xlabel('x (m)'); ylabel('y (m)'); zlabel('振幅');
    axis equal;
    colormap(jet);
    colorbar;
end
sgtitle('矩形薄膜自由振动模态');

%% 2. 受迫振动分析 - 模态叠加法
fprintf('\n=== 受迫振动分析 (模态叠加法) ===\n');

% 选择参与计算的模态数量
M_max = 8;  % x方向最大模态数
N_max = 8;  % y方向最大模态数

% 初始化模态参数
omega_mn = zeros(M_max, N_max);
phi_mn = cell(M_max, N_max);
q_mn = zeros(M_max, N_max);  % 模态坐标

% 计算各模态参数
fprintf('计算模态参数...\n');
for m = 1:M_max
    for n = 1:N_max
        % 固有频率
        omega_mn(m,n) = pi * sqrt(T/rho) * sqrt((m/a)^2 + (n/b)^2);
        
        % 模态形状函数
        phi_mn{m,n} = sin(m*pi*X/a) .* sin(n*pi*Y/b);
        
        % 模态质量 (假设为单位质量)
        M_mn = 1.0;
        
        % 模态力 (外力在模态上的投影)
        phi_force = sin(m*pi*x_force/a) * sin(n*pi*y_force/b);
        F_mn = F0 * phi_force;
        
        % 稳态振幅 (忽略阻尼)
        if abs(omega_f^2 - omega_mn(m,n)^2) > 1e-6
            q_mn(m,n) = F_mn / (M_mn * (omega_mn(m,n)^2 - omega_f^2));
        else
            q_mn(m,n) = 0;  % 共振情况，需要阻尼模型
        end
    end
end

% 叠加各模态得到总响应
w_forced = zeros(size(X));
for m = 1:M_max
    for n = 1:N_max
        w_forced = w_forced + q_mn(m,n) * phi_mn{m,n};
    end
end

% 绘制受迫振动响应
figure('Position', [100, 100, 1000, 400]);

subplot(1, 2, 1);
surf(X, Y, w_forced, 'EdgeColor', 'none');
title('受迫振动稳态响应幅值');
xlabel('x (m)'); ylabel('y (m)'); zlabel('振幅 (m)');
axis equal;
colormap(jet);
colorbar;

subplot(1, 2, 2);
contourf(X, Y, abs(w_forced), 20, 'LineColor', 'none');
hold on;
plot(x_force, y_force, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
title('响应幅值等高线 (红点为外力作用点)');
xlabel('x (m)'); ylabel('y (m)');
colorbar;
axis equal;

%% 3. 受迫振动时程分析
fprintf('\n=== 时程响应分析 ===\n');

% 选择观察点
x_obs = 0.7; y_obs = 0.4;

% 时间向量
t = 0:dt:t_end;

% 计算观察点的响应 (包含瞬态和稳态)
w_obs = zeros(size(t));
damping_ratio = 0.02;  % 阻尼比

for m = 1:min(5, M_max)  % 只计算前几阶模态以减少计算量
    for n = 1:min(5, N_max)
        phi_obs = sin(m*pi*x_obs/a) * sin(n*pi*y_obs/b);
        phi_force = sin(m*pi*x_force/a) * sin(n*pi*y_force/b);
        
        F_mn = F0 * phi_force;
        M_mn = 1.0;
        
        omega_d = omega_mn(m,n) * sqrt(1 - damping_ratio^2);  % 有阻尼固有频率
        
        % 稳态响应
        if abs(omega_f^2 - omega_mn(m,n)^2) > 1e-6
            H = 1 / (M_mn * (omega_mn(m,n)^2 - omega_f^2 + 1i*2*damping_ratio*omega_mn(m,n)*omega_f));
            w_steady = abs(H) * F_mn * phi_obs * sin(omega_f*t + angle(H));
        else
            w_steady = zeros(size(t));
        end
        
        % 瞬态响应 (初始条件为0)
        A = -real(H) * F_mn * phi_obs;
        B = (damping_ratio*omega_mn(m,n)*real(H) - omega_f*imag(H)) * F_mn * phi_obs / omega_d;
        
        w_transient = exp(-damping_ratio*omega_mn(m,n)*t) .* ...
                     (A*cos(omega_d*t) + B*sin(omega_d*t));
        
        w_obs = w_obs + w_steady + w_transient;
    end
end

% 绘制时程响应
figure('Position', [100, 100, 800, 600]);

subplot(2, 1, 1);
plot(t, w_obs, 'b-', 'LineWidth', 2);
xlabel('时间 (s)'); ylabel('位移 (m)');
title(sprintf('观察点 (%.1f, %.1f) 的位移时程响应', x_obs, y_obs));
grid on;

% 计算并显示各阶模态的贡献
subplot(2, 1, 2);
modal_contributions = zeros(1, num_modes);
for i = 1:num_modes
    m = modes_to_show(i,1);
    n = modes_to_show(i,2);
    phi_obs = sin(m*pi*x_obs/a) * sin(n*pi*y_obs/b);
    phi_force = sin(m*pi*x_force/a) * sin(n*pi*y_force/b);
    
    F_mn = F0 * phi_force;
    M_mn = 1.0;
    
    if abs(omega_f^2 - omega_mn(m,n)^2) > 1e-6
        q_mn_val = F_mn / (M_mn * (omega_mn(m,n)^2 - omega_f^2));
        modal_contributions(i) = abs(q_mn_val * phi_obs);
    else
        modal_contributions(i) = 0;
    end
end

bar(modal_contributions);
xlabel('模态阶数'); ylabel('贡献幅值');
title('各阶模态对响应的贡献');
set(gca, 'XTickLabel', compose('(%d,%d)', modes_to_show));

%% 4. 共振分析 - 扫描频率响应
fprintf('\n=== 频率响应分析 ===\n');

% 频率扫描范围
omega_range = linspace(0.1, 150, 200);
response_amp = zeros(size(omega_range));

for k = 1:length(omega_range)
    omega_current = omega_range(k);
    w_temp = 0;
    
    for m = 1:min(6, M_max)
        for n = 1:min(6, N_max)
            phi_obs = sin(m*pi*x_obs/a) * sin(n*pi*y_obs/b);
            phi_force = sin(m*pi*x_force/a) * sin(n*pi*y_force/b);
            
            F_mn = F0 * phi_force;
            M_mn = 1.0;
            
            % 考虑阻尼的频率响应函数
            H = 1 / (M_mn * (omega_mn(m,n)^2 - omega_current^2 + ...
                            1i*2*damping_ratio*omega_mn(m,n)*omega_current));
            w_temp = w_temp + abs(H * F_mn * phi_obs);
        end
    end
    
    response_amp(k) = w_temp;
end

% 绘制频率响应
figure('Position', [100, 100, 800, 400]);
plot(omega_range, response_amp, 'b-', 'LineWidth', 2);
xlabel('激励频率 (rad/s)'); ylabel('响应幅值 (m)');
title('频率响应函数');
grid on;

% 标记前几个固有频率
hold on;
for i = 1:min(6, num_modes)
    m = modes_to_show(i,1);
    n = modes_to_show(i,2);
    omega_natural = omega_mn(m,n);
    plot([omega_natural, omega_natural], [0, max(response_amp)], 'r--', 'LineWidth', 1);
    text(omega_natural, max(response_amp)*0.9, sprintf('(%d,%d)', m, n), ...
         'HorizontalAlignment', 'center');
end
legend('频率响应', '固有频率');

%% 5. 动画显示受迫振动
fprintf('\n=== 生成受迫振动动画 ===\n');

% 选择要动画显示的时间段
t_animate = 0:0.05:1;  % 前1秒
w_animate = zeros(Ny, Nx, length(t_animate));

% 计算动画帧
fprintf('计算动画帧...');
for idx = 1:length(t_animate)
    t_current = t_animate(idx);
    w_frame = zeros(size(X));
    
    for m = 1:min(5, M_max)
        for n = 1:min(5, N_max)
            phi_xy = sin(m*pi*X/a) .* sin(n*pi*Y/b);
            phi_force = sin(m*pi*x_force/a) * sin(n*pi*y_force/b);
            
            F_mn = F0 * phi_force;
            M_mn = 1.0;
            
            % 稳态响应
            if abs(omega_f^2 - omega_mn(m,n)^2) > 1e-6
                H = 1 / (M_mn * (omega_mn(m,n)^2 - omega_f^2 + 1i*2*damping_ratio*omega_mn(m,n)*omega_f));
                w_frame = w_frame + real(H * F_mn * phi_xy * exp(1i*omega_f*t_current));
            end
        end
    end
    
    w_animate(:,:,idx) = w_frame;
end
fprintf('完成\n');

% 创建动画
figure('Position', [100, 100, 800, 600]);
for idx = 1:length(t_animate)
    surf(X, Y, w_animate(:,:,idx), 'EdgeColor', 'none');
    title(sprintf('受迫振动响应 t=%.2fs', t_animate(idx)));
    xlabel('x (m)'); ylabel('y (m)'); zlabel('位移 (m)');
    zlim([-max(abs(w_animate(:))), max(abs(w_animate(:)))]);
    colormap(jet);
    colorbar;
    
    % 标记外力作用点
    hold on;
    plot3(x_force, y_force, 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
    
    drawnow;
end

%% 输出总结信息
fprintf('\n=== 分析总结 ===\n');
fprintf('薄膜尺寸: %.1f x %.1f m\n', a, b);
fprintf('张力: %.0f N/m, 面密度: %.1f kg/m²\n', T, rho);
fprintf('外力: %.0f N, 频率: %.1f rad/s, 作用点: (%.1f, %.1f)\n', ...
        F0, omega_f, x_force, y_force);
fprintf('前6阶固有频率:\n');
for i = 1:6
    m = modes_to_show(i,1);
    n = modes_to_show(i,2);
    fprintf('  模态(%d,%d): %.2f Hz\n', m, n, omega_mn(m,n)/(2*pi));
end

fprintf('\n分析完成！\n');