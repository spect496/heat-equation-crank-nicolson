
function [w, errors, convergence_rate] = heatfd_Crank_Nicolson(varargin)
% 实现热方程Crank-Nicolson解法, 可视化，分析误差和收敛性
% 热方程: ut = D * uxx xl <= x <= xr, t0 <= t <= t1, D>0
% 初始条件: u(x, t0) = f(x)
% 边界条件: u(xl, t) = g1(t), u(xr, t)=g2(t)
% 该程序中，f(x)=sin(pi*x)，边界温度恒为0
%
% 输入参数:
% 空间边界xl, xr; 时间边界t0, t1; 空间步数M，时间步数N, 扩散系数D
% 输出：解矩阵w, size(w) = M+1 N+1, 可视化图表, 误差分析，收敛分析。

% 处理输入参数
if nargin == 1 && isstruct(varargin{1})
    % 如果输入是结构体，直接使用
    args = varargin{1};
else
    % 否则使用默认参数
    args = struct();
end

% 设置默认参数
if ~isfield(args, 'xl'), args.xl = 0; end
if ~isfield(args, 'xr'), args.xr = 1; end
if ~isfield(args, 't0'), args.t0 = 0; end
if ~isfield(args, 't1'), args.t1 = 1; end
if ~isfield(args, 'M'), args.M = 50; end
if ~isfield(args, 'N'), args.N = 1000; end
if ~isfield(args, 'D'), args.D = 0.01; end


% 创建网格
x = linspace(args.xl, args.xr, args.M+1);
t = linspace(args.t0, args.t1, args.N+1);

% 主计算
w = compute_cn_solution(args);

% 可视化
visualize_results(x, t, w, args.M, args.N);

% 误差分析和收敛性分析
errors = error_analysis(x, t, w, args.D);
convergence_rate = analyze_convergence(args);
end

function visualize_results(x, t, w, M, N)
% 可视化CN方法的解

figure('Position', [200, 200, 1200, 800]);

% 1. 3D图
subplot(2,2,1);
[X, T] = meshgrid(x, t);
surf(X, T, w', 'EdgeColor', 'none');
xlabel('时间 t');
ylabel('位置 x');
zlabel('温度 u(x,t)');
title('Crank-Nicolson方法求解热方程');
colorbar;
view(45, 30);

% 2. 不同时刻的温度分布
subplot(2,2,2);
time_indices = [1, round(N/4), round(N/2), round(3*N/4), N+1];
colors = ['r', 'g', 'b', 'm', 'k'];
legends = {};
hold on;
for i = 1:length(time_indices)
    idx = time_indices(i);
    plot(x, w(:,idx), colors(i), 'LineWidth', 1.5);
    legends{i} = sprintf('t = %.3f', t(idx));
end
hold off;
xlabel('位置 x');
ylabel('温度 u(x,t)');
title('不同时刻的温度分布');
legend(legends, 'Location', 'best');
grid on;

% 3. 热图
subplot(2,2,3);
imagesc(x, t, w);
xlabel('时间 t');
ylabel('位置 x');
title('温度分布热图');
colorbar;
axis xy;

% 4. 动画
subplot(2,2,4);
h_plot = plot(x, w(:,1), 'b-', 'LineWidth', 2);
xlabel('位置 x');
ylabel('温度 u(x,t)');
title(sprintf('时间 t = %.3f', t(1)));
ylim([min(w(:)), max(w(:))]);
grid on;

% 动画循环
for j = 1:10:length(t)
    set(h_plot, 'YData', w(:,j));
    title(sprintf('时间 t = %.3f', t(j)));
    drawnow;
    pause(0.01);
end
end

function errors = error_analysis(x, t, w, D)

% 计算最终时刻数值解和精确解之间的误差
errors = struct();

w_exact = exact_solution(x, t(end), D);
w_cn = w(:, end)';

% L2误差（均方误差）
errors.L2 = sqrt(trapz(x, (w_cn - w_exact).^2));

% 绘制误差分布图
figure('Position', [300, 300, 800, 600]);
subplot(2,1,1);
plot(x, w_cn, 'r-', 'LineWidth', 2); hold on;
plot(x, w_exact, 'b--', 'LineWidth', 1.5);
xlabel('位置 x');
ylabel('温度 u(x,t)');
title(sprintf('最终时刻CN解与精确解对比'));
legend('CN解', '精确解', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(x, abs(w_cn - w_exact), 'k-', 'LineWidth', 2);
xlabel('位置 x');
ylabel('绝对误差');
title(sprintf('绝对误差分布 (L2误差 = %.2e)', errors.L2));
grid on;
end

function convergence_rate = analyze_convergence(args)
% 收敛性分析函数, 验证是否2阶收敛, 绘制关于h的收敛图

refinement_times = 3; % 网格细化次数
refinement_ratios = zeros(refinement_times, 1);

errors_L2 = zeros(refinement_times, 1);
h_values = zeros(refinement_times, 1);

% 计算不同网格密度下的误差
for i = 1 : refinement_times
    ratio = 2 ^ -i; % 2倍细化
    refinement_ratios(i) = ratio; 

    % 第i次细化的参数
    current_args = args;
    current_args.M = args.M / ratio; % 空间，时间等比例细化
    current_args.N = args.N / ratio;

    h = (args.xr - args.xl) / current_args.M; 
    h_values(i) = h;

    w_current = compute_cn_solution(current_args); % 计算网格 
    x_current = linspace(args.xl, args.xr, current_args.M+1);
    w_exact = exact_solution(x_current, args.t1, args.D);

    % 计算误差
    errors_L2(i) = sqrt(trapz(x_current, (w_current(:, end) - w_exact').^2));
end

% 计算收敛阶
convergence_orders = zeros(refinement_times-1, 1);
for i = 1 : refinement_times-1
    convergence_orders(i) = log(errors_L2(i)/errors_L2(i+1)) / log(refinement_ratios(i)/refinement_ratios(i+1));
end
convergence_rate = mean(convergence_orders);

% 绘制收敛图(关于h)
figure;
loglog(h_values, errors_L2, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;

% 理论收敛线
ref_slope = 2; % 2阶
x_ref = [min(h_values), max(h_values)];
y_ref = errors_L2(1) * (x_ref / h_values(1)).^ref_slope;
loglog(x_ref, y_ref, 'k--', 'DisplayName', '理论O(h^2)');

xlabel('空间步长 h');
ylabel('L2误差');
title(sprintf('收敛性分析 (平均收敛阶 = %.3f)', convergence_rate));
legend('数值误差', '理论O(h^2)', 'Location', 'best');
grid on;
end

function w = compute_cn_solution(args)
    % Crank-Nicolson解法主要计算过程

    % 创建网格
    x = linspace(args.xl, args.xr, args.M+1);
    t = linspace(args.t0, args.t1, args.N+1);

    w = zeros(args.M+1, args.N+1); % 初始化解矩阵

    % 初始条件
    w(:,1) = sin(pi * x); 

    % 边界条件
    w(1,:)=0;
    w(end,:)=0;

    h = (args.xr - args.xl) / args.M; % 空间步长
    k = (args.t1 - args.t0) / args.N; % 时间步长
    r = args.D * k / (h^2); % 稳定性参数

    % 构造系数矩阵

    % 内部格点数
    m = args.M - 1; 
    n = args.N;

    % 构造三对角矩阵
    a = diag((2+2*r)*ones(m,1)) + diag(-r*ones(m-1,1), 1); 
    a = a + diag(-r*ones(m-1,1), -1);
    b = diag((2-2*r)*ones(m,1)) + diag(r*ones(m-1,1), 1); 
    b = b + diag(r*ones(m-1,1), -1);

    % 内部求解
    for j = 1:n
        rhs = b * w(2:end-1,j) + r * [w(1,j) + w(1,j+1); zeros(m-2,1); 
            w(end,j) + w(end,j+1)];
        w(2:end-1,j+1) = a \ rhs;
    end
end

function u_exact = exact_solution(x, t, D)
% 精确解函数
u_exact = sin(pi * x) * exp(-D * pi^2 * t);
end



