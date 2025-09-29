
% 这个示例演示如何使用heatfd_Crank_Nicolson函数求解1D热方程

clear; clc; close all;

function_path = 'C:\Users\yuki\Documents\MATLAB\heat-equation-crank-nicolson';
addpath(function_path); % 将项目根目录添加到路径

fprintf('=== 1D热方程Crank-Nicolson求解器示例 ===\n\n');

% 示例1: 使用默认参数
fprintf('1. 使用默认参数运行...\n');
[w1, errors1, rate1] = heatfd_Crank_Nicolson();
fprintf('默认参数收敛阶: %.3f\n\n', rate1);

% 示例2: 自定义参数 - 更密的网格
fprintf('2. 使用更密的网格运行...\n');
args.M=100;
args.N=2000;

[w2, errors2, rate2] = heatfd_Crank_Nicolson(args);
fprintf('密网格收敛阶: %.3f\n\n', rate2);


fprintf('\n示例运行完成！\n');
fprintf('查看生成的图表来观察数值解的行为。\n');