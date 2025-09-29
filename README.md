\# 一维热方程求解器（Crank-Nicolson方法）

用MATLAB实现一维热方程的数值求解，使用了Crank-Nicolson方法。这个方法是无条件稳定的，在时间和空间上都是二阶精度。

\## 文件结构

-heat-equation-crank-nicolson/
- heatfd_Crank_Nicolson.m                 (主程序文件)
- examples/                   (示例脚本)
   - example.m           (使用示例)
- images/                     (结果图片)
- convergence_plot.png    (收敛性分析图)
- error_analysis.png      (误差分析图)
- solution_visualization.png (整体可视化图)
- README.md                   (项目说明（就是这个文件）)

\## 背景

解决下面的热方程
```
∂u/∂t = D · ∂²u/∂x²,  x ∈ [xl, xr], t ∈ [t0, t1]
```
初始条件： `u(x, t0) = sin(πx)` 和 Dirichlet 边界条件 `u(xl, t) = u(xr, t) = 0`.


\## 我能做什么？


\- \*\*求解热方程\*\*：计算一根杆子上的温度随时间的变化

\- \*\*多种可视化\*\*：提供3D图、热力图、动画等多种方式看结果

\- \*\*误差分析\*\*：比较数值解和精确解的误差

\- \*\*收敛性验证\*\*：验证方法确实达到二阶收敛



\## 怎么使用？


在MATLAB里运行：

```matlab

% 最简单的方式（使用默认参数）

\[w, errors, convergence\_rate] = heatfd\_cn();

% 或者自定义参数

params.M = 100;

params.N = 2000;

\[w, errors, convergence\_rate] = heatfd\_cn(params);

```



\## 参数说明



\-  M: 空间网格数（默认50）

\-  N: 时间步数（默认1000）

\-  D: 扩散系数（默认0.01）

\-  xl，xr: 空间边界（默认0, 1）

\-  t0, t1: 时间边界（默认0, 1）



\## 输出结果



运行后会生成：

1\. 3D温度演化图

2\. 不同时刻的温度分布

3\. 热力图

4\. 温度变化动画

5\. 与精确解的对比图

6\. 误差分布图

7\. 收敛性分析图

\## 代码结构

heatfd\_cn.m - 主文件，包含所有功能

compute\_cn\_solution - 核心求解器

visualize\_results - 画图功能

error\_analysis - 误差分析

analyze\_convergence - 收敛性分析

exact\_solution - 精确解函数

