% 清除工作区变量和命令行窗口
clear;
clc;

load('lolp.mat');
load('lolpess.mat');
data_without_ress = lolp;
data_with_ress = lolpess;
% --- 绘制直方图 ---

figure; % 创建一个新的图形窗口

% 为两个直方图定义共同的 bin 边界，以确保正确比较
% 找到所有数据的最大值，以设置合适的 bin
max_val = max([data_without_ress; data_with_ress]);
bin_edges = 0:1:ceil(max_val + 5); % bin 从0开始，步长为1，直到最大值加一些缓冲

% 绘制 'Without RESS' 直方图
h1 = histogram(data_without_ress, 'BinEdges', bin_edges, 'Normalization', 'probability');
h1.FaceColor = [0.2 0.5 0.8]; % 'Without RESS' 的蓝色
h1.DisplayName = 'Without RESS'; % 图例名称
hold on; % 保持当前图，以便在同一坐标轴上绘制第二个直方图

% 绘制 'With RESS' 直方图
h2 = histogram(data_with_ress, 'BinEdges', bin_edges, 'Normalization', 'probability');
h2.FaceColor = [0.8 0.3 0.3]; % 'With RESS' 的红色
h2.DisplayName = 'With RESS'; % 图例名称

% --- 添加 NERC 限制线 ---
nerc_limit = 3; % 根据图像观察，大约是3小时
xline(nerc_limit, 'g-', 'LineWidth', 2, 'DisplayName', 'NERC limit'); % 绘制垂直线

% --- 自定义图表 ---
xlabel('每年负载损失小时数', 'FontSize', 12); % X轴标签
ylabel('概率', 'FontSize', 12); % Y轴标签
title('负载损失小时数概率分布', 'FontSize', 14); % 图表标题
grid on; % 添加网格线，提高可读性
legend('Location', 'NorthEast'); % 显示图例，位置在右上角
set(gca, 'XLim', [0 50]); % 设置X轴范围，与图像相似
set(gca, 'YLim', [0 1]); % 设置Y轴范围，与图像相似
set(gca, 'FontSize', 10); % 调整刻度标签字体大小