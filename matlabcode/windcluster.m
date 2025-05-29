% function [p,P] = windcluster(n)
clear;clc;
warning('off');n = 12;
%% 读取数据

Boise_2007 = readtable('E:\复现2\matlabcode\data2\Boise\Boise_2007.csv');
Boise_2008 = readtable('E:\复现2\matlabcode\data2\Boise\Boise_2008.csv');
Boise_2009 = readtable('E:\复现2\matlabcode\data2\Boise\Boise_2009.csv');
Boise_2010 = readtable('E:\复现2\matlabcode\data2\Boise\Boise_2010.csv');
Boise_2011 = readtable('E:\复现2\matlabcode\data2\Boise\Boise_2011.csv');
Boise_2012 = readtable('E:\复现2\matlabcode\data2\Boise\Boise_2012.csv');
Boise_2013 = readtable('E:\复现2\matlabcode\data2\Boise\Boise_2013.csv');
Boise_2014 = readtable('E:\复现2\matlabcode\data2\Boise\Boise_2014.csv');

Boise_all = [Boise_2007.windSpeedAt80m_m_s_;Boise_2008.windSpeedAt80m_m_s_;Boise_2009.windSpeedAt80m_m_s_;
    Boise_2010.windSpeedAt80m_m_s_;Boise_2011.windSpeedAt80m_m_s_;Boise_2012.windSpeedAt80m_m_s_;
    Boise_2013.windSpeedAt80m_m_s_;Boise_2014.windSpeedAt80m_m_s_];
% Boise_all = Boise_all';
CO_2007 = readtable('E:\复现2\matlabcode\data2\CO\CO_2007.csv');
CO_2008 = readtable('E:\复现2\matlabcode\data2\CO\CO_2008.csv');
CO_2009 = readtable('E:\复现2\matlabcode\data2\CO\CO_2009.csv');
CO_2010 = readtable('E:\复现2\matlabcode\data2\CO\CO_2010.csv');
CO_2011 = readtable('E:\复现2\matlabcode\data2\CO\CO_2011.csv');
CO_2012 = readtable('E:\复现2\matlabcode\data2\CO\CO_2012.csv');
CO_2013 = readtable('E:\复现2\matlabcode\data2\CO\CO_2013.csv');
CO_2014 = readtable('E:\复现2\matlabcode\data2\CO\CO_2014.csv');

CO_all = [CO_2007.windSpeedAt80m_m_s_;CO_2008.windSpeedAt80m_m_s_;CO_2009.windSpeedAt80m_m_s_;
    CO_2010.windSpeedAt80m_m_s_;CO_2011.windSpeedAt80m_m_s_;CO_2012.windSpeedAt80m_m_s_;
    CO_2013.windSpeedAt80m_m_s_;CO_2014.windSpeedAt80m_m_s_];
% CO_all = CO_all';

Fort_2007 = readtable('E:\复现2\matlabcode\data2\Fort\Fort_2007.csv');
Fort_2008 = readtable('E:\复现2\matlabcode\data2\Fort\Fort_2008.csv');
Fort_2009 = readtable('E:\复现2\matlabcode\data2\Fort\Fort_2009.csv');
Fort_2010 = readtable('E:\复现2\matlabcode\data2\Fort\Fort_2010.csv');
Fort_2011 = readtable('E:\复现2\matlabcode\data2\Fort\Fort_2011.csv');
Fort_2012 = readtable('E:\复现2\matlabcode\data2\Fort\Fort_2012.csv');
Fort_2013 = readtable('E:\复现2\matlabcode\data2\Fort\Fort_2013.csv');
Fort_2014 = readtable('E:\复现2\matlabcode\data2\Fort\Fort_2014.csv');

Fort_all = [Fort_2007.windSpeedAt80m_m_s_;Fort_2008.windSpeedAt80m_m_s_;Fort_2009.windSpeedAt80m_m_s_;
    Fort_2010.windSpeedAt80m_m_s_;Fort_2011.windSpeedAt80m_m_s_;Fort_2012.windSpeedAt80m_m_s_;
    Fort_2013.windSpeedAt80m_m_s_;Fort_2014.windSpeedAt80m_m_s_];
% Fort_all = Fort_all';
%% 统计

Fortmean = mean(Fort_all)
COmean = mean(CO_all)
Boisemean = mean(Boise_all)

Fortvar = var(Fort_all)
COvar = var(CO_all)
Boisevar = var(Boise_all)

rFC = corr(Fort_all, CO_all)
rFB = corr(Fort_all,Boise_all)
%% 模糊C聚类

maxIter = 100;   % 最大迭代次数
errorTol = 1e-5; % 收敛阈值
exponent = 2;    % 模糊系数
options = fcmOptions(NumClusters=n, ...
    MaxNumIteration=10000,Verbose = false);
% 执行模糊 C 均值聚类
[Fcenters, FU] = fcm(Fort_all, options);
[Ccenters, CU] = fcm(CO_all, options);
[Bcenters, BU] = fcm(Boise_all, options);

% 将每个数据点归类
[maxFU,IF] = max(FU);
[maxCU,IC] = max(CU);
[maxBU,IB] = max(BU);

%保存数据
save("windstates.mat", "Fcenters", "Ccenters","Bcenters","IF","IC","IB");
%% 计算状态概率转移矩阵
function p = state(I,n)
T = zeros(n,n);
D = zeros(n,1);
P = zeros(n,n);

for k = 1:length(I)
    state = I(k);
    D(state) = D(state)+1;
    if k ~= length(I)
    statenext = I(k+1);
    T(state,statenext) = T(state,statenext) + 1;
    end
end

for i = 1:n
    for j = 1:n
        if D(i) > 0
            P(i, j) = T(i, j) / D(i);
        else
            P(i, j) = 1 / 12; % 拉普拉斯平滑
        end
    end
end
p = P;
end

pF = state(IF,n);
pC = state(IC,n);
pB = state(IB,n);
save("p.mat", "pF", "pC","pB");
