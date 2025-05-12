clear
clc
%%数据导入
warning('off')
% 获取当前脚本的完整路径
current_script_path = mfilename('fullpath');
% 获取当前脚本所在的目录
[current_dir, ~, ~] = fileparts(current_script_path);
% 获取项目根目录（当前目录的上一级）
project_root = fileparts(current_dir);
% 构建数据目录的路径
source_data_dir = fullfile(project_root, 'RTS-GMLC', 'RTS_Data', 'SourceData');

% 使用相对路径读取数据文件
gen = readtable(fullfile(source_data_dir, 'gen.csv'));
branch = readtable(fullfile(source_data_dir, 'branch.csv'));
dc_branch = readtable(fullfile(source_data_dir, 'dc_branch.csv'));
bus = readtable(fullfile(source_data_dir, 'bus.csv'));

% 使用相对路径加载本地mat文件
load(fullfile(current_dir, 'windstates.mat'));
load(fullfile(current_dir, 'p.mat'));
load(fullfile(current_dir, 'power.mat'));

%% 数据处理
%减去储能；加入一条直流支路
baseMVA = 100;
MTTR_b = branch.Duration;%修复时间
MTTR_g = gen.MTTRHr(1:end-1);
Line_fail = branch.PermOutRate;%故障率
Unit_fail = gen.FOR(1:end-1);
busLoad = (bus.MWLoad-bus.MVARLoad)/baseMVA;%节点负荷
nb = length(MTTR_b)+1;%线路数
nn = length(busLoad);%节点数
PF_MAX = [branch.ContRating/baseMVA;1];%线路最大容量
FrBranch = [branch.FromBus;dc_branch.FromBus];%首节点号
FrBranch = 24*((FrBranch-mod(FrBranch,100))./100-1) + mod(FrBranch,100);%修改编号
ToBranch = [branch.ToBus;dc_branch.ToBus];%末节点号
ToBranch = 24*((ToBranch-mod(ToBranch,100))./100-1) + mod(ToBranch,100);%修改编号
Branch_X = [branch.X;0.01];%电抗%%%%%%%%%%%%%%%%%%%%%%%%%%%直流线路如何解决？
Node_Branch = zeros(nn,nb);%节点-线路关联矩阵
for i = 1:nb
    Node_Branch(FrBranch(i),i) = 1;
    Node_Branch(ToBranch(i),i) = -1;
end
% 机组参数
unit_bus = gen.BusID(1:end-1);%机组对应节点
unit_bus = 24*((unit_bus-mod(unit_bus,100))./100-1) + mod(unit_bus,100);%修改编号
Pmax_unit = gen.PMaxMW(1:end-1);%机组最大出力
Pmin_unit = gen.PMinMW(1:end-1)/baseMVA;%机组最小出力
ng = length(unit_bus);%机组数量
Node_Unit = zeros(nn,ng);
for i = 1:ng
    Node_Unit(unit_bus(i),i) = 1;
end

Fail = [Line_fail;dc_branch.LineFORPerm;Unit_fail];
MTTR = [MTTR_b;dc_branch.MTTRLineHours;MTTR_g];

%% 优化模型建模
S = sdpvar(nn,1);%负荷削减
PF = sdpvar(nb,1);%线路潮流
P_unit = sdpvar(ng,1);%机组出力
Theta = sdpvar(nn,1);%相角
xg = sdpvar(ng,1);%机组故障标志
xl = sdpvar(nb,1);%线路故障标志
Lreal = sdpvar(nn,1);%负荷实际功率
Preal = sdpvar(ng,1); %机组实际最大出力
%L是节点负荷，A是节点-线路关联矩阵，node_number是节点数，branch_number是线路数，unit_number是机组数
%Trans_unit是节点-机组关联矩阵,X是线路电抗
F = [
    P_unit >= xg.*Pmin_unit;
    P_unit <= xg.*Preal;
    %相角约束
    Theta(1) == 0;
    Theta >= -pi;
    Theta <= pi;
    %线路约束
    PF >= -PF_MAX.*xl;
    PF <= PF_MAX.*xl;
    %潮流计算约束
    PF-pinv(diag(Branch_X))*Node_Branch'*Theta >= -10000*(1-xl);
    PF-pinv(diag(Branch_X))*Node_Branch'*Theta <= 10000*(1-xl);
    %节点功率约束
    Node_Unit*P_unit == Lreal-S+Node_Branch*PF;
    %削减量约束
    S <= Lreal;
    S >= 0;
    ];
ops = sdpsettings('solver','gurobi','verbose',0);
object = sum(S);
solve = optimizer(F,object,ops,[xl;xg;Lreal;Preal],{S;P_unit;Theta;PF});
%% SMC
state = ones(nb+ng,1);
VarianceCoefficient = 1;
s = 1;
t = 0;
statetime = zeros(nb+ng,1);
r = rand(nb+ng,1);
statetime(:,1) = ceil((-log(r)./Fail)*8760);  % 生成机组线路状态持续时间(单位小时)
windB = 1;%初始风速
wind = [Bcenters(windB);Bcenters(windB);Bcenters(windB);Bcenters(windB)];
Lreal = busLoad;
SystemEDNS = 0;SystemDNS = [];LOLP = 0;LOLF = 0;

% 初始化计时器和进度显示
sim_start_time = tic;
last_display_time = toc(sim_start_time);
disp('开始蒙特卡洛模拟...')
disp(['目标变异系数: ', num2str(0.01), ', 最小模拟年数: ', num2str(50)])

tic
while (VarianceCoefficient> 0.01||s<50)
    % 找到下一个发生状态转换的元件及其转换时间
    [Dk,TranstionElementIndex] = min(statetime(:,1));
    if t + Dk > 8760
        Dk = 8760-t;%不超过一年之内
    end
    xl = state(1:nb);
    xg = state(nb+1:nb+ng);
    DNS = 0;
    %某机组状态下每个小时的最小削负荷
    for i = 1:Dk
        windpower = 3*wind.^3;
        Preal = [Pmax_unit(1:74);HydroP(t+i,1:7)';Pmax_unit(82);HydroP(t+i,8:16)';Pmax_unit(92);
            HydroP(t+i,17:20)';PVP(t+i,1:20)';CSPP(t+i,1);PVP(t+i,21:25)';RTPVP(t+i,1:31)';
            min(windpower,800)]/baseMVA;
        tmpsol = solve{[xl;xg;Lreal;Preal]};
        DNS = DNS + baseMVA*sum(tmpsol{1})*1;%累加该状态下缺电量
        if DNS > 0
            LOLP = LOLP+1;%累加故障时间
        end

        rB = rand(12,1);
        windtime(:,1) = -log(rB)./pB(windB,:)';  % 切换风速状态
        [windmintime,windB] = min(windtime(:,1));
        wind = [Bcenters(windB);Bcenters(windB);Bcenters(windB);Bcenters(windB)];
    end

    if DNS > 0
       LOLF = LOLF+1;%累加故障次数
    end


    %计算可靠性指标
    SystemEDNS = SystemEDNS+DNS;
    SystemDNS = [SystemDNS DNS];
    if s >= 100%仿真时间大于100年
        EENS_mean = SystemEDNS/(s+t/8760);
        EENS_var = var(SystemDNS)/(s+t/8760);
        VarianceCoefficient = sqrt(EENS_var)/EENS_mean;
    end


    r = rand;% 切换机组线路状态
    if state(TranstionElementIndex,1) == 0
        statetime(TranstionElementIndex,1) = ceil((-log(r)./Fail(TranstionElementIndex,1))*8760);
    else
        statetime(TranstionElementIndex,1) = MTTR(TranstionElementIndex,1);
    end
    state(TranstionElementIndex,1) = 1-state(TranstionElementIndex,1);

    t = t + Dk;
    if t + Dk == 8760
        s = s+1;%年份加一
        t = 0;%重新仿真一年

        % 显示当前进度和统计信息

        % 计算当前可靠性指标（无论s是多少都计算）
        current_EENS = SystemEDNS/(s+t/8760);
        current_LOLP = LOLP/(s*8760+t);
        current_LOLF = LOLF/(s+t/8760);

        % 根据s是否大于等于100，决定显示格式
        if s >= 100
            % 显示包含变异系数的进度信息
            progress_msg = sprintf(['已模拟%d年, 当前变异系数: %.4f (目标: %.4f)\n' ...
                               '当前EENS: %.2f, LOLP: %.6f, LOLF: %.4f'], ...
                               s, VarianceCoefficient, 0.01, ...
                               current_EENS, current_LOLP, current_LOLF);
        else
            % 显示不包含变异系数的进度信息
            progress_msg = sprintf(['已模拟%d年, 最小需模拟年数: %d\n' ...
                               '当前EENS: %.2f, LOLP: %.6f, LOLF: %.4f'], ...
                               s, 50, current_EENS, current_LOLP, current_LOLF);
        end
        disp(progress_msg)
        disp('----------------------------------------')
    end

end

% 显示模拟完成信息和总耗时
total_sim_time = toc(sim_start_time);
if total_sim_time > 3600
    time_str = sprintf('%.2f小时', total_sim_time/3600);
elseif total_sim_time > 60
    time_str = sprintf('%.2f分钟', total_sim_time/60);
else
    time_str = sprintf('%.2f秒', total_sim_time);
end
disp('----------------------------------------')
disp('蒙特卡洛模拟完成')
disp(['总模拟时间: ', time_str])
disp(['总模拟年数: ', num2str(s)])
disp('----------------------------------------')

disp('MC year')
disp(s)
disp('EENS')
disp(EENS_mean)
disp('LOLP')
LOLP = LOLP/(s*8760+t);
disp(LOLP)
disp('LOLF')
LOLF = LOLF/(s+t/8760);
disp(LOLF)

toc