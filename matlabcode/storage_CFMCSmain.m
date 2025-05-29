clear
clc
%%数据导入
warning('off')
gen = readtable('E:\复现2\RTS-GMLC\RTS_Data\SourceData\gen.csv');
branch = readtable('E:\复现2\RTS-GMLC\RTS_Data\SourceData\branch.csv');
dc_branch = readtable('E:\复现2\RTS-GMLC\RTS_Data\SourceData\dc_branch.csv');
bus = readtable('E:\复现2\RTS-GMLC\RTS_Data\SourceData\bus.csv');

load('windstates.mat');
load('p.mat');
load('power.mat');

%% 数据处理
%减去储能；加入一条直流支路
baseMVA = 100;
MTTR_b = branch.Duration;%修复时间
MTTR_g = gen.MTTRHr(1:end-1);
Line_fail = branch.PermOutRate;%故障率
Unit_fail = gen.FOR(1:end-1);

busLoad = (bus.MWLoad-bus.MVARLoad);%节点负荷比例
sumL1 = sum(busLoad(1:24));
sumL2 = sum(busLoad(25:48));
sumL3 = sum(busLoad(49:73));
propL1 = busLoad(1:24)./sumL1;
propL2 = busLoad(25:48)./sumL2;
propL3 = busLoad(49:73)./sumL3;

nb = length(MTTR_b)+1;%线路数
nn = length(busLoad);%节点数
PF_MAX = [branch.ContRating/baseMVA;1];%线路最大容量
FrBranch = [branch.FromBus;dc_branch.FromBus];%首节点号
FrBranch = 24*((FrBranch-mod(FrBranch,100))./100-1) + mod(FrBranch,100);%修改编号
ToBranch = [branch.ToBus;dc_branch.ToBus];%末节点号
ToBranch = 24*((ToBranch-mod(ToBranch,100))./100-1) + mod(ToBranch,100);%修改编号
Branch_X = [branch.X];%电抗%%%%%%%%%%%%%%%%%%%%%%%%%%%直流线路已解决
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
%储能
ess_max = 148.7/baseMVA; 
essp_max = 102.5723/baseMVA;
%% 优化模型建模
S = sdpvar(nn,1);%负荷削减
PF = sdpvar(nb,1);%线路潮流
PF_ac = PF(1:end-1,1);%交流
PF_dc = PF(end,1);%直流
P_unit = sdpvar(ng,1);%机组出力
Theta = sdpvar(nn,1);%相角
xg = sdpvar(ng,1);%机组故障标志
xl = sdpvar(nb,1);%线路故障标志
Lreal = sdpvar(nn,1);%负荷实际功率
Preal = sdpvar(ng,1); %机组实际最大出力
ess = sdpvar(nn,1);%储能容量
essp = sdpvar(nn,1);%储能出力,负数表示放电，正数表示充电

%L是节点负荷，A是节点-线路关联矩阵，node_number是节点数，branch_number是线路数，unit_number是机组数
%Trans_unit是节点-机组关联矩阵,X是线路电抗
F = [
    essp(1:60,:)==0;essp(62:end,:)==0;
    essp(61,:)>=-essp_max;essp(61,:)<=essp_max;
    ess(61,:)+essp(61,:)<=ess_max;
    ess(61,:)+essp(61,:)>=0;
    %储能相关约束
    P_unit >= xg.*Pmin_unit;
    P_unit <= xg.*Preal;
    %相角约束
    Theta(13) == 0;
    Theta >= -pi;
    Theta <= pi;
    %线路约束
    PF >= -PF_MAX.*xl;
    PF <= PF_MAX.*xl;
    %潮流计算约束
    PF_ac-pinv(diag(Branch_X))*Node_Branch(:,1:end-1)'*Theta >= -10000*(1-xl(1:end-1,1));
    PF_ac-pinv(diag(Branch_X))*Node_Branch(:,1:end-1)'*Theta <= 10000*(1-xl(1:end-1,1));
    %节点功率约束
    Node_Unit*P_unit - essp == Lreal-S+Node_Branch*PF;
    %削减量约束
    S <= Lreal;
    S >= 0;
    ];
ops = sdpsettings('solver','gurobi','verbose',0);
object = 1000*sum(S)-essp(61,1);
solve = optimizer(F,object,ops,[xl;xg;Lreal;Preal;ess],{S;P_unit;Theta;PF;essp});
%% SMC
state = ones(nb+ng,1);
VarianceCoefficient = 1;
s = 1;
t = 0;
statetime = zeros(nb+ng,1);
r = rand(nb+ng,1);
statetime(:,1) = ceil((-log(r)./Fail)*8760);  % 生成机组线路状态持续时间(单位小时）

ess(61,1) = ess_max;%初始储能
ess(1:60,:) = 0;
ess(62:end,:) = 0;

windC = 1;windF = 1;%初始风速
wind = [Ccenters(windC);Fcenters(windF);Fcenters(windF);Fcenters(windF)];
SystemEDNS = 0;SystemDNS = [];LOLP = 0;LOLF = 0;
DNS2 = zeros(2,1);%用于判断故障发生的时刻
tic
while (VarianceCoefficient> 0.01 && s<300)
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
        windpower = 1*(wind).^3;
        Preal = [Pmax_unit(1:74);HydroP(t+i,1:7)';Pmax_unit(82);HydroP(t+i,8:16)';Pmax_unit(92);
            HydroP(t+i,17:20)';PVP(t+i,1:20)';CSPP(t+i,1);PVP(t+i,21:25)';RTPVP(t+i,1:31)';
            min(windpower,800)]./baseMVA;
        %Lreal = [propL1*Load(t+i,1);propL2*Load(t+i,2);propL3*Load(t+i,3)]./baseMVA;
        Lreal = busLoad./baseMVA;
        Lreal = Lreal*1.27;
        tmpsol = solve{[xl;xg;Lreal;Preal;ess]};
        ess = ess+tmpsol{5};%更新储能状态

        DNS2(1,1) = DNS2(2,1);DNS2(2,1) = sum(tmpsol{1});
        if DNS2(2,1) - DNS2(1,1) == DNS2(2,1) && DNS2(2,1)~=0
            LOLF = LOLF+1;%累加故障次数
        end
        
        rC = rand(12,1);
        windtimeC(:,1) = -log(rC)./pC(windC,:)';  % 切换风速状态
        [~,windC] = min(windtimeC(:,1));
        rF = rand(12,1);
        windtimeF(:,1) = -log(rF)./pF(windF,:)';  % 切换风速状态
        [~,windF] = min(windtimeF(:,1));

        wind = [Ccenters(windC);Fcenters(windF);Fcenters(windF);Fcenters(windF)];
    end

    DNS = DNS + baseMVA*sum(tmpsol{1})*1;%累加该状态下缺电量
    if sum(tmpsol{1}) > 0
        LOLP = LOLP+1;%累加故障时间
    end

    %计算可靠性指标
    SystemEDNS = SystemEDNS+DNS;
    SystemDNS = [SystemDNS DNS];
    EENS_mean = SystemEDNS/(s+t/8760);


    if s >= 100%仿真时间大于100年
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
        disp('MC year')
        disp(s)
        disp('EENS')
        disp(EENS_mean)

        t = 0;%重新仿真一年
    end

end
Press = SystemEDNS/LOLF;

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
disp('Press')
disp(Press)
toc