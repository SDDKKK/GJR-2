warning('off');
% 获取当前脚本的完整路径
current_script_path = mfilename('fullpath');
% 获取当前脚本所在的目录
[current_dir, ~, ~] = fileparts(current_script_path);
% 获取项目根目录（当前目录的上一级）
project_root = fileparts(current_dir);
% 构建时间序列数据目录的路径
timeseries_dir = fullfile(project_root, 'RTS-GMLC', 'RTS_Data', 'timeseries_data_files');

% 使用相对路径读取数据文件
Hydro = readmatrix(fullfile(timeseries_dir, 'Hydro', 'REAL_TIME_hydro.csv'));
PV = readmatrix(fullfile(timeseries_dir, 'PV', 'REAL_TIME_pv.csv'));
CSP = readmatrix(fullfile(timeseries_dir, 'CSP', 'REAL_TIME_Natural_Inflow.csv'));
RTPV = readmatrix(fullfile(timeseries_dir, 'RTPV', 'REAL_TIME_rtpv.csv'));
Load = readmatrix(fullfile(timeseries_dir, 'Load', 'REAL_TIME_regional_Load.csv'));
Load = Load(1:end,:);
%% 处理为1h步长



function P = sumtoh(source,x,y)
n = 1;
for k = 1:8760
    for i = x:y
        P(k,i-x+1) = sum(source(n:n+11,i));
    end
    n = n+12;
end
end

HydroP = sumtoh(Hydro,5,24);
CSPP = sumtoh(CSP,5,5);
PVP = sumtoh(PV,5,29);
RTPVP = sumtoh(RTPV,5,35);
Load = sumtoh(Load,5,7);
%% 保存

% 使用相对路径保存mat文件
save(fullfile(current_dir, "power.mat"), "HydroP", "CSPP","PVP","RTPVP","Load");
