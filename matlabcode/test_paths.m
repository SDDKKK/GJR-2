% 测试脚本，验证所有相对路径是否正确
disp('开始测试相对路径...');

% 获取当前脚本的完整路径
current_script_path = mfilename('fullpath');
% 获取当前脚本所在的目录
[current_dir, ~, ~] = fileparts(current_script_path);
% 获取项目根目录（当前目录的上一级）
project_root = fileparts(current_dir);

% 构建数据目录的路径
source_data_dir = fullfile(project_root, 'RTS-GMLC', 'RTS_Data', 'SourceData');
timeseries_dir = fullfile(project_root, 'RTS-GMLC', 'RTS_Data', 'timeseries_data_files');
wind_data_dir = fullfile(current_dir, 'data2');

% 测试SourceData目录下的文件
try
    disp('测试SourceData目录下的文件...');
    gen = readtable(fullfile(source_data_dir, 'gen.csv'));
    branch = readtable(fullfile(source_data_dir, 'branch.csv'));
    dc_branch = readtable(fullfile(source_data_dir, 'dc_branch.csv'));
    bus = readtable(fullfile(source_data_dir, 'bus.csv'));
    disp('SourceData文件加载成功！');
catch e
    disp('SourceData文件加载失败！');
    disp(e.message);
end

% 测试timeseries_data_files目录下的文件
try
    disp('测试timeseries_data_files目录下的文件...');
    Hydro = readmatrix(fullfile(timeseries_dir, 'Hydro', 'REAL_TIME_hydro.csv'));
    PV = readmatrix(fullfile(timeseries_dir, 'PV', 'REAL_TIME_pv.csv'));
    CSP = readmatrix(fullfile(timeseries_dir, 'CSP', 'REAL_TIME_Natural_Inflow.csv'));
    RTPV = readmatrix(fullfile(timeseries_dir, 'RTPV', 'REAL_TIME_rtpv.csv'));
    Load = readmatrix(fullfile(timeseries_dir, 'Load', 'REAL_TIME_regional_Load.csv'));
    disp('时间序列数据文件加载成功！');
catch e
    disp('时间序列数据文件加载失败！');
    disp(e.message);
end

% 测试wind_data_dir目录下的文件
try
    disp('测试wind_data_dir目录下的文件...');
    if exist(fullfile(wind_data_dir, 'Boise', 'Boise_2007.csv'), 'file')
        disp('风速数据文件存在！');
        % 尝试加载一个文件作为示例
        Boise_2007 = readtable(fullfile(wind_data_dir, 'Boise', 'Boise_2007.csv'));
        disp('风速数据文件加载成功！');
    else
        disp('风速数据文件不存在！');
    end
catch e
    disp('风速数据文件加载失败！');
    disp(e.message);
end

% 测试本地mat文件
try
    disp('测试本地mat文件...');
    if exist(fullfile(current_dir, 'windstates.mat'), 'file') && ...
       exist(fullfile(current_dir, 'p.mat'), 'file') && ...
       exist(fullfile(current_dir, 'power.mat'), 'file')
        disp('本地mat文件存在！');
        % 尝试加载mat文件
        load(fullfile(current_dir, 'windstates.mat'));
        load(fullfile(current_dir, 'p.mat'));
        load(fullfile(current_dir, 'power.mat'));
        disp('本地mat文件加载成功！');
    else
        disp('部分本地mat文件不存在！');
    end
catch e
    disp('本地mat文件加载失败！');
    disp(e.message);
end

disp('路径测试完成！');
