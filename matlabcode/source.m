warning('off');
Hydro = readmatrix('E:\复现2\RTS-GMLC\RTS_Data\timeseries_data_files\Hydro\REAL_TIME_hydro.csv');
PV = readmatrix('E:\复现2\RTS-GMLC\RTS_Data\timeseries_data_files\PV\REAL_TIME_pv.csv');
CSP = readmatrix('E:\复现2\RTS-GMLC\RTS_Data\timeseries_data_files\CSP\REAL_TIME_Natural_Inflow.csv');
RTPV = readmatrix('E:\复现2\RTS-GMLC\RTS_Data\timeseries_data_files\RTPV\REAL_TIME_rtpv.csv');
Load = readmatrix('E:\复现2\RTS-GMLC\RTS_Data\timeseries_data_files\Load\REAL_TIME_regional_Load.csv');
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

save("power.mat", "HydroP", "CSPP","PVP","RTPVP","Load");
