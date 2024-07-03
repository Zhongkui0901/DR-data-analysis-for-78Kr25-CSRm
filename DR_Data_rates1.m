%% 选择至少两个TXT文件，获取电压
clear;
normalizaed = 0;%是否用冷却点归一化调制点
current_path = default_folder(true);
[filename,pathname] = uigetfile('*.txt','MultiSelect','on');
if isequal(pathname,0)
    cd(current_path)
    return;
end
file_numbers = length(filename);
detuning_voltages = cellfun(@get_voltage,filename').*[-1 1];

% start是NIM信号标志后的第start个点，例如5表示取NIM信号后第5个点开始参与计算
[Am,qcharge,Eion,Ie,radius_e,start,stop] = deal(56,15,4.35,0.04672,2.688,2,8);
%第1列DCCT，第2列Counts，默认第3列NIM3，第5列NIM1。根据实际NIM信号更改星号位置，%*f表示丢弃这一列的数据。
data_format = '%f %f %f %*f %*f %f %*f %*f %*f %*f';
%% 默认参数
[e,c,epsinon,lc,rad_tube,cooling_counts] = deal(1.602e-19,2.9979246e8,8.854187818e-12,4.0,12.5,30);
mu = vpa(931.4940954e6);%vpa函数保证计算更精确
me = vpa(0.5109989461e6);% 能量单位选择 eV
mion = Am * mu;
[gamaion,betaion]=Lorentz_factor(Eion * 1e6,mu);
%% 空间电荷效应及主高压，计算碰撞能量
Usp0 = -0.93 * (Ie / (4 * pi * epsinon * betaion * c)) * (1 + 2 * log(rad_tube / radius_e));
Ucath = -(gamaion - 1) * me + Usp0;
[gamaee,betaee]=Lorentz_factor(- (Ucath + detuning_voltages - Usp0),me);
Usp = Usp0 * betaion./ betaee;
[gamae,betae]=Lorentz_factor(- (Ucath + detuning_voltages - Usp),me); 
Erel = (sqrt(me^2 + mion^2 + 2*me*mion*gamaion.*gamae.*(1 - betaion.*betae)) - me - mion);
%% 参与计算的起始点
count_over_dcct = zeros(file_numbers,2);
Alpha_c = zeros(file_numbers,2);
sum_counts = zeros(file_numbers,2);
detuning_number = start:stop;
cooling_number = 5:cooling_counts;
absolute_factor = gamaion.^2 * qcharge * e^2 * c^2 * pi * radius_e^2 * 1e11 / Ie / lc;
%% 速率系数
for i = 1:file_numbers
    fid = fopen(fullfile(pathname,filename{i}));
    Data = textscan(fid,data_format);
    fclose(fid);
    DCCT = Data{1}; Count = Data{2};
% NIM3和NIM5分别是正负调制标志点
    NIM3 = Data{4}; NIM5 = Data{3};
    row_numbers = size(DCCT,1);
% 对DCCT线性插值，并使用外推
    nonzero_index = find(DCCT>0);
    F = griddedInterpolant(nonzero_index, DCCT(nonzero_index));
    disp(filename{i})%显示可能插值报错的文件名
    DCCT = F((1:row_numbers)');% 所有DCCT的行进行插值
% 寻找有效的行索引,flag变量是为了防止索引越界
    flag = row_numbers  - stop -1;
    %flag = row_numbers  - 6000;
    [NIM3_cooling, NIM3_index] = valid_index(NIM3,flag,cooling_number,detuning_number,cooling_counts);
    [NIM5_cooling, NIM5_index] = valid_index(NIM5,flag,cooling_number,detuning_number,cooling_counts);    
% 计算Counts除以DCCT的平均值
    count_over_dcct(i,:) = [mean(Count(NIM3_index)./DCCT(NIM3_index),'all'),...
    mean(Count(NIM5_index)./DCCT(NIM5_index),'all')];
    Alpha_c(i,:) = [mean(Count(NIM3_cooling)./DCCT(NIM3_cooling),'all'),...
    mean(Count(NIM5_cooling)./DCCT(NIM5_cooling),'all')]*betaion^2 * absolute_factor;
    sum_counts(i,:) = [sum(Count(NIM3_index),'all'),sum(Count(NIM5_index),'all')];
end
Alpha = count_over_dcct.*betae.*betaion * absolute_factor;
if normalizaed == true
    Alpha = Alpha./Alpha_c.*mean(Alpha_c);
end
ALL_DATA = [detuning_voltages Erel Alpha  Alpha_c sum_counts];
% 电压，碰撞能量，速率系数，统计计数等数据汇总并作图
columns = [1 3 5 7 9];
ALL_DATA = double(ALL_DATA(:,[columns columns+1]));
plot(ALL_DATA(:,[2 7 2 7]),ALL_DATA(:,[3 8 4 9]),'-*');
paras = sprintf('Am=%d,qcharge=%d,Ei=%.2f,start_stop_cooling_counts=%d-%d-%d,re=%.2f,Ie=%.3f,folder=%s',...
    Am,qcharge,Eion,start,stop,cooling_counts,radius_e,Ie*1e3,pathname);
cd(current_path);%返回程序的当前文件夹
%clearvars -except paras ALL_DATA
%save([paras,'.mat'])
%%
function detuning_voltage = get_voltage(cell_element)
    voltage_cell = strsplit(cell_element,'V');
    detuning_voltage = str2double(voltage_cell{1});
end
function [gamma,beta]=Lorentz_factor(kinetic,me)
    gamma = 1 + kinetic / me;
    beta = sqrt(1-1./ gamma.^2);
end
function [cooling, detuning] = valid_index(NIM,flag,cooling_number,detuning_number,cooling_counts)
NIM_index = find(NIM(1:flag)==1);
    if NIM_index(1) < cooling_counts + 1
        NIM_index(1) = [];
    end
    cooling = NIM_index - cooling_number;
    detuning = NIM_index + detuning_number;
end
function current_path = default_folder(is_default)
    fullpath = mfilename('fullpath');
    [current_path,~] = fileparts(fullpath);
    if is_default
        start_path = ['C:\Users\',getenv('username'),'\Desktop'];
        cd(start_path)%默认从桌面打开，也可以直接更改默认的起始路径
    end
end