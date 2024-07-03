%% ѡ����������TXT�ļ�����ȡ��ѹ
clear;
normalizaed = 0;%�Ƿ�����ȴ���һ�����Ƶ�
current_path = default_folder(true);
[filename,pathname] = uigetfile('*.txt','MultiSelect','on');
if isequal(pathname,0)
    cd(current_path)
    return;
end
file_numbers = length(filename);
detuning_voltages = cellfun(@get_voltage,filename').*[-1 1];

% start��NIM�źű�־��ĵ�start���㣬����5��ʾȡNIM�źź��5���㿪ʼ�������
[Am,qcharge,Eion,Ie,radius_e,start,stop] = deal(56,15,4.35,0.04672,2.688,2,8);
%��1��DCCT����2��Counts��Ĭ�ϵ�3��NIM3����5��NIM1������ʵ��NIM�źŸ����Ǻ�λ�ã�%*f��ʾ������һ�е����ݡ�
data_format = '%f %f %f %*f %*f %f %*f %*f %*f %*f';
%% Ĭ�ϲ���
[e,c,epsinon,lc,rad_tube,cooling_counts] = deal(1.602e-19,2.9979246e8,8.854187818e-12,4.0,12.5,30);
mu = vpa(931.4940954e6);%vpa������֤�������ȷ
me = vpa(0.5109989461e6);% ������λѡ�� eV
mion = Am * mu;
[gamaion,betaion]=Lorentz_factor(Eion * 1e6,mu);
%% �ռ���ЧӦ������ѹ��������ײ����
Usp0 = -0.93 * (Ie / (4 * pi * epsinon * betaion * c)) * (1 + 2 * log(rad_tube / radius_e));
Ucath = -(gamaion - 1) * me + Usp0;
[gamaee,betaee]=Lorentz_factor(- (Ucath + detuning_voltages - Usp0),me);
Usp = Usp0 * betaion./ betaee;
[gamae,betae]=Lorentz_factor(- (Ucath + detuning_voltages - Usp),me); 
Erel = (sqrt(me^2 + mion^2 + 2*me*mion*gamaion.*gamae.*(1 - betaion.*betae)) - me - mion);
%% ����������ʼ��
count_over_dcct = zeros(file_numbers,2);
Alpha_c = zeros(file_numbers,2);
sum_counts = zeros(file_numbers,2);
detuning_number = start:stop;
cooling_number = 5:cooling_counts;
absolute_factor = gamaion.^2 * qcharge * e^2 * c^2 * pi * radius_e^2 * 1e11 / Ie / lc;
%% ����ϵ��
for i = 1:file_numbers
    fid = fopen(fullfile(pathname,filename{i}));
    Data = textscan(fid,data_format);
    fclose(fid);
    DCCT = Data{1}; Count = Data{2};
% NIM3��NIM5�ֱ����������Ʊ�־��
    NIM3 = Data{4}; NIM5 = Data{3};
    row_numbers = size(DCCT,1);
% ��DCCT���Բ�ֵ����ʹ������
    nonzero_index = find(DCCT>0);
    F = griddedInterpolant(nonzero_index, DCCT(nonzero_index));
    disp(filename{i})%��ʾ���ܲ�ֵ������ļ���
    DCCT = F((1:row_numbers)');% ����DCCT���н��в�ֵ
% Ѱ����Ч��������,flag������Ϊ�˷�ֹ����Խ��
    flag = row_numbers  - stop -1;
    %flag = row_numbers  - 6000;
    [NIM3_cooling, NIM3_index] = valid_index(NIM3,flag,cooling_number,detuning_number,cooling_counts);
    [NIM5_cooling, NIM5_index] = valid_index(NIM5,flag,cooling_number,detuning_number,cooling_counts);    
% ����Counts����DCCT��ƽ��ֵ
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
% ��ѹ����ײ����������ϵ����ͳ�Ƽ��������ݻ��ܲ���ͼ
columns = [1 3 5 7 9];
ALL_DATA = double(ALL_DATA(:,[columns columns+1]));
plot(ALL_DATA(:,[2 7 2 7]),ALL_DATA(:,[3 8 4 9]),'-*');
paras = sprintf('Am=%d,qcharge=%d,Ei=%.2f,start_stop_cooling_counts=%d-%d-%d,re=%.2f,Ie=%.3f,folder=%s',...
    Am,qcharge,Eion,start,stop,cooling_counts,radius_e,Ie*1e3,pathname);
cd(current_path);%���س���ĵ�ǰ�ļ���
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
        cd(start_path)%Ĭ�ϴ�����򿪣�Ҳ����ֱ�Ӹ���Ĭ�ϵ���ʼ·��
    end
end