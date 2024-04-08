function [Timer, Reboot] = Main_Catia(Dir_Total, Mode, Num_Infill) %Manual_1,DOE_2,Infill_3

%**********************************
%|IGP取值范围|
%     c1  [0.15, 0.50]     0-100% 
%     c2  [0.10, 0.70]     0-100%
%     c3  [-0.01, 0.10]    1-90%
%     c4  [-0.04, 0.08]    1-90%
% 
%    X_T  [0.20, 0.45]     1-99%
%      T  [0.04, 0.21]     1-99%
%   Rho0  [0.18, 0.71]     1-95%
%Beta_TE  [1.00, 2.53]     25-95%
%**********************************

if Mode == 1
    FileName_AddOn = '_Manual';
elseif Mode == 2
    FileName_AddOn = '_DOE';
elseif Mode == 3
    FileName_AddOn = '_Infill';
end

if Mode == 1                                                                %计算节点数量
    Num_Thread = 1;
elseif Mode == 2
    Num_Thread = 4;
elseif Mode == 3
    Num_Thread = Num_Infill;
end
% p = gcp('nocreate');
% if isempty(p)
%     parpool(8);
% elseif p.NumWorkers ~= 8
%     delete(p);
%     parpool(8);
% end


%% 文件目录
Dir.Total = Dir_Total; 
Dir.Catia = [Dir.Total,'\','Catia'];


%% 参数_程序
Para_Catia.Common.FileName_Load_DesignVar_UnPre_NoS = ['DesignVar_UnPre',FileName_AddOn];
Para_Catia.Common.FileName_LoSa_DesignVar = ['CatiaR_DesignVar',FileName_AddOn];

Para_Catia.Common.FileName_LoSa_Geometry = ['CatiaR_Geometry',FileName_AddOn];

Para_Catia.Common.FileName_Save_CATPart = 'zFlyingWing';

Para_Catia.GenWingPlan.Act_Draw = 0;                                        %是否生成机翼平面图

Para_Catia.GenSection.Density_Airfoil = 100;                                %翼型弦向坐标数量
Para_Catia.GenSection.Diameter_Fillet = 1;                                  %机翼后缘圆角直径(不后处理时自动为0)

Para_Catia.GenCATPart.Act_Debug = 0;                                        %是否开启Debug模式(遇错不停止并保存成.CATPart.Error)
Para_Catia.GenCATPart.Act_Demo = 0;                                         %是否开启演示模式
Para_Catia.GenCATPart.Index_SurfConnect = [2,2,1];                          %曲面建模展向连接索引
Para_Catia.GenCATPart.Mode_GuideLine = 1;                                   %前后缘引导线模式,线性_1,曲线_2(连接索引同Index_SurfConnect)

Para_Catia.GenCATPart.Mode_Post = 1;                                        %是否后处理,不处理_0,光顺+ICEM格式_1
if Para_Catia.GenCATPart.Mode_Post == 0                                     %***不处理***
    Para_Catia.GenSection.Diameter_Fillet = 0;
    Para_Catia.GenCATPart.Diameter_Fillet = 0;
end
if Para_Catia.GenCATPart.Mode_Post == 1                                     %***光顺+ICEM格式***
    Para_Catia.GenCATPart.Diameter_Fillet...
    = Para_Catia.GenSection.Diameter_Fillet;
    
    Para_Catia.GenCATPart.Rate_Fluid = 30;                                  %延伸倍率(相对根弦)_外流场
    Para_Catia.GenCATPart.Rate_BOI_Fluid = 14;                              %延伸倍率(相对MAC)_加密场

    Para_Catia.GenCATPart.ScaleRate = 0.4;                                  %模型缩放比例
    Para_Catia.GenCATPart.Index_NonActSection_ICEM = 4 ;                    %不参与生成ICEM拓扑的截面
    Para_Catia.GenCATPart.Index_Layer2Connect_ICEM = [2,2,1];               %ICEM的Layer2拓扑的展向连接索引
end
if Para_Catia.GenCATPart.Mode_Post == 2
end

Para_Catia.GenCATPart.Act_SaveStl = 1;                                      %是否保存为.stl格式
Para_Catia.GenCATPart.Act_SaveIgs = 1;                                      %是否保存为.igs格式


%% 参数_机翼_Manual
%********************************************************************************************
%|单位|
%长度：mm
%角度：Deg
%
%|规范|
%根弦：全部须赋值
%各个翼段：
% 1)展向位置：比例参数定义为(优先比例[]+1e5,其余[]+1e10),最外侧须赋值
% 2)梢弦长：自动填充为NaN,相同参数定义为1e30(与外段翼一致,可连续设置),最外侧须赋值
% 3)梢部扭转角：自动填充为NaN,相似参数定义为[]+1e2,相同参数定义为1e30(与内段翼一致)
% 4)前缘后掠角：自动填充为NaN(与外段翼一致),最外侧须赋值
% 5)上反角：自动填充为NaN(与外段翼一致),最外侧须赋值
% 6)翼型：相似参数定义为[]+1e2(不可连续设置),相同参数定义为1e30(与内段翼一致,弯度厚度可分开设置)
%********************************************************************************************

if Mode == 1
    Mode_DataSource = 1;                                                    %参数化建模_1,引用外部数据_2
    
    if Mode_DataSource == 1
        %根弦
        %                  前缘坐标(x, y, z)  根弦长  扭转角  翼型(c1, c2, c3, c4,   X_T, T, Rho0, Beta_TE)	
        Data_Import.Root = [0, 0, 0,          6500   0       0.3348, 0.6184, 0.0334, 0.0163,...
                                                             0.2964, 0.1600, 0.3778, 1.8354];

        %每个翼段
        %                   展向位置    梢弦长    梢部扭转角  前缘后掠角  上反角  翼型(c1, c2, c3, c4,   X_T, T, Rho0, Beta_TE)
        Data_Import.Panel = [1/2+1e10   NaN       0          NaN        NaN     0.3348, 0.6184, 0.0334, 0.0163,...
                                                                            0.2964, 0.1400, 0.3778, 1.8354;
                                                                           
                             3250       2000      0          NaN        -6      0.3348, 0.6184, 0.0334, 0.0163,...
                                                                                0.2964, 0.1500, 0.3778, 1.8354;
                                                                               
                             1/2+1e10   NaN       NaN        NaN        NaN     0.3348, 0.6184, 0.0334, 0.0163,...
                                                                                0.2964, 0.1400, 0.3778, 1.8354;
                                                                               
                             7555.6     2000      -2         NaN        -2      0.3348, 0.6184, 0.0334, 0.0163,...
                                                                                0.2964, 0.1200, 0.3778, 1.8354;
                                                                               
                             8927.8     100       -2          35         9      0.3348, 0.6184, 0.0000, 0.0000,...
                                                                            0.2964, 0.1200, 0.3778, 1.8354];
        
    elseif Mode_DataSource == 2
        load('ToolR_Catia_Inverse.mat');
        
        Data_Import = ToolR_Catia_Inverse.Data_Import;
        Data_Section = ToolR_Catia_Inverse.Data_Section;
        Data_Airfoil = ToolR_Catia_Inverse.Data_Airfoil;
    end
    
    Data_Import.FlightState = ones(1,3) * NaN;
end


%% 初始化
if ~(Mode == 1 && Mode_DataSource == 2)
    Mode_DataSource = 1;
    Data_Section = [];
    Data_Airfoil = [];
end

if Mode == 2
    load(['DOER_',Para_Catia.Common.FileName_Load_DesignVar_UnPre_NoS,'.mat'])
    eval(['Data_Import=', 'DOER_',Para_Catia.Common.FileName_Load_DesignVar_UnPre_NoS,';']);
end

if Mode == 3
    %加载DOE数据
    Para_Catia.Common.FileName_Load_DesignVar_UnPre_NoS = replace(Para_Catia.Common.FileName_Load_DesignVar_UnPre_NoS,'_Infill','_DOE');
    load(['DOER_',Para_Catia.Common.FileName_Load_DesignVar_UnPre_NoS,'.mat']);
    eval(['Data_Import_DOE=', 'DOER_',Para_Catia.Common.FileName_Load_DesignVar_UnPre_NoS,';']);
    
    Para_Catia.Common.FileName_LoSa_Geometry = replace(Para_Catia.Common.FileName_LoSa_Geometry,'_Infill','_DOE');
    load([Para_Catia.Common.FileName_LoSa_Geometry,'.mat']);
    eval(['Data_Geometry_DOE=', Para_Catia.Common.FileName_LoSa_Geometry,';']);
    
    Para_Catia.Common.FileName_LoSa_DesignVar = replace(Para_Catia.Common.FileName_LoSa_DesignVar,'_Infill','_DOE');
    load([Para_Catia.Common.FileName_LoSa_DesignVar,'.mat']);
    eval(['Data_DesignVar_DOE=', Para_Catia.Common.FileName_LoSa_DesignVar,';']);
    
    %加载Infill数据
    Para_Catia.Common.FileName_Load_DesignVar_UnPre_NoS = replace(Para_Catia.Common.FileName_Load_DesignVar_UnPre_NoS,'_DOE','_Infill');
    load(['OptimizeR_',Para_Catia.Common.FileName_Load_DesignVar_UnPre_NoS,'.mat']);
    eval(['Data_Import_Infill=', 'OptimizeR_',Para_Catia.Common.FileName_Load_DesignVar_UnPre_NoS,';']);
    
    Para_Catia.Common.FileName_LoSa_Geometry = replace(Para_Catia.Common.FileName_LoSa_Geometry,'_DOE','_Infill');
    if exist([Para_Catia.Common.FileName_LoSa_Geometry,'.mat'],'file')
        load([Para_Catia.Common.FileName_LoSa_Geometry,'.mat']);
        eval(['Data_Geometry_Infill=', Para_Catia.Common.FileName_LoSa_Geometry,';']);
    else
        Data_Geometry_Infill = '';
    end
    
    Para_Catia.Common.FileName_LoSa_DesignVar = replace(Para_Catia.Common.FileName_LoSa_DesignVar,'_DOE','_Infill');
    if exist([Para_Catia.Common.FileName_LoSa_DesignVar,'.mat'],'file')
        load([Para_Catia.Common.FileName_LoSa_DesignVar,'.mat']);
        eval(['Data_DesignVar_Infill=', Para_Catia.Common.FileName_LoSa_DesignVar,';']);
    else
        Data_DesignVar_Infill = '';
    end
    
    %合并
    Data_Import = [Data_Import_DOE, Data_Import_Infill];
    Data_Geometry = [Data_Geometry_DOE, Data_Geometry_Infill]; 
    Data_DesignVar = [Data_DesignVar_DOE, Data_DesignVar_Infill];
end


%% 数据校验
if ~strcmp(pwd, Dir.Total) ||...
   ~exist(Dir_Total,'dir') || strcmp(Dir_Total(end),'\') || ~ischar(Dir_Total)
    error('*** Error：Dir_Total ***');
end

% if Mode == 1
%     FileName_Mat = dir('*.mat');
%     for nMat = 1:length(FileName_Mat)
%         if contains(FileName_Mat(nMat).name, '_Manual')
%             error('*** Error：_Manual is already exist ***');
%         end
%     end
% end

if Mode ~= 3 && nargin > 2
    error('*** Error：Mode & Input ***');
end

if Mode == 3 && nargin < 3
    error('*** Error：Mode & Input ***');
end

if Mode ~= 1
    if mod(length(Data_Import),Num_Thread) ~= 0
        error('*** Error：Num_Thread & Num_Catia ***');
    end
end

if contains(Para_Catia.Common.FileName_Save_CATPart, ["_A","_B"])
    error('*** Error：Common.FileName_Save_CATPart ***');
end

if Mode == 2 || Mode == 3
    if ~contains(Dir_Total, Para_Catia.Common.FileName_Save_CATPart(2:end))
        error('*** Error：Common.FileName_Save_CATPart ***');
    end
end

if Para_Catia.GenCATPart.Mode_Post == 1 || Para_Catia.GenCATPart.Mode_Post == 2
    if Mode == 1 && Mode_DataSource == 2
        if (Data_Section(1).Up(end,3) - Data_Section(1).Low(end,3)) > Para_Catia.GenSection.Diameter_Fillet
            error('*** Error：GenSection.Diameter_Fillet ***');
        end
    end
    
    % if Mode == 1 || Mode == 2
    %     Span_Half = Data_Import(1).Panel(end,1);
    % elseif Mode == 3
    %     Span_Half = Data_Geometry(1).b/2 * 1000;
    % end
    % if Span_Half > 4000 && Para_Catia.GenSection.Diameter_Fillet <= 2 ||...
    %    Span_Half < 4000 && Para_Catia.GenSection.Diameter_Fillet > 2
    %     error('*** Error：GenSection.Diameter_Fillet ***');
    % end
end

if Para_Catia.GenSection.Diameter_Fillet ~= 0 && Para_Catia.GenCATPart.Mode_Post == 0  ||...
   Para_Catia.GenSection.Diameter_Fillet == 0 && Para_Catia.GenCATPart.Mode_Post ~= 0
    error('*** Error：GenSection.Diameter_Fillet & GenCATPart.Mode_Post ***');
end

if sum(Para_Catia.GenCATPart.Index_SurfConnect) ~= size(Data_Import(1).Panel,1)
    error('*** Error：GenCATPart.Index_SurfConnect ***')
end

if Para_Catia.GenCATPart.Mode_Post == 1
    if Mode == 1 || Mode == 2
        Span_Half = Data_Import(1).Panel(end,1);
    elseif Mode == 3
        Span_Half = Data_Geometry(1).b/2 * 1000;
    end
    if (Span_Half < 4000 && Para_Catia.GenCATPart.ScaleRate ~= 1) ||...
        Span_Half * Para_Catia.GenCATPart.ScaleRate > 4000
        error('*** Error：GenCATPart.ScaleRate ***');
    end
end

if Para_Catia.GenCATPart.Mode_Post == 1
    if min(Para_Catia.GenCATPart.Index_NonActSection_ICEM) < 0 ||...
       ismember(1,Para_Catia.GenCATPart.Index_NonActSection_ICEM) ||...
       max(Para_Catia.GenCATPart.Index_NonActSection_ICEM) > size(Data_Import(1).Panel,1)
        error('*** Error：GenCATPart.Index_NonActSection_ICEM ***');
    end
    if sum(Para_Catia.GenCATPart.Index_NonActSection_ICEM) ~= 0 &&...
       (size(Data_Import(1).Panel,1) - length(Para_Catia.GenCATPart.Index_NonActSection_ICEM) + 1 ~= 5)
        error('*** Error：Data_Import.Panel & GenCATPart.Index_NonActSection_ICEM ***');
    end
end

if Para_Catia.GenCATPart.Mode_Post == 1
    if length(Para_Catia.GenCATPart.Index_Layer2Connect_ICEM) > 3
        error('*** Error：GenCATPart.Index_Layer2Connect_ICEM ***');
    end
    if sum(Para_Catia.GenCATPart.Index_Layer2Connect_ICEM) ~= size(Data_Import(1).Panel,1)
        error('*** Error：GenCATPart.Index_Layer2Connect_ICEM ***');
    end
end
    

%% 主程序
if Mode == 1 || Mode == 2
    nStart_Catia = 1;
    nStart_Data = 1;
end
if Mode == 3
    nStart_Catia = length(Data_Import) - Num_Infill + 1;
    nStart_Data = length(Data_DesignVar_DOE) + 1;
end

T0 = clock;
fprintf('******************** Catia ********************\n');
fprintf('Process Data...\n');

% parfor nCatia = nStart_Catia:length(Data_Import)
for nCatia = nStart_Catia:length(Data_Import)
    %输入数据处理
    Data_DesignVar(nCatia) = Catia_PreImport(Data_Import(nCatia));

    %计算机翼平面数据
    Data_Geometry(nCatia) = Catia_CalGeometry(Data_DesignVar(nCatia));
   
    %生成各截面坐标
    if Mode == 1 && Mode_DataSource == 2
        Data_Geometry(nCatia).Section = Data_Section;
        Data_Geometry(nCatia).Airfoil = Data_Airfoil;
    else
        [Data_Geometry(nCatia).Section, Data_Geometry(nCatia).Airfoil, Data_DesignVar(nCatia)]...
        = Catia_GenSection(Para_Catia, Data_DesignVar(nCatia));
    end
    
    Data_Geometry(nCatia).Diameter_Fillet = Para_Catia.GenSection.Diameter_Fillet;
    if Para_Catia.GenCATPart.Mode_Post == 1 ||...
       Para_Catia.GenCATPart.Mode_Post == 2
        Data_Geometry(nCatia).Rate.Fluid = Para_Catia.GenCATPart.Rate_Fluid;
        Data_Geometry(nCatia).Rate.BOI_Fluid = Para_Catia.GenCATPart.Rate_BOI_Fluid;
    end
    if Para_Catia.GenCATPart.Mode_Post == 1
        Data_Geometry(nCatia).ScaleRate = Para_Catia.GenCATPart.ScaleRate;
    end
end
eval([Para_Catia.Common.FileName_LoSa_DesignVar, '=Data_DesignVar(nStart_Data:end);']);
save(Para_Catia.Common.FileName_LoSa_DesignVar, Para_Catia.Common.FileName_LoSa_DesignVar)
eval([Para_Catia.Common.FileName_LoSa_Geometry, '=Data_Geometry(nStart_Data:end);']);
save(Para_Catia.Common.FileName_LoSa_Geometry, Para_Catia.Common.FileName_LoSa_Geometry);

%绘制平面图
if Para_Catia.GenWingPlan.Act_Draw == 1
    close all;
    for nCatia = nStart_Catia:length(Data_Import)
        Catia_DrawWingPlan(Data_Import(nCatia), Data_Geometry(nCatia), nCatia);
    end
end

%Catia自动生成
fprintf('Generate CATPart...\n');
for nIter = ((nStart_Catia-1)/Num_Thread+1):length(Data_Import)/Num_Thread
    
    [~,~] = dos('taskkill /s SERVER0 /f /t /im EXCEL.EXE');
    
    fprintf('剩余：%d\n',length(Data_Import)-(nIter-1)*Num_Thread);
    
    try
        % parfor nThread = 1:Num_Thread
        for nThread = 1:Num_Thread
            nCatia = (nIter-1)*Num_Thread + nThread;
            Reboot(nThread) = Catia_GenCATPart(Dir, Para_Catia, Data_Geometry(nCatia), nCatia, nThread, Mode);
        end
    catch                                                                   %Excel偶发Bug
        [~,~] = dos('taskkill /s SERVER0 /f /t /im EXCEL.EXE');
        fprintf('\nRetry...\n');
        
        % parfor nThread = 1:Num_Thread
        for nThread = 1:Num_Thread
            nCatia = (nIter-1)*Num_Thread + nThread;
            Reboot(nThread) = Catia_GenCATPart(Dir, Para_Catia, Data_Geometry(nCatia), nCatia, nThread, Mode);
        end
    end
    
    [~,~] = dos('taskkill /s SERVER0 /f /t /im EXCEL.EXE');
    
    %计算体积
    if Para_Catia.GenCATPart.Act_SaveStl == 1
        for nThread = 1:Num_Thread
            nCatia = (nIter-1)*Num_Thread + nThread;
            if Mode == 1
                FileName_CATPart = Para_Catia.Common.FileName_Save_CATPart;
            else
                FileName_CATPart = [Para_Catia.Common.FileName_Save_CATPart,'_',num2str(nCatia,'%04d')];
            end
            Dir_Stl = [Dir.Catia,'\',FileName_CATPart,'.stl'];
            
            if exist(Dir_Stl,'file')
                Data_Stl = stlread([Dir.Catia,'\',FileName_CATPart,'.stl']);
                Data_Geometry(nCatia).Volume = StlVolume(Data_Stl.Points, Data_Stl.ConnectivityList) / 1e9;
            end
        end
        
        eval([Para_Catia.Common.FileName_LoSa_Geometry, '=Data_Geometry(nStart_Data:end);']);
        save(Para_Catia.Common.FileName_LoSa_Geometry, Para_Catia.Common.FileName_LoSa_Geometry);
    end
    
    %后处理
    for nThread = 1:Num_Thread
        nCatia = (nIter-1)*Num_Thread + nThread;
        % parfeval(@Catia_Post, 0, Dir, Para_Catia, Data_DesignVar(nCatia), Data_Geometry(nCatia), nCatia, Mode, 1);
        Catia_Post(Dir, Para_Catia, Data_DesignVar(nCatia), Data_Geometry(nCatia), nCatia, Mode, 1);
    end
end
 
fprintf('***********************************************\n');
Timer = round(etime(clock,T0));
fprintf('耗时：%d秒\n\n',Timer);


end


















