function Reboot = Catia_GenCATPart(Dir, Para_Catia, Data_Geometry, nCatia, nThread, Mode)


Dir_Catia = Dir.Catia;

FileName_Save_CATPart = Para_Catia.Common.FileName_Save_CATPart;
if Mode == 1
    FileName_CATPart = FileName_Save_CATPart;
else
    FileName_CATPart = [FileName_Save_CATPart,'_',num2str(nCatia,'%04d')];
end

MAC_Length = Data_Geometry.MAC_Length;
MAC_x = Data_Geometry.MAC_x;
Data_Section = Data_Geometry.Section;
Data_Airfoil = Data_Geometry.Airfoil;

Act_Debug = Para_Catia.GenCATPart.Act_Debug;
Act_Demo = Para_Catia.GenCATPart.Act_Demo;
Index_SurfConnect = Para_Catia.GenCATPart.Index_SurfConnect;
Mode_GuideLine = Para_Catia.GenCATPart.Mode_GuideLine;

Mode_Post = Para_Catia.GenCATPart.Mode_Post;
Diameter_Fillet = Para_Catia.GenCATPart.Diameter_Fillet;
if Mode_Post == 1 || Mode_Post == 2
    Rate_Fluid = Para_Catia.GenCATPart.Rate_Fluid;
    Rate_BOI_Fluid = Para_Catia.GenCATPart.Rate_BOI_Fluid;
end
if Mode_Post == 1
    ScaleRate = Para_Catia.GenCATPart.ScaleRate;
    Index_NonActSection_ICEM = Para_Catia.GenCATPart.Index_NonActSection_ICEM;
    Index_Layer2Connect_ICEM = Para_Catia.GenCATPart.Index_Layer2Connect_ICEM;
end

Act_SaveStl = Para_Catia.GenCATPart.Act_SaveStl;
Act_SaveIgs = Para_Catia.GenCATPart.Act_SaveIgs;


%% 附加数据
%CatiaPara
Data{1,1} = [Dir_Catia,'\',FileName_CATPart];                                   
Data{1,2} = Act_Debug;
Data{1,3} = Act_Demo;
Data{1,4} = Act_SaveStl;
Data{1,5} = Act_SaveIgs;

%WingPara
Data{3,1} = length(Data_Section);                                           %截面数量
Data{3,2} = length(Data_Section(1).Up);                                     %翼型弦向坐标数量
Data{3,3} = join(string(Index_SurfConnect), ';');
Data{3,4} = Mode_GuideLine;

Data{4,1} = MAC_Length * 1000;                                              %MAC,m转成mm
Data{4,2} = MAC_x * 1000;

%PostPara
Data{6,1} = Mode_Post;
if Mode_Post == 0
    Data{6,2} = Diameter_Fillet;
    Data{6,3} = 999;
    Data{6,4} = join(string([999,999]), ';');
    Data{6,5} = 999;
    Data{6,6} = 999;
    Data{6,7} = 999;
end
if Mode_Post == 1
    Data{6,2} = Diameter_Fillet;
    Data{6,3} = Rate_Fluid;
    Data{6,4} = join(string([Rate_BOI_Fluid(1),999]), ';');
    Data{6,5} = ScaleRate;
    Data{6,6} = join(string(Index_NonActSection_ICEM), ';');
    Data{6,7} = join(string(Index_Layer2Connect_ICEM), ';');
end
if Mode_Post == 2
    Data{6,2} = Diameter_Fillet;
    Data{6,3} = Rate_Fluid;
    Data{6,4} = join(string([Rate_BOI_Fluid(1),Rate_BOI_Fluid(2)]), ';');
    Data{6,5} = 999;
    Data{6,6} = 999;
    Data{6,7} = 999;
end

%LeadPoint
for nSection = 1:length(Data_Airfoil)                                       %求前缘下表面曲率最大点
    if Mode_Post == 2 
        x_Raw = Data_Section(nSection).Low(:,1);
        z_Raw = Data_Section(nSection).Low(:,3);
        
        Length_Chord = abs(max(x_Raw) - min(x_Raw));
        Index = find(x_Raw <= min(x_Raw) + 0.2 * Length_Chord);
        x_Raw = x_Raw(Index);
        z_Raw = z_Raw(Index);

        x_Interp = linspace(x_Raw(1),x_Raw(end), 10000)';
        z_Interp = interp1(x_Raw,z_Raw, x_Interp, 'spline');

        dz1 = gradient(z_Interp) ./ gradient(x_Interp);
        dz2 = gradient(dz1) ./ gradient(x_Interp);
        Curv = abs(dz2) ./ ((1+dz1.^2).^(3/2));
        
        x_Interp = x_Interp(3:end-2);
        z_Interp = z_Interp(3:end-2);
        Curv = Curv(3:end-2);

        Index_Max = find(Curv==max(Curv),1);
        Data{6+nSection,1} = x_Interp(Index_Max);
        Data{6+nSection,2} = Data_Section(nSection).Low(1,2);
        Data{6+nSection,3} = z_Interp(Index_Max);
    else
        Data{6+nSection,1} = 999;
        Data{6+nSection,2} = 999;
        Data{6+nSection,3} = 999;
    end
end

Data{end+1,1}= [];


%% SectionPoint
%上表面
for nSection = 1:length(Data_Section)
    Num_Point = size(Data_Section(nSection).Up,1);
    Data(end + 1:end + Num_Point,1:3) = num2cell(Data_Section(nSection).Up);
end

%下表面
for nSection = length(Data_Section):-1:1
    Num_Point = size(Data_Section(nSection).Low,1);
    Data(end + 1:end + Num_Point,1:3) = num2cell(Data_Section(nSection).Low);
end


%% 写入Excel
Dir_Excel = [Dir_Catia,'\',FileName_CATPart,'.xls'];

if exist(Dir_Excel,'file')
    delete(Dir_Excel);
end
xlswrite(Dir_Excel, Data)


%% 运行VBA
if mod(nThread, 2) == 1
    CatiaName = 'CNEXT';
else
    CatiaName = 'DELMIA';
end

[Status,CMDOut] = dos(['cd /D ',Dir_Catia,'\','Template',' && ',...
                       'Set Dir_Excel_ForCatiaVBA=',Dir_Excel,' && ',...
                       'call ',CatiaName,' -batch -macro Catia_Template.catvba MainFun']);


%% 检查数据
Reboot = 0;

if Act_Debug == 0 && Mode_Post == 1
    
    Error_GenCATPart = 0;
    Error_WriteExcel = 0;
    
    %检查.CATPart以判定是否存在错误
    if Mode == 3
        Dir_CATPart_Error = [Dir_Catia,'\',FileName_CATPart,'.CATPart.Error'];
        if exist(Dir_CATPart_Error,'file')
            Error_GenCATPart = 1;
        end
    end
    
    %检查.xls以判定是否存在错误
    Data_Excel = readcell(Dir_Excel);
    Data_Excel{13,5} = 999;
    if isempty(Data_Excel{14,5}) || ismissing(Data_Excel{14,5})
        Error_WriteExcel = 1;
    end
    
    %生成重启符
    if Error_GenCATPart == 1 || Error_WriteExcel == 1
        Reboot = 1;
        
        %删除.CATPart.Error
        if Mode == 3
            Dir_CATPart_Error = [Dir_Catia,'\',FileName_CATPart,'.CATPart.Error'];
            if exist(Dir_CATPart_Error,'file')
                delete(Dir_CATPart_Error); 
            end
        end
    end
    
end


end














