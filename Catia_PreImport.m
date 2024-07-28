function Data_DesignVar = Catia_PreImport(Data_Import)


%% 数据校验
if sum(Data_Import.Root == 999) ~= 8
    if length(Data_Import.FlightState) ~= 3 || length(Data_Import.Root) ~= 13
        error('*** Error：Please Check Import ***')
    end

    if length(Data_Import.Root) ~= 13
        error('*** Error：Please Check Import ***')
    end
    for nPanel = 1:size(Data_Import.Panel,1)
        if length(Data_Import.Panel(nPanel,:)) ~= 13
            error('*** Error：Please Check Import ***')
        end
    end

    if ~isempty( find(isnan(Data_Import.Root), 1) ) || ~isempty( find(Data_Import.Root(6:end)/1e2>=1, 1) ) 
        error('*** Error：Please Check Import ***')
    end

    if ~isempty( find(isnan(Data_Import.Panel(:,1)), 1) ) ||...
       Data_Import.Panel(end,1)/1e5 >= 1
        error('*** Error：Please Check Import(展向位置) ***')
    end

    if isnan(Data_Import.Panel(end,2)) || Data_Import.Panel(end,2)/1e30 > 1
        error('*** Error：Please Check Import(梢弦长) ***')
    end

    if isnan(Data_Import.Panel(end,3))
        error('*** Error：Please Check Import(梢部扭转角) ***')
    end

    if isnan(Data_Import.Panel(end,4))
        error('*** Error：Please Check Import(前缘后掠角) ***')
    end

    if isnan(Data_Import.Panel(end,5))
        error('*** Error：Please Check Import(上反角) ***')
    end

    for nPanel = 1:size(Data_Import.Panel,1)
        if ~(sum(isnan( Data_Import.Panel(nPanel, 6:13) )) == 7 && abs( Data_Import.Panel(nPanel, 6)) == 1e30) &&...
           ~(sum(isnan( Data_Import.Panel(nPanel, 6:9) )) == 3 && abs( Data_Import.Panel(nPanel, 6)) == 1e30) &&...
           ~(sum(isnan( Data_Import.Panel(nPanel, 10:13) )) == 3 && abs( Data_Import.Panel(nPanel, 10)) == 1e30) &&...
           ~(sum(isnan( Data_Import.Panel(nPanel, 6:13) )) == 7 && abs( Data_Import.Panel(nPanel, 6)/1e2 ) >= 1) &&...
           isempty( find(Data_Import.Panel(nPanel, 6:13) < 1, 1) )
            error('*** Error：Please Check Import(翼型) ***')
        end
    end
end


%% 数据初始化
%飞行状态
Data_DesignVar.FlightState.Height = [];
Data_DesignVar.FlightState.Velocity = [];
Data_DesignVar.FlightState.Alpha = [];

%根弦
Data_DesignVar.Root.Dxyz = [];
Data_DesignVar.Root.Length = [];
Data_DesignVar.Root.Twist = [];

Data_DesignVar.Root.c1 = [];
Data_DesignVar.Root.c2 = [];
Data_DesignVar.Root.c3 = [];
Data_DesignVar.Root.c4 = [];

Data_DesignVar.Root.X_T = [];
Data_DesignVar.Root.T = [];
Data_DesignVar.Root.Rho0 = [];
Data_DesignVar.Root.Beta_TE = [];

%每个翼段
Num_Section = size(Data_Import.Panel,1);

Data_DesignVar.Panel(Num_Section).Span = [];
Data_DesignVar.Panel(Num_Section).Length = [];
Data_DesignVar.Panel(Num_Section).Twist = [];
Data_DesignVar.Panel(Num_Section).Sweepback = [];
Data_DesignVar.Panel(Num_Section).Dihedral = [];

Data_DesignVar.Panel(Num_Section).c1 = [];
Data_DesignVar.Panel(Num_Section).c2 = [];
Data_DesignVar.Panel(Num_Section).c3 = [];
Data_DesignVar.Panel(Num_Section).c4 = [];

Data_DesignVar.Panel(Num_Section).X_T = [];
Data_DesignVar.Panel(Num_Section).T = [];
Data_DesignVar.Panel(Num_Section).Rho0 = [];
Data_DesignVar.Panel(Num_Section).Beta_TE = [];


%% 飞行状态
Data_DesignVar.FlightState.Height = Data_Import.FlightState(1);
Data_DesignVar.FlightState.Velocity = Data_Import.FlightState(2);
Data_DesignVar.FlightState.Alpha = Data_Import.FlightState(3);


%% 根弦
Data_DesignVar.Root.Dxyz = [Data_Import.Root(1), Data_Import.Root(2), Data_Import.Root(3)];
Data_DesignVar.Root.Length = Data_Import.Root(4);
Data_DesignVar.Root.Twist = Data_Import.Root(5);

Data_DesignVar.Root.c1 = Data_Import.Root(6);
Data_DesignVar.Root.c2 = Data_Import.Root(7);
Data_DesignVar.Root.c3 = Data_Import.Root(8);
Data_DesignVar.Root.c4 = Data_Import.Root(9);

Data_DesignVar.Root.X_T = Data_Import.Root(10);
Data_DesignVar.Root.T = Data_Import.Root(11);
Data_DesignVar.Root.Rho0 = Data_Import.Root(12);
Data_DesignVar.Root.Beta_TE = Data_Import.Root(13);


%% 每个翼段_展向位置
Span(1) = 0;
for nPanel = 1:size(Data_Import.Panel,1)
    nSection = nPanel + 1;
    Span(nSection) = Data_Import.Panel(nPanel, 1);
end

%定位非全参数定义的翼段
nPanel_Auto_Main = [];
nPanel_Auto_Sub = [];
for nPanel = 1:size(Data_Import.Panel,1)
    if abs( Data_Import.Panel(nPanel,1)/1e5 ) >= 1 && abs( Data_Import.Panel(nPanel,1)/1e10 ) < 1
        nPanel_Auto_Main = [nPanel_Auto_Main, nPanel];
    end
    if abs( Data_Import.Panel(nPanel,1)/1e10 ) >= 1
        nPanel_Auto_Sub = [nPanel_Auto_Sub, nPanel];
    end
end
        
for nPanel = 1:size(Data_Import.Panel,1)
    
    nSection = nPanel + 1;
    
    %比例参数定义(优先)
    if ~isempty(find(nPanel_Auto_Main == nPanel, 1))
        
        i = nSection;
        j = nSection;
        while abs( Span(i)/1e5 ) >= 1
            i = i - 1;
            nStart = i;
        end
        while abs( Span(j)/1e5) >= 1
            j = j + 1;
            nEnd = j;
        end
            
        t_Span = Data_Import.Panel(nPanel, 1) - 1e5;
        Data_DesignVar.Panel(nPanel).Span = Span(nStart) + (Span(nEnd) - Span(nStart)) * t_Span;
        
        Span(nSection) = Data_DesignVar.Panel(nPanel).Span;
        
    end
    
end
 
for nPanel = 1:size(Data_Import.Panel,1)
    
    nSection = nPanel + 1;
    
    %比例参数定义(其余)
    if ~isempty(find(nPanel_Auto_Sub == nPanel, 1))
        
        i = nSection;
        j = nSection;
        while abs( Span(i)/1e5 ) >= 1
            i = i - 1;
            nStart = i;
        end
        while abs( Span(j)/1e5 ) >= 1
            j = j + 1;
            nEnd = j;
        end
            
        t_Span = Data_Import.Panel(nPanel, 1) - 1e10;
        Data_DesignVar.Panel(nPanel).Span = Span(nStart) + (Span(nEnd) - Span(nStart)) * t_Span;
        
        Span(nSection) = Data_DesignVar.Panel(nPanel).Span;
        
    %全参数定义
    else
        Data_DesignVar.Panel(nPanel).Span = Span(nSection);
    end
    
end


%% 每个翼段_梢弦长
Length(1) = Data_DesignVar.Root.Length;
for nPanel = 1:size(Data_Import.Panel,1)
    nSection = nPanel + 1;
    Length(nSection) = Data_Import.Panel(nPanel, 2);
end

%定位非全参数定义的翼段
nPanel_Auto = [];
for nPanel = 1:size(Data_Import.Panel,1)
    if isnan( Data_Import.Panel(nPanel,2) ) || abs( Data_Import.Panel(nPanel,2) ) == 1e30
        nPanel_Auto = [nPanel_Auto, nPanel];
    end
end

for nPanel = 1:size(Data_Import.Panel,1)
    
    nSection = nPanel + 1;
    
    if ~isempty(find(nPanel_Auto == nPanel, 1))
        %相同参数定义
        if abs( Data_Import.Panel(nPanel,2) ) == 1e30
            
            j = nSection;
            while isnan( Length(j) ) || abs( Length(j) ) == 1e30
                j = j + 1;
                nEnd = j;
            end
            
            Data_DesignVar.Panel(nPanel).Length = Length(nEnd);
            
            Length(nSection) = Data_DesignVar.Panel(nPanel).Length;
        end
        
    %全参数定义
    else
        Data_DesignVar.Panel(nPanel).Length = Length(nSection);
    end
    
end

%自动填充
for nPanel = 1:size(Data_Import.Panel,1)
    
    nSection = nPanel + 1;
    
    if ~isempty(find(nPanel_Auto == nPanel, 1))
        
        i = nSection;
        j = nSection;
        while isnan( Length(i) )
            i = i - 1;
            nStart = i;
        end
        while isnan( Length(j) )
            j = j + 1;
            nEnd = j;
        end
        
        t_Length = (Span(nSection) - Span(nStart))/(Span(nEnd) - Span(nStart));
        Data_DesignVar.Panel(nPanel).Length = Length(nStart) * (1-t_Length) + Length(nEnd) * t_Length;
        
        Length(nSection) = Data_DesignVar.Panel(nPanel).Length;
        
    end
    
end


%% 每个翼段_梢部扭转角
Twist(1) = Data_DesignVar.Root.Twist;
for nPanel = 1:size(Data_Import.Panel,1)
    nSection = nPanel + 1;
    Twist(nSection) = Data_Import.Panel(nPanel, 3);
end

%定位非全参数定义的翼段
nPanel_Auto = [];
for nPanel = 1:size(Data_Import.Panel,1)
    if isnan( Data_Import.Panel(nPanel,3) ) || abs( Data_Import.Panel(nPanel,3)/1e2 ) >= 1
        nPanel_Auto = [nPanel_Auto, nPanel];
    end
end

for nPanel = 1:size(Data_Import.Panel,1)
    
    nSection = nPanel + 1;
    
    if ~isempty(find(nPanel_Auto == nPanel, 1))
        %相同参数定义    
        if Data_Import.Panel(nPanel, 3)/1e30 == 1
            
            i = nSection;
            while abs( Twist(i)/1e2 ) >= 1 || isnan( Twist(i) )
                i = i - 1;
                nStart = i;
            end
            
            Data_DesignVar.Panel(nPanel).Twist = Twist(nStart);
            
            Twist(nSection) = Data_DesignVar.Panel(nPanel).Twist;
            
        %相似参数定义    
        elseif abs( Data_Import.Panel(nPanel, 3)/1e2 ) >= 1 && Data_Import.Panel(nPanel, 3)/1e30 ~= 1
            
            i = nSection;
            j = nSection;
            while abs( Twist(i)/1e2 ) >= 1 || isnan( Twist(i) )
                i = i - 1;
                nStart = i;
            end
            while abs( Twist(j)/1e2 ) >= 1 || isnan( Twist(j) )
                j = j + 1;
                nEnd = j;
            end
        
            t_Airfoil =  Data_Import.Panel(nPanel, 3) - 1e2;
            Data_DesignVar.Panel(nPanel).Twist = Twist(nStart) * (1-t_Airfoil) + Twist(nEnd) * t_Airfoil;
            
            Twist(nSection) = Data_DesignVar.Panel(nPanel).Twist;
         
        %自动填充
        elseif isnan( Data_Import.Panel(nPanel,3) )
            Data_DesignVar.Panel(nPanel).Twist = NaN;                       %在Catia_GenSection.m中计算扭转角
        end
    
    %全参数定义
    else
        Data_DesignVar.Panel(nPanel).Twist = Twist(nSection);
    end
    
end


%% 每个翼段_前缘后掠角
for nPanel = 1:size(Data_Import.Panel,1)
    Sweepback(nPanel) = Data_Import.Panel(nPanel, 4);
end

%定位非全参数定义的翼段
nPanel_Auto = [];
for nPanel = 1:size(Data_Import.Panel,1)
    if isnan( Data_Import.Panel(nPanel,4) )
        nPanel_Auto = [nPanel_Auto, nPanel];
    end
end

for nPanel = 1:size(Data_Import.Panel,1)
    
    %自动填充
    if ~isempty(find(nPanel_Auto == nPanel, 1))
        
        j = nPanel;
        while isnan( Sweepback(j) )
            j = j + 1;
            nEnd = j;
        end

        Data_DesignVar.Panel(nPanel).Sweepback = Sweepback(nEnd);
        
    %全参数定义
    else
        Data_DesignVar.Panel(nPanel).Sweepback = Sweepback(nPanel);
    end
    
end


%% 每个翼段_上反角
for nPanel = 1:size(Data_Import.Panel,1)
    Dihedral(nPanel) = Data_Import.Panel(nPanel, 5);
end

%定位非全参数定义的翼段
nPanel_Auto = [];
for nPanel = 1:size(Data_Import.Panel,1)
    if isnan( Data_Import.Panel(nPanel,5) )
        nPanel_Auto = [nPanel_Auto, nPanel];
    end
end

for nPanel = 1:size(Data_Import.Panel,1)
    
    %自动填充
    if ~isempty(find(nPanel_Auto == nPanel, 1))
        
        j = nPanel;
        while isnan( Dihedral(j) )
            j = j + 1;
            nEnd = j;
        end

        Data_DesignVar.Panel(nPanel).Dihedral = Dihedral(nEnd);
        
    %全参数定义
    else
        Data_DesignVar.Panel(nPanel).Dihedral = Dihedral(nPanel);
    end
    
end


%% 每个翼段_翼型
c1(1) = Data_DesignVar.Root.c1;
c2(1) = Data_DesignVar.Root.c2;
c3(1) = Data_DesignVar.Root.c3;
c4(1) = Data_DesignVar.Root.c4;

X_T(1) = Data_DesignVar.Root.X_T;
T(1) = Data_DesignVar.Root.T;
Rho0(1) = Data_DesignVar.Root.Rho0;
Beta_TE(1) = Data_DesignVar.Root.Beta_TE;

for nPanel = 1:size(Data_Import.Panel,1)
    nSection = nPanel + 1;
    
    c1(nSection) = Data_Import.Panel(nPanel, 6);
    c2(nSection) = Data_Import.Panel(nPanel, 7);
    c3(nSection) = Data_Import.Panel(nPanel, 8);
    c4(nSection) = Data_Import.Panel(nPanel, 9);
    
    X_T(nSection) = Data_Import.Panel(nPanel, 10);
    T(nSection) = Data_Import.Panel(nPanel, 11);
    Rho0(nSection) = Data_Import.Panel(nPanel, 12);
    Beta_TE(nSection) = Data_Import.Panel(nPanel, 13);
end

for nPanel = 1:size(Data_Import.Panel,1)
    
    nSection = nPanel + 1;
    
    %相同参数定义_弯度厚度均相同
    if sum(isnan( Data_Import.Panel(nPanel, 6:13) )) == 7 && abs( Data_Import.Panel(nPanel, 6)) == 1e30
        
        i = nSection;
        while isnan( c2(i) )
            i = i - 1;
            nStart = i;
        end
        
        Data_DesignVar.Panel(nPanel).c1 = c1(nStart);
        Data_DesignVar.Panel(nPanel).c2 = c2(nStart);
        Data_DesignVar.Panel(nPanel).c3 = c3(nStart);
        Data_DesignVar.Panel(nPanel).c4 = c4(nStart);
        
        Data_DesignVar.Panel(nPanel).X_T = X_T(nStart);
        Data_DesignVar.Panel(nPanel).T = T(nStart);
        Data_DesignVar.Panel(nPanel).Rho0 = Rho0(nStart);
        Data_DesignVar.Panel(nPanel).Beta_TE = Beta_TE(nStart);
        
    %相同参数定义_弯度相同
    elseif sum(isnan( Data_Import.Panel(nPanel, 6:9) )) == 3 && abs( Data_Import.Panel(nPanel, 6)) == 1e30

        i = nSection;
        while isnan( c2(i) )
            i = i - 1;
            nStart = i;
        end

        Data_DesignVar.Panel(nPanel).c1 = c1(nStart);
        Data_DesignVar.Panel(nPanel).c2 = c2(nStart);
        Data_DesignVar.Panel(nPanel).c3 = c3(nStart);
        Data_DesignVar.Panel(nPanel).c4 = c4(nStart);

        Data_DesignVar.Panel(nPanel).X_T = X_T(nSection);
        Data_DesignVar.Panel(nPanel).T = T(nSection);
        Data_DesignVar.Panel(nPanel).Rho0 = Rho0(nSection);
        Data_DesignVar.Panel(nPanel).Beta_TE = Beta_TE(nSection);
        
    %相同参数定义_厚度相同
    elseif sum(isnan( Data_Import.Panel(nPanel, 10:13) )) == 3 && abs( Data_Import.Panel(nPanel, 10)) == 1e30

        i = nSection;
        while isnan( T(i) )
            i = i - 1;
            nStart = i;
        end

        Data_DesignVar.Panel(nPanel).c1 = c1(nSection);
        Data_DesignVar.Panel(nPanel).c2 = c2(nSection);
        Data_DesignVar.Panel(nPanel).c3 = c3(nSection);
        Data_DesignVar.Panel(nPanel).c4 = c4(nSection);
        
        Data_DesignVar.Panel(nPanel).X_T = X_T(nStart);
        Data_DesignVar.Panel(nPanel).T = T(nStart);
        Data_DesignVar.Panel(nPanel).Rho0 = Rho0(nStart);
        Data_DesignVar.Panel(nPanel).Beta_TE = Beta_TE(nStart);
    
    %相似参数定义
    elseif sum(isnan( Data_Import.Panel(nPanel, 6:13) )) == 7 && abs( Data_Import.Panel(nPanel, 6)/1e2 ) >= 1
        
        t_Airfoil = Data_Import.Panel(nPanel, 6) - 1e2;
        
        Data_DesignVar.Panel(nPanel).c1 = c1(nSection-1) * (1-t_Airfoil) + c1(nSection+1) * t_Airfoil;
        Data_DesignVar.Panel(nPanel).c2 = c2(nSection-1) * (1-t_Airfoil) + c2(nSection+1) * t_Airfoil;
        Data_DesignVar.Panel(nPanel).c3 = c3(nSection-1) * (1-t_Airfoil) + c3(nSection+1) * t_Airfoil;
        Data_DesignVar.Panel(nPanel).c4 = c4(nSection-1) * (1-t_Airfoil) + c4(nSection+1) * t_Airfoil;
        
        Data_DesignVar.Panel(nPanel).X_T = X_T(nSection-1) * (1-t_Airfoil) + X_T(nSection+1) * t_Airfoil;
        Data_DesignVar.Panel(nPanel).T = T(nSection-1) * (1-t_Airfoil) + T(nSection+1) * t_Airfoil;
        Data_DesignVar.Panel(nPanel).Rho0 = Rho0(nSection-1) * (1-t_Airfoil) + Rho0(nSection+1) * t_Airfoil;
        Data_DesignVar.Panel(nPanel).Beta_TE = Beta_TE(nSection-1) * (1-t_Airfoil) + Beta_TE(nSection+1) * t_Airfoil;
        
    %全参数定义
    else
    
        Data_DesignVar.Panel(nPanel).c1 = c1(nSection);
        Data_DesignVar.Panel(nPanel).c2 = c2(nSection);
        Data_DesignVar.Panel(nPanel).c3 = c3(nSection);
        Data_DesignVar.Panel(nPanel).c4 = c4(nSection);
        
        Data_DesignVar.Panel(nPanel).X_T = X_T(nSection);
        Data_DesignVar.Panel(nPanel).T = T(nSection);
        Data_DesignVar.Panel(nPanel).Rho0 = Rho0(nSection);
        Data_DesignVar.Panel(nPanel).Beta_TE = Beta_TE(nSection);
        
    end
    
end


end















