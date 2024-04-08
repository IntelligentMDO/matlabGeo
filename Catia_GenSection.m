function [Data_Section, Data_Airfoil, Data_DesignVar] = Catia_GenSection(Para_Catia, Data_DesignVar)


Density_Airfoil = Para_Catia.GenSection.Density_Airfoil;
Diameter_Fillet = Para_Catia.GenSection.Diameter_Fillet;


%% 设定每个截面弦向坐标密度
% x0 = 1-cos(linspace(0,1,Density_Airfoil).*pi);
% x0 = x0./(max(x0));

%定义分界点及基函数
t = 0.3;
F1 = @(x) x.^4;
F2 = @(x) 1-cos(x.*pi);

%求解1阶导数相同点
x_Total = linspace(0,1,1000);
y_F1 = F1(x_Total);
y_F2 = F2(x_Total);
dy1_F1 = gradient(y_F1)./gradient(x_Total);
dy1_F2 = gradient(y_F2)./gradient(x_Total);

F_Err = @(x) abs( interp1(x_Total, dy1_F1, t, 'spline') - interp1(x_Total, dy1_F2, x, 'spline') );
InitValue = 0;
Range = [0,0.5];
Options_Grad = optimoptions(@fmincon, 'Display','off','StepTolerance',1e-7);
A = fmincon(F_Err, InitValue, [],[],[],[], Range(1),Range(2), [], Options_Grad);

%整合函数
F3 = @(x) F2(A+x*(1-A));
F_Total_Airfoil = @(x) F1(x).*(x>=0 & x<=t) + (F3((x-t)/(1-t))-F3(0)+F1(t)).*(x>t & x<=1);

%生成弦向坐标密度
x0 = F_Total_Airfoil(linspace(0,1,Density_Airfoil));
x0 = x0./(max(x0));


%% 根弦
Length(1) = Data_DesignVar.Root.Length;
Twist(1) = deg2rad(Data_DesignVar.Root.Twist);

c1 = Data_DesignVar.Root.c1;
c2 = Data_DesignVar.Root.c2;
c3 = Data_DesignVar.Root.c3;
c4 = Data_DesignVar.Root.c4;

X_T = Data_DesignVar.Root.X_T;
T = Data_DesignVar.Root.T;
Rho0 = Data_DesignVar.Root.Rho0;
Beta_TE = Data_DesignVar.Root.Beta_TE;

%求解翼型
z_Camber_Basic = Catia_IGP_Camber(c1, c2, c3, c4, x0);

z_TE = Diameter_Fillet / Length(1);
[z_Up_Basic, z_Low_Basic] = Catia_IGP_Thickness(X_T, T, Rho0, Beta_TE, z_TE, x0, z_Camber_Basic);

Data_Airfoil(1).x = x0';
Data_Airfoil(1).Camber = z_Camber_Basic';
Data_Airfoil(1).Up = z_Up_Basic';
Data_Airfoil(1).Low = z_Low_Basic';

%生成截面三维坐标(无修正)
x_Local = x0 * Length(1);
y_Local = zeros(1,length(x0));
z_Camber_Local = z_Camber_Basic * Length(1);
z_Up_Local = z_Up_Basic * Length(1);
z_Low_Local = z_Low_Basic * Length(1);

%坐标修正_翼根位置
Dxyz(1).Root = Data_DesignVar.Root.Dxyz;

%坐标修正_扭转
Axis = 1/4 * Length(1);
Ly = [cos(Twist(1))    0    -sin(Twist(1));
      0                1     0
      sin(Twist(1))    0     cos(Twist(1))];

for n = 1:length(x_Local)
    Dxyz(1).Twist.Camber(n,:) = ( Ly'*[x_Local(n)-Axis; 0; z_Camber_Local(n)] - [x_Local(n)-Axis; 0; z_Camber_Local(n)] )';
    Dxyz(1).Twist.Up(n,:) = ( Ly'*[x_Local(n)-Axis; 0; z_Up_Local(n)] - [x_Local(n)-Axis; 0; z_Up_Local(n)] )';
    Dxyz(1).Twist.Low(n,:) = ( Ly'*[x_Local(n)-Axis; 0; z_Low_Local(n)] - [x_Local(n)-Axis; 0; z_Low_Local(n)] )';
end

Data_Section(1).Camber = [x_Local', y_Local', z_Camber_Local'] + Dxyz(1).Root + Dxyz(1).Twist.Camber;
Data_Section(1).Up = [x_Local', y_Local', z_Up_Local'] + Dxyz(1).Root + Dxyz(1).Twist.Up;
Data_Section(1).Low = [x_Local', y_Local', z_Low_Local'] + Dxyz(1).Root + Dxyz(1).Twist.Low;


%% 每个翼段
nSection_Auto = [];

for nPanel = 1:length(Data_DesignVar.Panel)
    
    nSection = nPanel + 1;
    
    Span(nSection) = Data_DesignVar.Panel(nPanel).Span;
    Length(nSection) = Data_DesignVar.Panel(nPanel).Length;
    Sweepback = deg2rad(Data_DesignVar.Panel(nPanel).Sweepback);
    Dihedral = deg2rad(Data_DesignVar.Panel(nPanel).Dihedral);
    Twist(nSection) = deg2rad(Data_DesignVar.Panel(nPanel).Twist);
    
    c1 = Data_DesignVar.Panel(nPanel).c1;
    c2 = Data_DesignVar.Panel(nPanel).c2;
    c3 = Data_DesignVar.Panel(nPanel).c3;
    c4 = Data_DesignVar.Panel(nPanel).c4;
    
    X_T = Data_DesignVar.Panel(nPanel).X_T;
    T = Data_DesignVar.Panel(nPanel).T;
    Rho0 = Data_DesignVar.Panel(nPanel).Rho0;
    Beta_TE = Data_DesignVar.Panel(nPanel).Beta_TE;

    %求解翼型
    z_Camber_Basic = Catia_IGP_Camber(c1, c2, c3, c4, x0);
    
    z_TE = Diameter_Fillet / Length(nSection);
    [z_Up_Basic, z_Low_Basic] = Catia_IGP_Thickness(X_T, T, Rho0, Beta_TE, z_TE, x0, z_Camber_Basic);
    
    Data_Airfoil(nSection).x = x0';
    Data_Airfoil(nSection).Camber = z_Camber_Basic';
    Data_Airfoil(nSection).Up = z_Up_Basic';
    Data_Airfoil(nSection).Low = z_Low_Basic';
    
    %生成截面三维坐标(无修正)
    x_Local = x0 * Length(nSection);
    y_Local = ones(1,length(x0)) * Span(nSection);
    z_Camber_Local = z_Camber_Basic * Length(nSection);
    z_Up_Local = z_Up_Basic * Length(nSection);
    z_Low_Local = z_Low_Basic * Length(nSection);
        
    %坐标修正_翼根位置
    Dxyz(nSection).Root = Dxyz(1).Root;
    
    %坐标修正_后掠
    if nPanel < 2
        Dxyz(nSection).Sweepback = [Span(nSection)*tan(Sweepback), 0, 0];
    else
        Dxyz(nSection).Sweepback = Dxyz(nSection-1).Sweepback + [(Span(nSection)-Span(nSection-1))*tan(Sweepback), 0, 0];
    end
    
    %坐标修正_下反角
    if nPanel < 2
        Dxyz(nSection).Dihedral = [0, 0, Span(nSection)*tan(Dihedral)];
    else
        Dxyz(nSection).Dihedral = Dxyz(nSection-1).Dihedral + [0, 0, (Span(nSection)-Span(nSection-1))*tan(Dihedral)];
    end
    
    %坐标修正_扭转_全参数定义
    if ~isnan(Twist(nSection))
        Axis = 1/4 * Length(nSection);
        Ly = [cos(Twist(nSection))    0    -sin(Twist(nSection));
              0                       1     0
              sin(Twist(nSection))    0     cos(Twist(nSection))];
        
        for n = 1:length(x_Local)
            Dxyz(nSection).Twist.Camber(n,:) = ( Ly'*[x_Local(n)-Axis; 0; z_Camber_Local(n)] - [x_Local(n)-Axis; 0; z_Camber_Local(n)] )';
            Dxyz(nSection).Twist.Up(n,:) = ( Ly'*[x_Local(n)-Axis; 0; z_Up_Local(n)] - [x_Local(n)-Axis; 0; z_Up_Local(n)] )';
            Dxyz(nSection).Twist.Low(n,:) = ( Ly'*[x_Local(n)-Axis; 0; z_Low_Local(n)] - [x_Local(n)-Axis; 0; z_Low_Local(n)] )';
        end
    else
        nSection_Auto = [nSection_Auto, nSection];
        Dxyz(nSection).Twist.Camber = [0 0 0];
        Dxyz(nSection).Twist.Up = [0 0 0];
        Dxyz(nSection).Twist.Low = [0 0 0];
    end
  
    Data_Section(nSection).Camber = [x_Local', y_Local', z_Camber_Local']...
                                    + Dxyz(nSection).Root + Dxyz(nSection).Sweepback + Dxyz(nSection).Dihedral + Dxyz(nSection).Twist.Camber;
                                
    Data_Section(nSection).Up = [x_Local', y_Local', z_Up_Local']...
                                + Dxyz(nSection).Root + Dxyz(nSection).Sweepback + Dxyz(nSection).Dihedral + Dxyz(nSection).Twist.Up;
                            
    Data_Section(nSection).Low = [x_Local', y_Local', z_Low_Local']...
                                 + Dxyz(nSection).Root + Dxyz(nSection).Sweepback + Dxyz(nSection).Dihedral + Dxyz(nSection).Twist.Low;
    
end


%% 坐标修正_扭转_自动定义
for nSection = nSection_Auto
    
    nPanel = nSection - 1;
    
    %生成截面三维坐标(无修正)
    x0 = Data_Airfoil(nSection).x';
    z_Camber_Basic = Data_Airfoil(nSection).Camber';
    z_Up_Basic = Data_Airfoil(nSection).Up';
    z_Low_Basic = Data_Airfoil(nSection).Low';
    
    x_Local = x0 * Length(nSection);
    z_Camber_Local = z_Camber_Basic * Length(nSection);
    z_Up_Local = z_Up_Basic * Length(nSection);
    z_Low_Local = z_Low_Basic * Length(nSection);
    
    %计算扭转角
    i = nSection;
    j = nSection;
    while isnan( Twist(i) )
        i = i - 1;
        nStart = i;
    end
    while isnan( Twist(j) )
        j = j + 1;
        nEnd = j;
    end
        
    Coord_In_Lead = [Data_Section(nStart).Camber(1,1), Data_Section(nStart).Camber(1,3)];
    Coord_In_Trail = [Data_Section(nStart).Camber(end,1), Data_Section(nStart).Camber(end,3)];
    
    Coord_Out_Lead = [Data_Section(nEnd).Camber(1,1), Data_Section(nEnd).Camber(1,3)];
    Coord_Out_Trail = [Data_Section(nEnd).Camber(end,1), Data_Section(nEnd).Camber(end,3)];
    
    t_Length = (Span(nSection) - Span(nStart))/(Span(nEnd) - Span(nStart));

    Coord_Mid_Lead = Coord_In_Lead * (1-t_Length) + Coord_Out_Lead * t_Length;
    Coord_Mid_Trail = Coord_In_Trail * (1-t_Length) + Coord_Out_Trail * t_Length;
    
    Twist(nSection) = -atan( (Coord_Mid_Trail(2)-Coord_Mid_Lead(2))/(Coord_Mid_Trail(1)-Coord_Mid_Lead(1)) );
    Data_DesignVar.Panel(nPanel).Twist = rad2deg(Twist(nSection));
    
    Axis = 1/4 * Length(nSection);
    Ly = [cos(Twist(nSection))    0    -sin(Twist(nSection));
          0                       1     0
          sin(Twist(nSection))    0     cos(Twist(nSection))];
    
    for n = 1:length(x_Local)
        Dxyz(nSection).Twist.Camber(n,:) = ( Ly'*[x_Local(n)-Axis; 0; z_Camber_Local(n)] - [x_Local(n)-Axis; 0; z_Camber_Local(n)] )';
        Dxyz(nSection).Twist.Up(n,:) = ( Ly'*[x_Local(n)-Axis; 0; z_Up_Local(n)] - [x_Local(n)-Axis; 0; z_Up_Local(n)] )';
        Dxyz(nSection).Twist.Low(n,:) = ( Ly'*[x_Local(n)-Axis; 0; z_Low_Local(n)] - [x_Local(n)-Axis; 0; z_Low_Local(n)] )';
    end
    
    Data_Section(nSection).Camber = Data_Section(nSection).Camber + Dxyz(nSection).Twist.Camber;
    Data_Section(nSection).Up = Data_Section(nSection).Up + Dxyz(nSection).Twist.Up;
    Data_Section(nSection).Low = Data_Section(nSection).Low + Dxyz(nSection).Twist.Low;
    
end


end











