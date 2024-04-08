function  Data_Geometry = Catia_CalGeometry(Data_DesignVar)


%% 各截面数据
Num_Panel = length(Data_DesignVar.Panel);

%弦长
Length(1) = Data_DesignVar.Root.Length;
for nPanel = 1:Num_Panel
    nSection = nPanel + 1;
    Length(nSection) = Data_DesignVar.Panel(nPanel).Length;
end

%展向位置
Span(1) = 0;
for nPanel = 1:Num_Panel
    nSection = nPanel + 1;
    Span(nSection) = Data_DesignVar.Panel(nPanel).Span;
end

%后掠角
for nPanel = 1:Num_Panel
    Sweepback(nPanel) = deg2rad(Data_DesignVar.Panel(nPanel).Sweepback);
end

%x坐标(无扭转)
x_Section(1) = Data_DesignVar.Root.Dxyz(1);
for nPanel = 1:Num_Panel
    nSection = nPanel + 1;
    x_Section(nSection) = x_Section(nSection-1) + (Span(nSection)-Span(nSection-1)) * tan(Sweepback(nPanel));
end


%% 机翼面积
for nPanel = 1:Num_Panel
    nSection = nPanel + 1;
    Area(nPanel) = (Length(nSection) + Length(nSection-1))*(Span(nSection) - Span(nSection-1))/2;
end
Area_Total = sum(Area);


%% MAC_弦长
F_MAC_Chord = '( Length(1) + (Length(2)-Length(1))/(Span(2)-Span(1)) * (y-Span(1)) )  .*(y>=Span(1) & y<=Span(2))';

for nPanel = 2:Num_Panel
    nSection = nPanel + 1;
    F_MAC_Chord = [F_MAC_Chord, '+', '( Length(' num2str(nSection-1) ') + (Length(' num2str(nSection) ')-Length(' num2str(nSection-1) '))/(Span(' num2str(nSection) ')-Span(' num2str(nSection-1) ')) * (y-Span(' num2str(nSection-1) ')) )  .*(y>Span(' num2str(nSection-1) ') & y<=Span(' num2str(nSection) '))'];
end

eval(['F_MAC_Chord = @(y)', '(' F_MAC_Chord ').^2;']);
MAC_Length = 2/(Area_Total*2) * integral(F_MAC_Chord, 0,Span(end));


%% MAC_x坐标
F_MAC_x = ['( x_Section(1) + (x_Section(2)-x_Section(1))/(Span(2)-Span(1)) * (y-Span(1)) ) .* ',...
           '( Length(1) + (Length(2)-Length(1))/(Span(2)-Span(1)) * (y-Span(1)) )  .*(y>=Span(1) & y<=Span(2))'];
         
for nPanel = 2:Num_Panel
    nSection = nPanel + 1;
    F_MAC_x = [F_MAC_x, '+', '( x_Section(' num2str(nSection-1) ') + (x_Section(' num2str(nSection) ')-x_Section(' num2str(nSection-1) '))/(Span(' num2str(nSection) ')-Span(' num2str(nSection-1) ')) * (y-Span(' num2str(nSection-1) ')) ) .* ',...
                             '( Length(' num2str(nSection-1) ') + (Length(' num2str(nSection) ')-Length(' num2str(nSection-1) '))/(Span(' num2str(nSection) ')-Span(' num2str(nSection-1) ')) * (y-Span(' num2str(nSection-1) ')) )  .*(y>Span(' num2str(nSection-1) ') & y<=Span(' num2str(nSection) '))'];
end

eval(['F_MAC_x = @(y)', F_MAC_x, ';']);
MAC_x = 2/(Area_Total*2) * integral(F_MAC_x, 0,Span(end));


%% 等效锥形机翼
%判定是否需要等效
Act_Equal_Lead = 0;
x_Lead_B = [x_Section(1),x_Section(end)];
y_Lead_B = [Span(1),Span(end)];
for nSection = 2:length(x_Section)-1
    Delta_x = abs( x_Section(nSection) - interp1(y_Lead_B,x_Lead_B,Span(nSection)) );
    if Delta_x > 1e-4
        Act_Equal_Lead = 1;
    end
end

Act_Equal_Trail = 0;
x_Trail_B = [x_Section(1)+Length(1),x_Section(end)+Length(end)];
y_Trail_B = [Span(1),Span(end)];
for nSection = 2:length(x_Section)-1
    Delta_x = abs( x_Section(nSection) + Length(nSection) - interp1(y_Trail_B,x_Trail_B,Span(nSection)) );
    if Delta_x > 1e-4
        Act_Equal_Trail = 1;
    end
end

%求解
Options_Grad = optimoptions(@fmincon, 'Display','off','StepTolerance',1e-7);

if Act_Equal_Lead == 1
    for nRun = 1:3
        
        x_Trail = x_Section;
        y_Trail = Span;
        Max = max(max(x_Trail),max(y_Trail)) * nRun;
        
        InitValue = x_Section(1) - max(Length);
        Range = x_Section(1) + [-1,1]*2*max(Length);
        
        [x_Divide_Lead, fval] = fmincon(@(x_Divide) F_Err_S(x_Divide, x_Trail/Max, y_Trail/Max),...
                                        InitValue/Max, [],[],[],[], Range(1)/Max,Range(2)/Max, [], Options_Grad);
        
        x_Divide_Lead = x_Divide_Lead * Max;
        
        if fval <= 1e-5
            break;
        else
            if nRun == 3
                error('*** Error：Calculate Equivalent Planform ***');
            end
        end
        
    end
else
    x_Divide_Lead = x_Section(1);
end

if Act_Equal_Trail == 1
    for nRun = 1:3
        
        x_Trail = x_Section + Length;
        y_Trail = Span;
        Max = max(max(x_Trail),max(y_Trail)) * nRun;
        
        InitValue = x_Section(1) + Length(1) + max(Length);
        Range = x_Section(1) + Length(1) + [-1,1]*2*max(Length);
        
        [x_Divide_Trail, fval] = fmincon(@(x_Divide) F_Err_S(x_Divide, x_Trail/Max, y_Trail/Max),...
                                         InitValue/Max, [],[],[],[], Range(1)/Max,Range(2)/Max, [], Options_Grad);
        
        x_Divide_Trail = x_Divide_Trail * Max;
        
        if fval <= 1e-5
            break;
        else
            if nRun == 3
                error('*** Error：Calculate Equivalent Planform ***');
            end
        end
        
    end
else
    x_Divide_Trail = x_Section(1) + Length(1);
end

%赋值
x_Section_Equal = [x_Divide_Lead, x_Section(end)];
Length_Equal = [x_Divide_Trail-x_Divide_Lead, Length(end)];
Sweepback_Equal = repmat(rad2deg(atan( (x_Section_Equal(end)-x_Section_Equal(1))/Span(end) )),1,2);
Span_Equal = Span([1,end]);

S_Equal = (Length_Equal(1) + Length_Equal(end)) * Span_Equal(end)/2;


%% 1/4弦线后掠角
x_Root_Quarter = x_Section_Equal(1) + Length_Equal(1)/4;
x_Tip_Quarter = x_Section_Equal(end) + Length_Equal(end)/4;

Sweepback_Quarter = rad2deg(atan( (x_Tip_Quarter-x_Root_Quarter)/Span_Equal(end) ));


%% 赋值(单位：m & m^2)
%基于原始机翼
Data_Geometry.WingPlan.Length = Length;
Data_Geometry.WingPlan.Span = Span;
Data_Geometry.WingPlan.Sweepback = rad2deg(Sweepback);
Data_Geometry.WingPlan.x_Section = x_Section;

Data_Geometry.S = Area_Total / 1e6;
Data_Geometry.b = Span(end)*2 / 1e3;
Data_Geometry.AR = (Span(end)*2)^2/(Area_Total*2);

Data_Geometry.MAC_Length = MAC_Length / 1e3;
Data_Geometry.MAC_x = MAC_x / 1e3;

%基于等效机翼
Data_Geometry.WingPlan_Equal.Length = Length_Equal;
Data_Geometry.WingPlan_Equal.Span = Span_Equal;
Data_Geometry.WingPlan_Equal.Sweepback = rad2deg(Sweepback_Equal);
Data_Geometry.WingPlan_Equal.x_Section = x_Section_Equal;

Data_Geometry.S_Equal = S_Equal / 1e6;

Data_Geometry.TaperR = Length_Equal(end)/Length_Equal(1);
Data_Geometry.Sweepback_Quarter = Sweepback_Quarter;

%在Main_Catia.m中写入
Data_Geometry.Section = [];
Data_Geometry.Airfoil = [];
Data_Geometry.Diameter_Fillet = [];
Data_Geometry.Volume = [];
Data_Geometry.Rate = [];
Data_Geometry.ScaleRate = [];


end


%%
function [Err, S_1, S_2] = F_Err_S(x_Divide, x_Trail, y_Trail)

%简化后缘
F_Trail_Simple = @(y) x_Divide + (x_Trail(end)-x_Divide)/(y_Trail(end)-y_Trail(1)) .* y;

%求简化后缘与实际后缘的交点
nCorss = 0;
Coord_Corss = [];
for nSection = 1:length(y_Trail)-2
    F_Trail = ['(x_Trail(',num2str(nSection),') +',...
               '(x_Trail(',num2str(nSection+1),')-x_Trail(',num2str(nSection),'))/(y_Trail(',num2str(nSection+1),')-y_Trail(',num2str(nSection),')) .* (y-y_Trail(',num2str(nSection),')))',...
               '.*(y>y_Trail(',num2str(nSection),') & y<=y_Trail(',num2str(nSection+1),'))'];
    F_Trail = eval(['@(y)', F_Trail]);
    
    Str_F_Total = [char(vpa(F_Trail_Simple)),'-',char(vpa(F_Trail))];
    Str_F_Total = replace(Str_F_Total,'*','.*');
    F = eval(['@(y) abs(',Str_F_Total,')']);
    
    Options_Grad = optimoptions(@fmincon, 'Display','off','StepTolerance',1e-7);
    InitValue = y_Trail(nSection);
    Range = [y_Trail(nSection), y_Trail(nSection+1)];
    [y, fval] = fmincon(F, InitValue, [],[],[],[], Range(1),Range(2), [], Options_Grad);
    
    if fval < 1e-4 && y > max(1e-4,1.001*Range(1)) && y < 0.999*Range(2)
        nCorss = nCorss + 1;
        Coord_Corss(nCorss,2) = y;
        Coord_Corss(nCorss,1) = F_Trail(y);
    end
end

%排列
Coord_Total(:,1) = x_Trail;
Coord_Total(:,2) = y_Trail;
if abs(x_Divide - Coord_Total(1,1)) > 1e-4
    Coord_Total = [[x_Divide,y_Trail(1)]; Coord_Total];
end

if ~isempty(Coord_Corss)
    for nCross = 1:size(Coord_Corss,1)
        Index = find(Coord_Total(:,2) > Coord_Corss(nCross,2),1);
        
        Num_Row = size(Coord_Total,1);
        Coord_Total([1:Index-1,Index+1:Num_Row+1],:) = Coord_Total;
        
        Coord_Total(Index,:) = Coord_Corss(nCross,:);
    end
    
    for nCross = 1:size(Coord_Corss,1)
        Index_Cross(nCross) = find(Coord_Total(:,2) == Coord_Corss(nCross,2));
    end
    Index_Cross = [1,Index_Cross,size(Coord_Total,1)];
    
else
    Index_Cross = [1,size(Coord_Total,1)];
end

%求面积差
S_1 = 0;
S_2 = 0;
for nPoly = 1:length(Index_Cross)-1
    x_Poly = Coord_Total(Index_Cross(nPoly):Index_Cross(nPoly+1),1);
    y_Poly = Coord_Total(Index_Cross(nPoly):Index_Cross(nPoly+1),2);
    
    warning off;
    Poly = polyshape(x_Poly, y_Poly);
    warning on;

    if mod(nPoly,2) == 1
        S_1 = S_1 + Poly.area;
    else
        S_2 = S_2 + Poly.area;
    end
end

Err = abs(S_1 - S_2);

end




























