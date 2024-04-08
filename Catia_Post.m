function Catia_Post(Dir, Para_Catia, Data_DesignVar, Data_Geometry, nCatia, Mode, Act_Save)


Dir_Catia = Dir.Catia;

FileName_Save_CATPart = Para_Catia.Common.FileName_Save_CATPart;

if Mode == 1
    FileName_Catia = FileName_Save_CATPart;
else
    FileName_Catia = [FileName_Save_CATPart,'_',num2str(nCatia,'%04d')];
end

FontSize = 10;
LineWidth = 1.6;

Data_FigurePara = importdata('FigurePara.txt');
RangeLimit_x_Total = Data_FigurePara(1,:);
RangeLimit_y_Total = Data_FigurePara(2,:);
RangeLimit_z_Total = Data_FigurePara(3,:);
Scale_Length_Twist = Data_FigurePara(4,1);
Scale_Length_Dihedral = Data_FigurePara(5,1);
Scale_Angle_Twist = Data_FigurePara(6,1);
Scale_Angle_Dihedral = Data_FigurePara(7,1);


%% 读取数据
Data_Stl_Total = stlread([Dir_Catia,'\',FileName_Catia,'.stl']);
Data_Stl.Points = Data_Stl_Total.Points;                              
Data_Stl.ConnectivityList = Data_Stl_Total.ConnectivityList;


%% 作图准备
hFigure = figure('NumberTitle','off', 'Name',FileName_Catia);
FigureSize = [1400,940];
set(hFigure, 'Position',[(1920-FigureSize(1))/2,(1200-FigureSize(2))/2,FigureSize]);
hAxis = tight_subplot(1, 1, 0, 0.07, 0.07*FigureSize(2)/FigureSize(1));
set(hAxis, 'FontSize',FontSize);
axis equal;
hold on;
grid on;
view(3);

hTextBox = annotation('textbox', [0.01,0.888,0.1,0.1],...
                      'String',replace(FileName_Catia,'_','\_'),...
                      'EdgeColor','none', 'FontSize',FontSize+1);
hTextBox.HorizontalAlignment = 'Left';
hTextBox.VerticalAlignment = 'Top';

xlim(RangeLimit_x_Total);
ylim(RangeLimit_y_Total);
zlim(RangeLimit_z_Total);


%% 作图_基准
%绘制扭转角
x_End_Total = max(RangeLimit_x_Total) - Scale_Length_Twist * (max(RangeLimit_x_Total)-min(RangeLimit_x_Total));

for nSection = 1:length(Data_Geometry.Section)
    Length = Scale_Length_Twist * (max(RangeLimit_x_Total)-min(RangeLimit_x_Total));
    Twist = 0;

    XYZ_Center = [x_End_Total, Data_Geometry.Section(nSection).Up(end,2), 0];
    XYZ(1,:) = XYZ_Center - [1/4*Length, 0, -1/4*Length*sin(Twist)];
    XYZ(2,:) = XYZ_Center + [3/4*Length, 0, -3/4*Length*sin(Twist)];
    plot3(XYZ(:,1), -XYZ(:,2), XYZ(:,3), 'k:','LineWidth',LineWidth);
end

%绘制下反角
for nSection = 1:length(Data_Geometry.Section)-1
    Span = Scale_Length_Dihedral * (max(RangeLimit_y_Total)-min(RangeLimit_y_Total));
    Dihedral = 0;

    XYZ(1,:) = [x_End_Total, Data_Geometry.Section(nSection).Up(end,2), 0];
    XYZ(2,:) = XYZ(1,:) + [0, Span, Span*sin(Dihedral)];
    plot3(XYZ(:,1), -XYZ(:,2), XYZ(:,3), 'k:','LineWidth',LineWidth);
end


%% 作图_当前
%绘制曲面
patch('Faces',Data_Stl.ConnectivityList, 'Vertices',Data_Stl.Points.* [1,-1,1],...
      'Edgecolor','none', 'Facecolor',[220,220,220]/255,...
      'FaceAlpha',0.4);
camlight('headlight');
material('dull');

%绘制翼型
for nSection = 1:length(Data_Geometry.Section)
    XYZ_Up = Data_Geometry.Section(nSection).Up;
    XYZ_Low = Data_Geometry.Section(nSection).Low;
    
    plot3(XYZ_Up(:,1), -XYZ_Up(:,2), XYZ_Up(:,3), 'k-','LineWidth',LineWidth);
    plot3(XYZ_Low(:,1), -XYZ_Low(:,2), XYZ_Low(:,3), 'k-','LineWidth',LineWidth)
end

%绘制扭转角
for nSection = 1:length(Data_Geometry.Section)
    Length = Scale_Length_Twist * (max(RangeLimit_x_Total)-min(RangeLimit_x_Total));
    if nSection == 1
        Twist = Scale_Angle_Twist * deg2rad(Data_DesignVar.Root.Twist);
    else
        nPanel = nSection - 1;
        Twist = Scale_Angle_Twist * deg2rad(Data_DesignVar.Panel(nPanel).Twist);
    end

    XYZ_Center = [x_End_Total, Data_Geometry.Section(nSection).Up(end,2), 0];
    XYZ(1,:) = XYZ_Center - [1/4*Length, 0, -1/4*Length*sin(Twist)];
    XYZ(2,:) = XYZ_Center + [3/4*Length, 0, -3/4*Length*sin(Twist)];
    plot3(XYZ(:,1), -XYZ(:,2), XYZ(:,3), 'k-','LineWidth',LineWidth);
end

%绘制下反角
for nSection = 1:length(Data_Geometry.Section)-1
    Span = Scale_Length_Dihedral * (max(RangeLimit_y_Total)-min(RangeLimit_y_Total));
    Dihedral = Scale_Angle_Dihedral * deg2rad(Data_DesignVar.Panel(nSection).Dihedral);

    XYZ(1,:) = [x_End_Total, Data_Geometry.Section(nSection).Up(end,2), 0];
    XYZ(2,:) = XYZ(1,:) + [0, Span, Span*sin(Dihedral)];
    plot3(XYZ(:,1), -XYZ(:,2), XYZ(:,3), 'k-','LineWidth',LineWidth);
end

%绘制俯视图
Length = Data_Geometry.WingPlan.Length;
Span = Data_Geometry.WingPlan.Span;
x_Section = Data_Geometry.WingPlan.x_Section;
Num_Panel = length(Span)-1;

z = [RangeLimit_z_Total(1), RangeLimit_z_Total(1)];

x = [x_Section(1), x_Section(1) + Length(1)];                               %根弦
y = [Span(1), Span(1)];
plot3(x,-y,z,'k-','LineWidth',LineWidth);

for nPanel = 1:Num_Panel                                                    %各个翼段
    nSection = nPanel + 1;

    %前缘
    x = [x_Section(nSection-1), x_Section(nSection)];
    y = [Span(nSection-1), Span(nSection)];
    plot3(x,-y,z,'k-','LineWidth',LineWidth);

    %后缘
    x = [x_Section(nSection) + Length(nSection), x_Section(nSection-1) + Length(nSection-1)];
    y = [Span(nSection), Span(nSection-1)];
    plot3(x,-y,z,'k-','LineWidth',LineWidth);

    %梢弦
    x = [x_Section(nSection), x_Section(nSection) + Length(nSection)];
    y = [Span(nSection), Span(nSection)];
    plot3(x,-y,z,'k-','LineWidth',LineWidth);
end


%% 保存
if Act_Save == 1
    %保存.fig
    % saveas(hFigure, [Dir_Catia,'\',FileName_Catia,'.fig']);

    %保存.png
    print(hFigure, [Dir_Catia,'\',FileName_Catia,'.png'],'-dpng','-r150');

    %删除.fig
    % Dir_Fig = [Dir_Catia,'\',FileName_Catia,'.fig'];
    % if exist(Dir_Fig,'file')
    %     dos(['del ',Dir_Fig]);
    % end
    
    %关闭figure
    close(hFigure);
end


end















