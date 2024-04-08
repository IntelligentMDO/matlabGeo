function  Catia_DrawWingPlan(Data_Import, Data_Geometry, nCatia)


hFigure = figure('NumberTitle','off', 'Name',num2str(nCatia));
FigureSize = [940,940];
set(hFigure, 'Position',[(1920-FigureSize(1))/2,(1200-FigureSize(2))/2,FigureSize]);
hAxis = tight_subplot(1, 1, 0, 0.07, 0.07*FigureSize(2)/FigureSize(1));


%% 平面图
LineStyle_Total = {'-','-.'};
LineWidth_Total = [1.6,1.6];

for nDraw = 1:2
    %赋值
    if nDraw == 1
        Length = Data_Geometry.WingPlan.Length;
        Span = Data_Geometry.WingPlan.Span;
        x_Section = Data_Geometry.WingPlan.x_Section;
    else
        Length = Data_Geometry.WingPlan_Equal.Length;
        Span = Data_Geometry.WingPlan_Equal.Span;
        x_Section = Data_Geometry.WingPlan_Equal.x_Section;
    end
    Num_Panel = length(Span)-1;

    %根弦
    x = [x_Section(1), x_Section(1) + Length(1)];
    y = [Span(1), Span(1)];
    plot(x, y, 'k', 'LineWidth',LineWidth_Total(nDraw));
    hold on;
    
    %各个翼段
    for nPanel = 1:Num_Panel
        nSection = nPanel + 1;

        %前缘
        x = [x_Section(nSection-1), x_Section(nSection)];
        y = [Span(nSection-1), Span(nSection)];
        plot(x, y, ['k',LineStyle_Total{nDraw}], 'LineWidth',LineWidth_Total(nDraw));

        %后缘
        x = [x_Section(nSection) + Length(nSection), x_Section(nSection-1) + Length(nSection-1)];
        y = [Span(nSection), Span(nSection-1)];
        plot(x, y, ['k',LineStyle_Total{nDraw}], 'LineWidth',LineWidth_Total(nDraw));

        %梢弦
        x = [x_Section(nSection), x_Section(nSection) + Length(nSection)];
        y = [Span(nSection), Span(nSection)];

        if nPanel < Num_Panel
            if sum(isnan( Data_Import.Panel(nPanel, 6:13) )) == 7
                plot(x, y, 'k--', 'LineWidth',LineWidth_Total(nDraw));
            else
                plot(x, y, ['k',LineStyle_Total{nDraw}], 'LineWidth',LineWidth_Total(nDraw));
            end
        else
            plot(x, y, ['k',LineStyle_Total{nDraw}], 'LineWidth',LineWidth_Total(nDraw));
        end
    end
end

%MAC
x = [Data_Geometry.MAC_x, Data_Geometry.MAC_x + Data_Geometry.MAC_Length] * 1000;
y = [0, 0];
plot(x, y, 'LineWidth',5);

axis equal;
grid on;


%% 平面几何参数
Message={['              S：',num2str(Data_Geometry.S,'%.4f'),' (',num2str(Data_Geometry.S_Equal,'%.4f'),')'];...
         ['              b：',num2str(Data_Geometry.b,'%.2f')];...
         ['           AR：',num2str(Data_Geometry.AR,'%.2f')];...
         ['        MAC：',num2str(Data_Geometry.MAC_Length,'%.4f')];
         ['    TaperR：',num2str(Data_Geometry.TaperR,'%.2f')];
         ['Sweep1/4：',num2str(Data_Geometry.Sweepback_Quarter,'%.1f')]};

annotation('textbox', [0.15,0.8,0.1,0.1], 'String',Message, 'LineStyle','none', 'HorizontalAlignment','left');

drawnow;


end












