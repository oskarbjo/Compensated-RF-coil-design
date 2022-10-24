function uiPlot(DATA1,DATA2,DATA3,DATA4,DATA5,DATA6)
    fig=uifigure;
    N=191;
    ax = uiaxes(fig);
    plot(ax,DATA1(1,:),DATA4(1,:));
    hold(ax,'on');
    plot(ax,DATA2(1,:),DATA5(1,:));
    plot(ax,DATA3(1,:),DATA6(1,:));
    sld=uislider(fig,'Limits',[1,300],'Value',191,'ValueChangedFcn',@(sld,event) updatePlot(sld,ax,DATA1,DATA2,DATA3,DATA4,DATA5,DATA6));
    sld.Value=191;

    
end

function updatePlot(sld,ax,DATA1,DATA2,DATA3,DATA4,DATA5,DATA6)
    sld.Value = round(sld.Value)
    plot(ax,DATA1(sld.Value,:),DATA4(sld.Value,:));
    hold(ax,'on');
    plot(ax,DATA2(sld.Value,:),DATA5(sld.Value,:));
    plot(ax,DATA3(sld.Value,:),DATA6(sld.Value,:));
    hold(ax,'off');
end