%% 绘制甘特图
function [] = printgante( a,title_str )
    %a中每一列的格式：蒸馏塔/管道 油罐 开始时间 结束时间 原油类型
    w=0.5;       %横条宽度
    color=['r','g','b','c','m','y','w'];    %颜色

    dset = [];   %蒸馏塔/管道的集合
    %统计蒸馏塔/管道的个数
    for i=1:size(a,1)
        if ~ismember(a(i,1), dset)
            dset = [dset, a(i,1)];
        end
    end

    xlabel('time'); 
    ylabel('distillating tower/pipe line'); 
    title(title_str);
    axis([0 270 0 10]); 
    set(gca,'Box','on'); 
    set(gca,'XTick',0:30:240);      %设置刻度
    set(gca,'YTick',0:10/(size(dset,2)+1):10);      %设置刻度
    set(gca,'YTickLabel',{'';num2str((1:size(dset,2)-1)','DS%d');'Pipe';''});       %设置标记

    for ii=1:size(a,1) 
        %计算矩形区域的边界点
        x=a(ii,[3 3 4 4]);
        y=a(ii,1)*10/(size(dset,2)+1)+[-w/2 w/2 w/2 -w/2];
        %使用不同的颜色填充封闭的矩形区域
        if a(ii,5)==1
            p=patch('xdata',x,'ydata',y,'facecolor',color(1),'edgecolor','k'); 
        elseif a(ii,5)==2
            p=patch('xdata',x,'ydata',y,'facecolor',color(2),'edgecolor','k'); 
        elseif a(ii,5)==3
            p=patch('xdata',x,'ydata',y,'facecolor',color(3),'edgecolor','k'); 
        elseif a(ii,5)==4
            p=patch('xdata',x,'ydata',y,'facecolor',color(4),'edgecolor','k'); 
        elseif a(ii,5)==5
            p=patch('xdata',x,'ydata',y,'facecolor',color(5),'edgecolor','k'); 
        elseif a(ii,5)==6
            p=patch('xdata',x,'ydata',y,'facecolor',color(6),'edgecolor','k'); 
        elseif a(ii,5)==0
            p=patch('xdata',x,'ydata',y,'facecolor',color(7),'edgecolor','k'); 
        end
        %设置矩形区域的文字
        text(a(ii,3)+0.5,a(ii,1)*10/(size(dset,2)+1),[num2str(a(ii,2))]); 
    end 
end

