function [  ] = Draw_bubble( data )
%% 绘制气泡图
    x = data(:,1)';
    y = data(:,2)';
    z = data(:,3)';
    
    %计算能耗等级
 	energe_level = unique(data(:,5));
    level = zeros(1,size(data,1));
    for i=1:size(data,1)
        for j=1:length(energe_level)
            if data(i,5) == energe_level(j)
                k=j;
                break;
            end
        end
        level(i) = k;
    end
    sz_base = 90;
    u = 50 + sz_base * level;
    
    %生成颜色
 	gnums = unique(data(:,4));
	c = parula(length(gnums));    
    c = fliplr(c')';
    
    %根据使用的罐的个数绘制图像
    figure;
    g_num = data(:,4);
    legendtext = cell(1,length(gnums));
    for i=1:length(gnums)
        ind = find(g_num== gnums(i));  %供油罐个数为tmp(i)的记录下标
        scatter3(x(ind),y(ind),z(ind),u(ind),c(i,:),'fill');%颜色是用罐个数，气泡大小是切换次数
        legendtext{i} = num2str(gnums(i));
        hold on;
    end
    legend(legendtext);
    title('Visible of Prato front');    
    range = [min(x), max(x), min(y), max(y), min(z), max(z)];
    axis(range);% 坐标轴的显示范围
    xlabel('Pipeline mixing cost');% x轴名称、字体及大小
    ylabel('Tank bottom mixing cost');% y轴名称、字体及大小
    zlabel('Number of switchs');% z轴名称、字体及大小
%    text(max(x), min(y), max(z),'气泡的面积越大，代表能耗越高','fontsize',12);
%     colorbar;%显示colorbar
%     caxis([6.5,10.5]);%设置colorbar的范围
%     colormap(c);%设置colorbar的颜色
end