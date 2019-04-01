%% 显示并保存PF结果
function [] = show_all_results(  )
    problem = 'DMOPs1_OS';                  %问题
    title = {'管道混合成本','罐底混合成本','供油罐切换次数','供油罐个数'};
    filepath = sprintf('F:\\多目标优化论文\\实验——代码\\PlatEMO v1.5 (2017-12)\\Data/%s/%s.xls',problem,problem);
    PF=load(sprintf('F:\\多目标优化论文\\实验——代码\\PlatEMO v1.5 (2017-12)\\Data/%s/%s.mat',problem,problem));
    
    %保存数据到excel
    if exist(filepath,'file')
        delete(filepath);
    end
    save2excel(PF.PF.objs,title,filepath);
        
    %绘制气泡图
    Draw_bubble(PF.PF.objs);
end

function [  ] = Draw_bubble( data )
%% 绘制气泡图
    x = data(:,2)';
    y = data(:,3)';
    z = data(:,4)';
    
    %生成颜色
 	gnums = unique(data(:,1));
	c = parula(length(gnums));    
    c = fliplr(c')';
    
    %根据使用的罐的个数绘制图像
    figure;
    g_num = data(:,1);
    legendtext = cell(1,length(gnums));
    for i=1:length(gnums)
        ind = find(g_num== gnums(i));  %供油罐个数为tmp(i)的记录下标
        scatter3(x(ind),y(ind),z(ind),80,c(i,:),'fill');%气泡大小是切换次数
        legendtext{i} = num2str(gnums(i));
        hold on;
    end
    legend(legendtext);
    title('四目标Pareto最优解集');
    range = [min(x), max(x), min(y), max(y), min(z), max(z)];
    axis(range);% 坐标轴的显示范围
    xlabel('管道混合次数');% x轴名称、字体及大小
    ylabel('罐底混合次数');% y轴名称、字体及大小
    zlabel('供油罐切换次数');% z轴名称、字体及大小
end