%% 显示并保存PF结果
function [] = show_all_results(  )
    problem = 'DMOPs4_OS';                  %问题
    title = {'管道混合成本','罐底混合成本','供油罐切换次数','供油罐个数','能耗成本'};
    filepath = sprintf('Data/%s/%s.xls',problem,problem);
    PF=load(sprintf('%s.mat',problem));
    
    %保存数据到excel
    if exist(filepath,'file')
        delete(filepath);
    end
    save2excel(PF.PF.objs,title,filepath);
        
    %绘制气泡图
    Draw_bubble(PF.PF.objs);
end