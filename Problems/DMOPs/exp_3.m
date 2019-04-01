%% 实验3（侯老师）
function [] = exp_3(runs)
    %% 实验参数设置
    N = 100;                                           %种群规模
    evaluation = 100*N;                          %评价次数
    problem = 'DMOPs3_OS';                  %问题
    algorithmsandoperators = {
        'NSGAII','DE',...
        'NSGAIII','DE',...
        'IBEA','EAreal',...
        'ARMOEA','FEP',...
        'RVEA','DE',...
        'MOEAD','DE',...
        'KnEA','DE',...
        'GrEA','DE'};                                  %对比算法及其操作算子(第一个算法代表待检验性能的算法)
    
    %% 进行实验
    do_exp(N,evaluation,problem,algorithmsandoperators,runs);
    show_save_pf(problem,false,false);                   %布尔变量控制是否显示详细调度图(只针对于调度问题)、pf盒图
end

%% 显示并保存PF结果
function [] = show_save_pf(problem,flag1,flag2)
    title = {'管道混合成本','罐底混合成本','供油罐切换次数','供油罐个数','能耗成本'};
    filepath = sprintf('Data/%s/%s.xls',problem,problem);
    PF=load(sprintf('%s.mat',problem));
    
    %保存数据到excel
    save2excel(PF.PF.objs,title,filepath);
    
    %显示详细调度甘特图
    if flag1
        [~, scheduleplans] = fat_3(PF.PF.decs);
        for i = 1:length(scheduleplans)
            figure;
            printgante(scheduleplans{i});
            pause(1);%暂停1秒后运行
        end
    end
    
    %显示盒图
    if flag2
        boxplot(PF.PF.objs);
    end
end