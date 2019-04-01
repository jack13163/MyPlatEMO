function [  ] = show_special_results( )
%显示3个特殊的解
    inds = [18 10 41];
    PF = load('F:\多目标优化论文\实验——代码\PlatEMO v1.5 (2017-12)\Data\DMOPs4_OS\DMOPs4_OS.mat');
    lim = [0,40;0,50;0,20;0,10;0,200;];
    labels = {'Pipeline mixing cost','Tank bottom mixing cost','Number of switchings','Number of tanks','energe cost'};
    for i=1:length(inds)
        figure;
        data = PF.PF(inds(i)).obj;
        Draw_radar(data,lim,labels);
        figure;
        [~,plan]=fat_4(PF.PF(inds(i)).dec);
        printgante(plan{1},sprintf('调度%d',i));
    end
end

