function [] = exp_4()
    runs = 1:30;                                   %10组实验，每组单独运行30次
    %% 实验参数设置
    N = [50 100 150 200 250 50 100 150 200 250];                                            %种群规模
    evaluation = [N(1:5) * 100, N(6:10) * 1000];                                 %评价次数
    problem = 'DMOPs4_OS';                              %问题
    algorithmsandoperators = {
        'KnEA','EAreal',...
        'NSGAII','EAreal',...
        'NSGAIII','EAreal',...
        'IBEA','EAreal',...
        'ARMOEA','FEP',...
        'RVEA','EAreal',...
        'MOEAD','EAreal',...
        'GrEA','EAreal'};                                  %对比算法及其操作算子(第一个算法代表待检验性能的算法)
    base_path = 'F:\多目标优化论文\实验——代码\PlatEMO v1.5 (2017-12)\Data\DMOPs4_OS\';
    
    %% 进行实验
    try
        for i=1:length(runs)
            for j=1:length(N)
                r = (runs(i)-1)*length(N) + j;
                dir_path = sprintf('%s%d',base_path,r);
                filenames = dir(sprintf('%s\\*.mat',dir_path));
                if ~exist(dir_path,'dir') || length(filenames) < length(algorithmsandoperators)/2
                    do_exp(N(j),evaluation(j),problem,algorithmsandoperators,r);
                end
            end
        end
    catch err
        fname = 'err_inf.txt';
    	fp = fopen(fname,'a');      %以追加的方式打开文件，若文件不存在，则创建
     	fprintf(fp,'%s\n',err.message);            
        exp_4();                        %重新执行程序
    end
end