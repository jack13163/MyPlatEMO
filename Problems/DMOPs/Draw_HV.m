%% 绘制趋势图function [] = Draw_HV()    runs = 7;    problem = 'DMOPs1_OS';                            %问题      % base_path = 'F:\多目标优化论文\实验——代码\PlatEMO v1.5 (2017-12)\Data\DMOPs1_OS\';    algorithms = {        'MPHH1',...        'BBMOPSO',...        'MOPSO',...        'NSGAII',...        'NSGAIII'};                           %对比算法及其操作算子(第一个算法代表待检验性能的算法)    % 绘制HV收敛曲线    result_compare3( problem,algorithms,runs );end%% 计算不同算法的c指标function [] = result_compare3( problem,algorithms,runs )    % 查找文件是否存在    base_path = 'F:\多目标优化论文\实验——代码\PlatEMO v1.5 (2017-12)\Data';    shape = {'+','d','x','^','v','s','<','>'};    m = 4;  %目标个数        %加载实验数据    tmp_result = cell(1,length(algorithms));    for i=1:length(algorithms)        filepath =sprintf('%s/%s/%d/%s.mat',base_path,problem,runs,algorithms{i});            if ~exist(filepath,'file')            disp('请确定实验数据文件存在！\n');            return;        end        tmp_result{i} = load(filepath);    end            %% 寻找最大值和最小值    Objs_Max = zeros(1,m);     % 找到最大的值    Objs_Min = ones(1,m) * inf;    for i=1:length(algorithms)        for j=1:size(tmp_result{i}.global_result,1)            Objs_Max = max([tmp_result{i}.global_result{j,2}.objs;Objs_Max]);            Objs_Min = min([tmp_result{i}.global_result{j,2}.objs;Objs_Min]);        end    end    Objs_Max = Objs_Max + 1;    Objs_Min = Objs_Min - 1;                    %% 计算HV指标    figure;    hold on;    title(sprintf('%s on HV',problem),'Interpreter','none');	xlabel('Evaluation');   	ylabel('HV');   	grid on;    xmax = 0;    for i=1:length(algorithms)        % 挑选100个数据        [t,~] = size(tmp_result{i}.global_result);%t为算法i的记录个数        sample_num = 100;%采样点的个数        if t > sample_num            ind = sort(randperm(t,sample_num));        else            ind = 1:t;        end        result_hv = zeros(1,length(ind));        %t：结果记录次数                try            for j=1:length(ind)                A_Objs = tmp_result{i}.global_result{ind(j),2}.objs;%取个体ind(j)                A_Objs = unique(A_Objs,'rows');                n = size(A_Objs,1);                A_Objs = (A_Objs - repmat(Objs_Min,n,1))./(repmat(Objs_Max - Objs_Min,n,1));                result_hv(j) = HV(A_Objs,ones(1,m));%记录第j个个体的hv值                fprintf('%f%%\n',round(j/length(ind)*100,2));            end        catch ex            fprintf(ex.message);        end                %绘制图像        x = [tmp_result{i}.global_result{:,1}];        h1 = plot(x(ind),result_hv);        set(h1,'marker',shape{i});        set(h1,'markersize',4);        xmax = max(max(x),xmax);    end   	axis([1 max(x) 0 1]);    legend(algorithms);end