%% 显示实验结果
function show_result_2()    
    problem = 'DMOPs2_OS';                  %问题
    runs = 1:10;
    algorithmsandoperators = {
        'MPHH2','DE',...
        'NSGAII','DE',...
        'NSGAIII','DE',...
        'IBEA','EAreal',...
        'ARMOEA','FEP',...
        'RVEA','DE',...
        'MOEAD','DE',...
        'KnEA','DE',...
        'GrEA','DE'};                                  %对比算法及其操作算子(第一个算法代表待检验性能的算法)
    metrics = 'IGD';                                 %对比的指标
    
    % 计算指标值
    for k=1:length(runs)
        figure;
        for i=1:length(algorithmsandoperators)/2
            subplot(3,3,i);
            result_compare( problem,algorithmsandoperators{2*i-1},metrics,runs(k) );
        end
    end
    
    % 计算序号runs实验结果的指标的均值和方差
    figure;
    result_compare2(problem,{algorithmsandoperators{1:2:end}},metrics,runs);
end

%% 绘制不同算法多次运行的均值方差性能比较图
function [ result_mean, result_std ] = result_compare2( problem,algorithms,metrics,runs )
    % 查找文件是否存在
    base_path = 'Data';
    method = str2func(metrics);
    
    % 加载pf数据
    pf_filepath = sprintf('%s/%s/%s.mat',base_path,problem,problem);
    if ~exist(pf_filepath,'file')
        disp('请确定实验数据文件存在！\n');
        return;
    end
    pf = load(pf_filepath);
	%标准化
  	PF_Objs = pf.PF.objs;
   	Objs_Max = max(PF_Objs);
   	Objs_Min = min(PF_Objs);
	m = size(PF_Objs,1);
	New_PF_Objs =  (PF_Objs - repmat(Objs_Min,m,1))./(repmat(Objs_Max - Objs_Min,m,1));
    
    % 加载实验数据
	result = zeros(length(runs), length(algorithms));        %result：length(runs)行，length(algorithms)列
    for i=1:length(runs)
        for j=1:length(algorithms)
            algo_filepath = sprintf('%s/%s/%d/%s.mat',base_path,problem,runs(i),algorithms{j});%dual为待比较的实验结果的集合
            if  ~exist(algo_filepath,'file')
                disp('请确定实验数据文件存在！\n');
                return;
            end
            algo = load(algo_filepath);

            %计算IGD值
            [t,~] = size(algo.global_result);
            A_Objs = algo.global_result{t,2}.objs;      %取最后一次实验结果
            n = size(A_Objs,1);
            New_A_Objs = (A_Objs - repmat(Objs_Min,n,1))./(repmat(Objs_Max - Objs_Min,n,1));
            result(i,j) = method(New_A_Objs,New_PF_Objs);
        end
    end
    

    if length(runs) > 1
        %对比两种算法的性能，其中第一个为待检验算法，其他为参考对比算法
        h = zeros(1,size(result,2)-1);
        p = zeros(1,size(result,2)-1);
        for i=1:length(h)
            [h(i),p(i)] = ttest2(result(:,1),result(:,i+1),'Alpha',0.05,'tail','both');    %双边检验，置信度为0.05
            if h(i) > 0
                if mean(result(:,1)) < mean(result(:,i+1))
                    fprintf('有95%%的把握认为%s优于%s的性能，概率为：%d\n', algorithms{1}, algorithms{i+1}, p(i));
                else
                    fprintf('有95%%的把握认为%s劣于%s的性能，概率为：%d\n', algorithms{1}, algorithms{i+1}, p(i));
                end
            end
        end
        %绘制图像
        result_mean = mean(result);
        result_std = std(result);
        shape = {'x','+','*','o','s','d','v','<','>','p','h','.'};      %最多有12种形状
        title(sprintf('%s on %s runs %d times',strrep(problem,'_','\_'),metrics,length(runs)));
        xlabel('mean');
        ylabel('std');
        hold on
        for i=1:length(result_mean)
            x = result_mean(i);
            y = result_std(i);
            plot(x,y,shape{i});
            text(x+0.004, y, algorithms{i}, 'fontsize', 10);
        end
        grid on;
   	else            
        %绘制图像
        x = 1:length(algorithms);
        plot(x,result);
        title(sprintf('%s on %s runs %d times',strrep(problem,'_','\_'),metrics,length(runs)));
        xlabel('algorithms');
        ylabel(metrics);
        axis([1,8,0,1]) % 设置x轴和y轴的绘制范围
        set(gca,'xticklabel',algorithms);
        grid on;
    end
end

%% 绘制不同算法第k次运行的指标对比图
function [ result ] = result_compare( problem,algorithm,metrics,k )
    % 查找文件是否存在
    base_path = 'Data';
    method = str2func(metrics);
    pf_filepath = sprintf('%s/%s/%s.mat',base_path,problem,problem);
 	algo_filepath = sprintf('%s/%s/%d/%s.mat',base_path,problem,k,algorithm);%k为实验次数
    if ~exist(pf_filepath,'file') || ~exist(algo_filepath,'file')
        disp('请确定实验数据文件存在！\n');
        return;
    end
    %加载实验数据
    pf = load(pf_filepath);
   	algo = load(algo_filepath);
    %标准化
    PF_Objs = pf.PF.objs;
    Objs_Max = max(PF_Objs);
    Objs_Min = min(PF_Objs);
    m = size(PF_Objs,1);
    New_PF_Objs =  (PF_Objs - repmat(Objs_Min,m,1))./(repmat(Objs_Max - Objs_Min,m,1));
    
    %计算IGD值
 	[t,~] = size(algo.global_result);
    result = zeros(1,t);        %t：结果记录次数
    for i=1:t
        A_Objs = algo.global_result{i,2}.objs;
        m = size(A_Objs,1);
        New_A_Objs = (A_Objs - repmat(Objs_Min,m,1))./(repmat(Objs_Max - Objs_Min,m,1));
        result(i) = method(New_A_Objs,New_PF_Objs);
    end
    %绘制图像
    x = [algo.global_result{:,1}];
    plot(x,result);
    title(sprintf('%s on %s',algorithm,strrep(problem,'_','\_')));
    xlabel('evaluation');
    ylabel(metrics);
    axis([1 max(x) 0 1]);
    %set(gca,'XTick',1:t);
    grid on;
end