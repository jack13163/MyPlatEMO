%% 显示实验结果
function show_result_4()    
    runs = 1:300;
    turns = 10;                                       %每10个实验构成一组
    problem = 'DMOPs4_OS';                  %问题
    base_path = 'F:\多目标优化论文\实验——代码\PlatEMO v1.5 (2017-12)\Data\DMOPs4_OS\';
    algorithmsandoperators = {
        'IBEA','EAreal',...
        'NSGAII','DE',...
        'NSGAIII','DE',...
        'ARMOEA','FEP',...
        'RVEA','DE',...
        'MOEAD','DE',...
        'KnEA','DE',...
        'GrEA','DE'};                                  %对比算法及其操作算子(第一个算法代表待检验性能的算法)
    dir_file = dir(base_path);
    k = 0;
    for i=1:length(dir_file)
        f_name = dir_file(i).name;
        if ~ismember('.',f_name)
            tmp = str2num(f_name);
            if tmp > max(runs)
                break;
            end
            if tmp > k
                k = tmp;
            end
        end
    end
    max_turns = floor(k/turns);
    result = zeros(turns,length(algorithmsandoperators)-2);
    for i=1:max_turns
        group_runs = (i-1)*turns+1:i*turns;
        % 计算序号runs实验结果的c指标
        result = result + get_C(problem,{algorithmsandoperators{1:2:end}},group_runs);
    end
    result = result/max_turns;
    c_filename = sprintf('%sc.xls',base_path);
    if exist(c_filename,'file')
        delete(c_filename);
    end
    save2excel(result,[],c_filename);
end

%% 计算不同算法的c指标
function [ result ] = get_C( problem,algorithms,runs )
    % 查找文件是否存在
    base_path = 'F:\多目标优化论文\实验——代码\PlatEMO v1.5 (2017-12)\Data';
    result = zeros(length(runs),(length(algorithms)-1)*2);
    
    for i=1:length(runs)
        %加载实验数据
        tmp_result = cell(1,length(algorithms));
        for j=1:length(algorithms)
            filepath =sprintf('%s/%s/%d/%s.mat',base_path,problem,runs(i),algorithms{j});    
            if ~exist(filepath,'file')
                disp('请确定实验数据文件存在！\n');
                return;
            end

            tmp_result{j} = load(filepath);
        end

        %% 计算C指标
        for j=2:length(algorithms)
            [t1,~] = size(tmp_result{1}.global_result);               %自己的算法
            [t2,~] = size(tmp_result{j}.global_result);               %对比算法
            A_Objs = tmp_result{1}.global_result{t1,2}.objs;      %取最后一次实验结果
            B_Objs = tmp_result{j}.global_result{t2,2}.objs;      %取最后一次实验结果
                
            result(i,2*(j-1)-1) = C(A_Objs,B_Objs);
            result(i,2*(j-1)) = C(B_Objs,A_Objs); 
        end
    end
end