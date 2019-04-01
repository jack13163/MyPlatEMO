function [ ] = RecordResult( Global )
    persistent global_result;
    
    % 设置静态变量为空
    if Global.evaluated <= Global.N
        global_result = [];%第一次运行时，设置上次函数运行留下的persistent变量为空
    end
    
    %% 判断目录是否存在
	method = func2str(Global.algorithm);
	question = func2str(Global.problem);
    dir_1 = sprintf('Data/%s/',question);                               %问题对应的第一级目录
    dir_2 = sprintf('%s%d/',dir_1,Global.run);                       %问题对应的第二级目录，单次运行结果存放在一个文件夹里面，文件夹的名称为运行序号
    %判断路径是否存在    
    if ~exist(dir_1,'dir')
        mkdir(dir_1);
    end
    if ~exist(dir_2,'dir')
        mkdir(dir_2);
    end
    
    %% 输出运行过程中的信息
    clc;%清空命令行
    fprintf('run %d\n',Global.run);
    fprintf('%s with %s on %s\n',method,func2str(Global.operator),question);
    fprintf('N: %d     evaluation: %d\n',Global.N,Global.evaluation);
    fprintf('status: (%6.2f%%), %.2fs passed...\n',Global.evaluated/Global.evaluation*100,Global.runtime);
    
    [r,~] = size(Global.result);
   	result = {Global.result{r,:}};                                           %每个算法最后一次运行结果
    global_result = [global_result;result];
    
    %% 运行完毕后，保存实验结果到本地
	if Global.evaluated >= Global.evaluation
        filename = sprintf('%s%s.mat',dir_1,question);              %存放汇总的运行结果（近似PF）
        filename2 = sprintf('%s%s.mat',dir_2,method);              %存放每个算法的所有运行结果
        if exist(filename2,'file')
            delete(filename2);
        end
        save(filename2,'global_result');
        
        if exist(filename,'file')
            PF = load(filename);
            PF = CombineArchive(result{1,2},PF.PF,500);       %设置储备集中元素的个数
        else
            PF = CombineArchive(result{1,2},[],500);     %个体组成的PF
        end
       	fprintf('\nnumber of Non-dominant Solutions :%d\n',size(PF.objs,1));
       	save(filename,'PF');        %保存结果
 	end
end

%%   合并储备集
function [ Archive ] = CombineArchive( Population,OldArchive,Na )
    %% 若种群为空，则直接返回原储备集
    if isempty(Population)
        Archive = OldArchive;
        return;
    end
    
    %% 原来储备集为空
    if isempty(OldArchive)
      	%找到第一个可行的解，作为储备集的初始解，并更新种群为剩余解
       	OldArchive = Population(1);
      	Population = Population(2:end);
    end
    
	M = size(OldArchive.objs,2);%目标个数
    
    %% 原储备集非空
    for i = 1:length(Population)
        %将每个候选者与外部储备集中的元素进行比较，更新储备集
        OldArchive = up_vac0(Population(i),OldArchive,M,Na);
    end
	Archive = OldArchive;
end

%% 更新外部储备集
function [ new_AC ] = up_vac0( par_eff,old_AC,M,Na )
    %% 更新储备集
    ss1 = length(old_AC);
    old_AC_Objs = old_AC.objs;
    par_eff_Obj = par_eff.objs;
    %判断粒子与外部储备集中粒子的支配关系
    for i = 1:ss1
        bb1 = 0;
        bb2 = 0;
        for j = 1:M     %目标函数的个数M
            aa1 = old_AC_Objs(i,j);
            aa2 = par_eff_Obj(1,j);
            if aa2 < aa1
                bb1 = bb1 + 1;
            elseif aa2 == aa1
                bb2 = bb2 + 1;
            end
        end
        %判断支配关系
        if bb1 == M         %候选者支配old_AC(i,:)，即候选者支配外部储备集中的当前粒子
            old_AC_Objs(i,:) = inf;          %删除外部储备集中被候选者支配的粒子
        elseif bb2 > 0 && bb1 == M - bb2       %候选者支配old_AC(i,:)，即候选者支配外部储备集中的当前粒子
            old_AC_Objs(i,:) = inf;
        elseif bb1 + bb2 == 0       %候选者不支配old_AC(i,:)
            par_eff_Obj(1,:) = inf;         %删除候选者粒子，候选者粒子被淘汰，不能进入外部储备集
            break;      %当前候选者支配于外部储备集中的当前粒子，退出
        elseif bb2 ~= 0 && bb1 == 0         %候选者不支配old_AC(i,:)
            par_eff_Obj(1,:) = inf;
            break;
        end
    end
    %将候选者加入到临时外部储备集中，临时外部储备集中被支配的粒子所在的行变为inf，即无穷大
    tmp_AC = [par_eff,old_AC]; 
    tmp_AC_Objs = [par_eff_Obj;old_AC_Objs]; 
    new_AC = [];
    %遍历合并后的临时外部储备集
    for i = 1:ss1+1
        %判断当前粒子是否被其他粒子所支配，若否，则将其纳入新的外部储备集
        if tmp_AC_Objs(i,1) ~= inf
            new_AC = [new_AC,tmp_AC(i)];
        end
    end
    
    %% 删除多余节点
    new_AC_Objs = new_AC.objs;
    ss2 = length(new_AC);
    %若外部储备集中元素个数超出其最大容量
    if ss2 > Na
        %遍历各个目标函数
        lim_f = zeros(5,2);
        for i = 1:M
            lim_f(i,2) = max(new_AC_Objs(:,i));       %求目标函数i的极大值
            lim_f(i,1) = min(new_AC_Objs(:,i));        %求目标函数i的极小值
        end
        [m,n] = size(new_AC_Objs);
        DD = zeros(m,n);        %拥挤度
        deep = zeros(m,1);
        for i = 1:M
           [~, ind] = sort(new_AC_Objs(:,i));        %将目标函数i的函数值按照从小到大的顺序排列
           %外部储备集中粒子个数ss2>50，又当前只有一个粒子进入外部储备集，故当前情况下，外部储备集中粒子个数为NA+1
           for j = 1:m
               %两端的拥挤度最大
               if lim_f(i,1) == new_AC_Objs(ind(j),i) || lim_f(i,2) == new_AC_Objs(ind(j),i)
                   DD(ind(j),i) = inf;
               else
                   DD(ind(j),i) = new_AC_Objs(ind(j+1),i) - new_AC_Objs(ind(j-1),i)/(lim_f(i,2) - lim_f(i,1));
               end
           end
        end
        %从外部储备集中删除拥挤距离值最小的粒子
        for jj = 1:m
            deep(jj) = sum(DD(jj,:));
        end
        [~, ind] = sort(deep);
        new_AC(ind(1:end-Na)) = [];
    end
end