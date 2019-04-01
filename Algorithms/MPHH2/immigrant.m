function [Population]=immigrant(Population,MP,Num,mode)
%% 移民算子
move = zeros(MP,Num);            %待移民个体
delete = zeros(MP,Num);           %待替换个体
POPNUM = fix(length(Population)/MP);

if mode ==1
    %% 模式1:好的相互交换（环形）
    for i=1:MP
        %种群个体
        bird = ((i-1)*POPNUM+1):i*POPNUM;
        [FrontNo] = NDSort(Population(bird).objs,Population(bird).cons,POPNUM);     %非支配排序
        [~,ind] = sort(FrontNo);
        move(i,:) = (i-1) * POPNUM + ind(1:Num);       %选择优秀个体移民
    end
    t = MP;
    while t>1
        Population(move(t,:)) = Population(move(t-1,:));
        t = t - 1;
    end
    Population(move(1,:)) = Population(move(MP,:));
elseif mode == 2
    %% 模式2:随机交换（环形）
    for i=1:MP
        move(i,:) = (i-1) * POPNUM + randperm(POPNUM,Num);       %选择优秀个体移民
    end
    t = MP;
    while t>1
        Population(move(t,:)) = Population(move(t-1,:));
        t = t - 1;
    end
    Population(move(1,:)) = Population(move(MP,:));
elseif mode == 3    
    %% 模式3:好的替换坏的（环形）
    for i=1:MP
        %种群个体
        bird = ((i-1)*POPNUM+1):i*POPNUM;
        [FrontNo] = NDSort(Population(bird).objs,Population(bird).cons,POPNUM);     %非支配排序
        [~,ind] = sort(FrontNo);
        move(i,:) = (i-1) * POPNUM + ind(1:Num);       %选择优秀个体移民
        delete(i,:) = (i-1) * POPNUM + ind(POPNUM-Num+1:end);       %选择优秀个体移民
    end
    t = MP;
    while t>1
        Population(delete(t,:)) = Population(move(t-1,:));
        t = t - 1;
    end
    Population(delete(1,:)) = Population(move(MP,:));
elseif mode == 4
    %% 模式4:好的替换坏的（单向）
    for i=1:MP
        %种群个体
        bird = ((i-1)*POPNUM+1):i*POPNUM;
        [FrontNo] = NDSort(Population(bird).objs,Population(bird).cons,POPNUM);     %非支配排序
        [~,ind] = sort(FrontNo);
        move(i,:) = (i-1) * POPNUM + ind(1:Num);       %选择优秀个体移民
        delete(i,:) = (i-1) * POPNUM + ind(POPNUM-Num+1:end);       %选择优秀个体移民
    end
    t = MP;
    while t>1
        Population(delete(t,:)) = Population(move(t-1,:));
        t = t - 2;
    end
end
