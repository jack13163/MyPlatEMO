function [Population]=immigrant(Population,MP,Num,mode)
%% ��������
move = zeros(MP,Num);            %���������
delete = zeros(MP,Num);           %���滻����
POPNUM = fix(length(Population)/MP);

if mode ==1
    %% ģʽ1:�õ��໥���������Σ�
    for i=1:MP
        %��Ⱥ����
        bird = ((i-1)*POPNUM+1):i*POPNUM;
        [FrontNo] = NDSort(Population(bird).objs,Population(bird).cons,POPNUM);     %��֧������
        [~,ind] = sort(FrontNo);
        move(i,:) = (i-1) * POPNUM + ind(1:Num);       %ѡ�������������
    end
    t = MP;
    while t>1
        Population(move(t,:)) = Population(move(t-1,:));
        t = t - 1;
    end
    Population(move(1,:)) = Population(move(MP,:));
elseif mode == 2
    %% ģʽ2:������������Σ�
    for i=1:MP
        move(i,:) = (i-1) * POPNUM + randperm(POPNUM,Num);       %ѡ�������������
    end
    t = MP;
    while t>1
        Population(move(t,:)) = Population(move(t-1,:));
        t = t - 1;
    end
    Population(move(1,:)) = Population(move(MP,:));
elseif mode == 3    
    %% ģʽ3:�õ��滻���ģ����Σ�
    for i=1:MP
        %��Ⱥ����
        bird = ((i-1)*POPNUM+1):i*POPNUM;
        [FrontNo] = NDSort(Population(bird).objs,Population(bird).cons,POPNUM);     %��֧������
        [~,ind] = sort(FrontNo);
        move(i,:) = (i-1) * POPNUM + ind(1:Num);       %ѡ�������������
        delete(i,:) = (i-1) * POPNUM + ind(POPNUM-Num+1:end);       %ѡ�������������
    end
    t = MP;
    while t>1
        Population(delete(t,:)) = Population(move(t-1,:));
        t = t - 1;
    end
    Population(delete(1,:)) = Population(move(MP,:));
elseif mode == 4
    %% ģʽ4:�õ��滻���ģ�����
    for i=1:MP
        %��Ⱥ����
        bird = ((i-1)*POPNUM+1):i*POPNUM;
        [FrontNo] = NDSort(Population(bird).objs,Population(bird).cons,POPNUM);     %��֧������
        [~,ind] = sort(FrontNo);
        move(i,:) = (i-1) * POPNUM + ind(1:Num);       %ѡ�������������
        delete(i,:) = (i-1) * POPNUM + ind(POPNUM-Num+1:end);       %ѡ�������������
    end
    t = MP;
    while t>1
        Population(delete(t,:)) = Population(move(t-1,:));
        t = t - 2;
    end
end
