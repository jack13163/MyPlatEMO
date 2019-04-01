function [ ] = RecordResult( Global )
    persistent global_result;
    
    % ���þ�̬����Ϊ��
    if Global.evaluated <= Global.N
        global_result = [];%��һ������ʱ�������ϴκ����������µ�persistent����Ϊ��
    end
    
    %% �ж�Ŀ¼�Ƿ����
	method = func2str(Global.algorithm);
	question = func2str(Global.problem);
    dir_1 = sprintf('Data/%s/',question);                               %�����Ӧ�ĵ�һ��Ŀ¼
    dir_2 = sprintf('%s%d/',dir_1,Global.run);                       %�����Ӧ�ĵڶ���Ŀ¼���������н�������һ���ļ������棬�ļ��е�����Ϊ�������
    %�ж�·���Ƿ����    
    if ~exist(dir_1,'dir')
        mkdir(dir_1);
    end
    if ~exist(dir_2,'dir')
        mkdir(dir_2);
    end
    
    %% ������й����е���Ϣ
    clc;%���������
    fprintf('run %d\n',Global.run);
    fprintf('%s with %s on %s\n',method,func2str(Global.operator),question);
    fprintf('N: %d     evaluation: %d\n',Global.N,Global.evaluation);
    fprintf('status: (%6.2f%%), %.2fs passed...\n',Global.evaluated/Global.evaluation*100,Global.runtime);
    
    [r,~] = size(Global.result);
   	result = {Global.result{r,:}};                                           %ÿ���㷨���һ�����н��
    global_result = [global_result;result];
    
    %% ������Ϻ󣬱���ʵ����������
	if Global.evaluated >= Global.evaluation
        filename = sprintf('%s%s.mat',dir_1,question);              %��Ż��ܵ����н��������PF��
        filename2 = sprintf('%s%s.mat',dir_2,method);              %���ÿ���㷨���������н��
        if exist(filename2,'file')
            delete(filename2);
        end
        save(filename2,'global_result');
        
        if exist(filename,'file')
            PF = load(filename);
            PF = CombineArchive(result{1,2},PF.PF,500);       %���ô�������Ԫ�صĸ���
        else
            PF = CombineArchive(result{1,2},[],500);     %������ɵ�PF
        end
       	fprintf('\nnumber of Non-dominant Solutions :%d\n',size(PF.objs,1));
       	save(filename,'PF');        %������
 	end
end

%%   �ϲ�������
function [ Archive ] = CombineArchive( Population,OldArchive,Na )
    %% ����ȺΪ�գ���ֱ�ӷ���ԭ������
    if isempty(Population)
        Archive = OldArchive;
        return;
    end
    
    %% ԭ��������Ϊ��
    if isempty(OldArchive)
      	%�ҵ���һ�����еĽ⣬��Ϊ�������ĳ�ʼ�⣬��������ȺΪʣ���
       	OldArchive = Population(1);
      	Population = Population(2:end);
    end
    
	M = size(OldArchive.objs,2);%Ŀ�����
    
    %% ԭ�������ǿ�
    for i = 1:length(Population)
        %��ÿ����ѡ�����ⲿ�������е�Ԫ�ؽ��бȽϣ����´�����
        OldArchive = up_vac0(Population(i),OldArchive,M,Na);
    end
	Archive = OldArchive;
end

%% �����ⲿ������
function [ new_AC ] = up_vac0( par_eff,old_AC,M,Na )
    %% ���´�����
    ss1 = length(old_AC);
    old_AC_Objs = old_AC.objs;
    par_eff_Obj = par_eff.objs;
    %�ж��������ⲿ�����������ӵ�֧���ϵ
    for i = 1:ss1
        bb1 = 0;
        bb2 = 0;
        for j = 1:M     %Ŀ�꺯���ĸ���M
            aa1 = old_AC_Objs(i,j);
            aa2 = par_eff_Obj(1,j);
            if aa2 < aa1
                bb1 = bb1 + 1;
            elseif aa2 == aa1
                bb2 = bb2 + 1;
            end
        end
        %�ж�֧���ϵ
        if bb1 == M         %��ѡ��֧��old_AC(i,:)������ѡ��֧���ⲿ�������еĵ�ǰ����
            old_AC_Objs(i,:) = inf;          %ɾ���ⲿ�������б���ѡ��֧�������
        elseif bb2 > 0 && bb1 == M - bb2       %��ѡ��֧��old_AC(i,:)������ѡ��֧���ⲿ�������еĵ�ǰ����
            old_AC_Objs(i,:) = inf;
        elseif bb1 + bb2 == 0       %��ѡ�߲�֧��old_AC(i,:)
            par_eff_Obj(1,:) = inf;         %ɾ����ѡ�����ӣ���ѡ�����ӱ���̭�����ܽ����ⲿ������
            break;      %��ǰ��ѡ��֧�����ⲿ�������еĵ�ǰ���ӣ��˳�
        elseif bb2 ~= 0 && bb1 == 0         %��ѡ�߲�֧��old_AC(i,:)
            par_eff_Obj(1,:) = inf;
            break;
        end
    end
    %����ѡ�߼��뵽��ʱ�ⲿ�������У���ʱ�ⲿ�������б�֧����������ڵ��б�Ϊinf���������
    tmp_AC = [par_eff,old_AC]; 
    tmp_AC_Objs = [par_eff_Obj;old_AC_Objs]; 
    new_AC = [];
    %�����ϲ������ʱ�ⲿ������
    for i = 1:ss1+1
        %�жϵ�ǰ�����Ƿ�����������֧�䣬�������������µ��ⲿ������
        if tmp_AC_Objs(i,1) ~= inf
            new_AC = [new_AC,tmp_AC(i)];
        end
    end
    
    %% ɾ������ڵ�
    new_AC_Objs = new_AC.objs;
    ss2 = length(new_AC);
    %���ⲿ��������Ԫ�ظ����������������
    if ss2 > Na
        %��������Ŀ�꺯��
        lim_f = zeros(5,2);
        for i = 1:M
            lim_f(i,2) = max(new_AC_Objs(:,i));       %��Ŀ�꺯��i�ļ���ֵ
            lim_f(i,1) = min(new_AC_Objs(:,i));        %��Ŀ�꺯��i�ļ�Сֵ
        end
        [m,n] = size(new_AC_Objs);
        DD = zeros(m,n);        %ӵ����
        deep = zeros(m,1);
        for i = 1:M
           [~, ind] = sort(new_AC_Objs(:,i));        %��Ŀ�꺯��i�ĺ���ֵ���մ�С�����˳������
           %�ⲿ�����������Ӹ���ss2>50���ֵ�ǰֻ��һ�����ӽ����ⲿ���������ʵ�ǰ����£��ⲿ�����������Ӹ���ΪNA+1
           for j = 1:m
               %���˵�ӵ�������
               if lim_f(i,1) == new_AC_Objs(ind(j),i) || lim_f(i,2) == new_AC_Objs(ind(j),i)
                   DD(ind(j),i) = inf;
               else
                   DD(ind(j),i) = new_AC_Objs(ind(j+1),i) - new_AC_Objs(ind(j-1),i)/(lim_f(i,2) - lim_f(i,1));
               end
           end
        end
        %���ⲿ��������ɾ��ӵ������ֵ��С������
        for jj = 1:m
            deep(jj) = sum(DD(jj,:));
        end
        [~, ind] = sort(deep);
        new_AC(ind(1:end-Na)) = [];
    end
end