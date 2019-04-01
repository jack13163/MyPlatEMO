%% ��ʾʵ����
function show_result_4()    
    runs = 1:300;
    turns = 10;                                       %ÿ10��ʵ�鹹��һ��
    problem = 'DMOPs4_OS';                  %����
    base_path = 'F:\��Ŀ���Ż�����\ʵ�顪������\PlatEMO v1.5 (2017-12)\Data\DMOPs4_OS\';
    algorithmsandoperators = {
        'NSGAIII','DE',...
        'IBEA','EAreal',...
        'NSGAII','DE',...
        'ARMOEA','FEP',...
        'RVEA','DE',...
        'MOEAD','DE',...
        'KnEA','DE',...
        'GrEA','DE'};                                  %�Ա��㷨�����������(��һ���㷨������������ܵ��㷨)
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
        % �������runsʵ������cָ��
        result = result + get_C(problem,{algorithmsandoperators{1:2:end}},group_runs);
    end
    result = result/max_turns;
    c_filename = sprintf('%sc.xls',base_path);
    if exist(c_filename,'file')
        delete(c_filename);
    end
    save2excel(result,[],c_filename);
end

%% ���㲻ͬ�㷨��cָ��
function [ result ] = get_C( problem,algorithms,runs )
    % �����ļ��Ƿ����
    base_path = 'F:\��Ŀ���Ż�����\ʵ�顪������\PlatEMO v1.5 (2017-12)\Data';
    result = zeros(length(runs),(length(algorithms)-1)*2);
    
    for i=1:length(runs)
        %����ʵ������
        tmp_result = cell(1,length(algorithms));
        for j=1:length(algorithms)
            filepath =sprintf('%s/%s/%d/%s.mat',base_path,problem,runs(i),algorithms{j});    
            if ~exist(filepath,'file')
                disp('��ȷ��ʵ�������ļ����ڣ�\n');
                return;
            end

            tmp_result{j} = load(filepath);
        end

        %% ����Cָ��
        for j=2:length(algorithms)
            [t1,~] = size(tmp_result{1}.global_result);               %�Լ����㷨
            [t2,~] = size(tmp_result{j}.global_result);               %�Ա��㷨
            A_Objs = tmp_result{1}.global_result{t1,2}.objs;      %ȡ���һ��ʵ����
            B_Objs = tmp_result{j}.global_result{t2,2}.objs;      %ȡ���һ��ʵ����
                
            result(i,2*(j-1)-1) = C(A_Objs,B_Objs);
            result(i,2*(j-1)) = C(B_Objs,A_Objs); 
        end
    end
end

%% ���㲻ͬ�㷨��cָ��
function [ result ] = result_compare3( problem,algorithms,runs )
    % �����ļ��Ƿ����
    base_path = 'Data';
    shape = {'+','d','x','^','v','s','<','>'};
    
    %����ʵ������
    tmp_result = cell(1,length(algorithms));
    for i=1:length(algorithms)
        filepath =sprintf('%s/%s/%d/%s.mat',base_path,problem,runs,algorithms{i});    
        if ~exist(filepath,'file')
            disp('��ȷ��ʵ�������ļ����ڣ�\n');
            return;
        end

        tmp_result{i} = load(filepath);
    end
    
    %����PF����
    pf_filepath = sprintf('%s/%s/%s.mat',base_path,problem,problem);
    if ~exist(pf_filepath,'file')
        disp('��ȷ��ʵ�������ļ����ڣ�\n');
        return;
    end
    pf = load(pf_filepath);
    
    %��׼��PF����
    PF_Objs = pf.PF.objs;
    Objs_Max = max(PF_Objs);
    Objs_Min = min(PF_Objs);
    m = size(PF_Objs,1);
    New_PF_Objs =  (PF_Objs - repmat(Objs_Min,m,1))./(repmat(Objs_Max - Objs_Min,m,1));
    
    %% ����IGDָ��
    figure;
    hold on;
    title(sprintf('%s on IGD',problem),'Interpreter','none');
	xlabel('Evaluation');
   	ylabel('IGD');
   	grid on;
    xmax = 0;
    for i=1:length(algorithms)
         %����IGDֵ
        [t,~] = size(tmp_result{i}.global_result);
        result_igd = zeros(1,t);        %t�������¼����
        for j=1:t
            A_Objs = tmp_result{i}.global_result{j,2}.objs;
            m = size(A_Objs,1);
            New_A_Objs = (A_Objs - repmat(Objs_Min,m,1))./(repmat(Objs_Max - Objs_Min,m,1));
            result_igd(j) = IGD(New_A_Objs,New_PF_Objs);
        end
        
        %����ͼ��
        ind = sort(randperm(t,100));
        x = [tmp_result{i}.global_result{:,1}];
        h1 = plot(x(ind),result_igd(ind));
        set(h1,'marker',shape{i});
        set(h1,'markersize',4);
        xmax = max(max(x),xmax);
    end
   	axis([1 max(x) 0 1]);
    legend(algorithms);
        
    %% ����HVָ��
    figure;
    hold on;
    title(sprintf('%s on HV',problem),'Interpreter','none');
	xlabel('Evaluation');
   	ylabel('HV');
   	grid on;
    xmax = 0;
    for i=1:length(algorithms)
         %����HVֵ������HV��ֵ������۽ϴ���ˣ���ȡ��������ķ�����
        [t,~] = size(tmp_result{i}.global_result);
        sample_num = 100;%������ĸ���
        if t > sample_num
            ind = sort(randperm(t,sample_num));
        else
            ind = 1:t;
        end
        result_hv = zeros(1,length(ind));        %t�������¼����
        for j=1:length(ind)
            A_Objs = tmp_result{i}.global_result{ind(j),2}.objs;%ȡ����ind(j)
            m = size(A_Objs,1);
            New_A_Objs = (A_Objs - repmat(Objs_Min,m,1))./(repmat(Objs_Max - Objs_Min,m,1));
            result_hv(j) = HV(New_A_Objs,New_PF_Objs);%��¼��j�������hvֵ
            fprintf('%f%%\n',round(j/length(ind)*100,2));
        end
        
        %����ͼ��
        x = [tmp_result{i}.global_result{:,1}];
        h1 = plot(x(ind),result_hv);
        set(h1,'marker',shape{i});
        set(h1,'markersize',4);
        xmax = max(max(x),xmax);
    end
   	axis([1 max(x) 0 1]);
    legend(algorithms);
    
    %% ����Cָ��
	result = cell(length(algorithms)+1);        %result��length(runs)�У�length(algorithms)��
    for i=1:length(algorithms)
        for j=1:length(algorithms)
            %���ñ���
            if i==1            
                result{i,j+1} = algorithms{j};
            end
            if j==1
                result{i+1,j} = algorithms{i};
            end
            if i~= j
                [t1,~] = size(tmp_result{i}.global_result);
                [t2,~] = size(tmp_result{j}.global_result);
                A_Objs = tmp_result{i}.global_result{t1,2}.objs;      %ȡ���һ��ʵ����
                B_Objs = tmp_result{j}.global_result{t2,2}.objs;      %ȡ���һ��ʵ����

                result{i+1,j+1} = C(A_Objs,B_Objs);
            else
                result{i+1,j+1} = 0;
            end
        end
    end
    xlswrite(sprintf('%s/%s/%d/cָ��.xls',base_path,problem,runs), result);                                      % ��resultд�뵽excel�ļ���
end