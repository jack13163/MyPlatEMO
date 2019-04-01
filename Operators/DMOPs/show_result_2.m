%% ��ʾʵ����
function show_result_2()    
    problem = 'DMOPs2_OS';                  %����
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
        'GrEA','DE'};                                  %�Ա��㷨�����������(��һ���㷨������������ܵ��㷨)
    metrics = 'IGD';                                 %�Աȵ�ָ��
    
    % ����ָ��ֵ
    for k=1:length(runs)
        figure;
        for i=1:length(algorithmsandoperators)/2
            subplot(3,3,i);
            result_compare( problem,algorithmsandoperators{2*i-1},metrics,runs(k) );
        end
    end
    
    % �������runsʵ������ָ��ľ�ֵ�ͷ���
    figure;
    result_compare2(problem,{algorithmsandoperators{1:2:end}},metrics,runs);
end

%% ���Ʋ�ͬ�㷨������еľ�ֵ�������ܱȽ�ͼ
function [ result_mean, result_std ] = result_compare2( problem,algorithms,metrics,runs )
    % �����ļ��Ƿ����
    base_path = 'Data';
    method = str2func(metrics);
    
    % ����pf����
    pf_filepath = sprintf('%s/%s/%s.mat',base_path,problem,problem);
    if ~exist(pf_filepath,'file')
        disp('��ȷ��ʵ�������ļ����ڣ�\n');
        return;
    end
    pf = load(pf_filepath);
	%��׼��
  	PF_Objs = pf.PF.objs;
   	Objs_Max = max(PF_Objs);
   	Objs_Min = min(PF_Objs);
	m = size(PF_Objs,1);
	New_PF_Objs =  (PF_Objs - repmat(Objs_Min,m,1))./(repmat(Objs_Max - Objs_Min,m,1));
    
    % ����ʵ������
	result = zeros(length(runs), length(algorithms));        %result��length(runs)�У�length(algorithms)��
    for i=1:length(runs)
        for j=1:length(algorithms)
            algo_filepath = sprintf('%s/%s/%d/%s.mat',base_path,problem,runs(i),algorithms{j});%dualΪ���Ƚϵ�ʵ�����ļ���
            if  ~exist(algo_filepath,'file')
                disp('��ȷ��ʵ�������ļ����ڣ�\n');
                return;
            end
            algo = load(algo_filepath);

            %����IGDֵ
            [t,~] = size(algo.global_result);
            A_Objs = algo.global_result{t,2}.objs;      %ȡ���һ��ʵ����
            n = size(A_Objs,1);
            New_A_Objs = (A_Objs - repmat(Objs_Min,n,1))./(repmat(Objs_Max - Objs_Min,n,1));
            result(i,j) = method(New_A_Objs,New_PF_Objs);
        end
    end
    

    if length(runs) > 1
        %�Ա������㷨�����ܣ����е�һ��Ϊ�������㷨������Ϊ�ο��Ա��㷨
        h = zeros(1,size(result,2)-1);
        p = zeros(1,size(result,2)-1);
        for i=1:length(h)
            [h(i),p(i)] = ttest2(result(:,1),result(:,i+1),'Alpha',0.05,'tail','both');    %˫�߼��飬���Ŷ�Ϊ0.05
            if h(i) > 0
                if mean(result(:,1)) < mean(result(:,i+1))
                    fprintf('��95%%�İ�����Ϊ%s����%s�����ܣ�����Ϊ��%d\n', algorithms{1}, algorithms{i+1}, p(i));
                else
                    fprintf('��95%%�İ�����Ϊ%s����%s�����ܣ�����Ϊ��%d\n', algorithms{1}, algorithms{i+1}, p(i));
                end
            end
        end
        %����ͼ��
        result_mean = mean(result);
        result_std = std(result);
        shape = {'x','+','*','o','s','d','v','<','>','p','h','.'};      %�����12����״
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
        %����ͼ��
        x = 1:length(algorithms);
        plot(x,result);
        title(sprintf('%s on %s runs %d times',strrep(problem,'_','\_'),metrics,length(runs)));
        xlabel('algorithms');
        ylabel(metrics);
        axis([1,8,0,1]) % ����x���y��Ļ��Ʒ�Χ
        set(gca,'xticklabel',algorithms);
        grid on;
    end
end

%% ���Ʋ�ͬ�㷨��k�����е�ָ��Ա�ͼ
function [ result ] = result_compare( problem,algorithm,metrics,k )
    % �����ļ��Ƿ����
    base_path = 'Data';
    method = str2func(metrics);
    pf_filepath = sprintf('%s/%s/%s.mat',base_path,problem,problem);
 	algo_filepath = sprintf('%s/%s/%d/%s.mat',base_path,problem,k,algorithm);%kΪʵ�����
    if ~exist(pf_filepath,'file') || ~exist(algo_filepath,'file')
        disp('��ȷ��ʵ�������ļ����ڣ�\n');
        return;
    end
    %����ʵ������
    pf = load(pf_filepath);
   	algo = load(algo_filepath);
    %��׼��
    PF_Objs = pf.PF.objs;
    Objs_Max = max(PF_Objs);
    Objs_Min = min(PF_Objs);
    m = size(PF_Objs,1);
    New_PF_Objs =  (PF_Objs - repmat(Objs_Min,m,1))./(repmat(Objs_Max - Objs_Min,m,1));
    
    %����IGDֵ
 	[t,~] = size(algo.global_result);
    result = zeros(1,t);        %t�������¼����
    for i=1:t
        A_Objs = algo.global_result{i,2}.objs;
        m = size(A_Objs,1);
        New_A_Objs = (A_Objs - repmat(Objs_Min,m,1))./(repmat(Objs_Max - Objs_Min,m,1));
        result(i) = method(New_A_Objs,New_PF_Objs);
    end
    %����ͼ��
    x = [algo.global_result{:,1}];
    plot(x,result);
    title(sprintf('%s on %s',algorithm,strrep(problem,'_','\_')));
    xlabel('evaluation');
    ylabel(metrics);
    axis([1 max(x) 0 1]);
    %set(gca,'XTick',1:t);
    grid on;
end