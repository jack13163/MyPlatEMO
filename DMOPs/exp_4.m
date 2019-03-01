function [] = exp_4()
    runs = 1:30;                                   %10��ʵ�飬ÿ�鵥������30��
    %% ʵ���������
    N = [50 100 150 200 250 50 100 150 200 250];                                            %��Ⱥ��ģ
    evaluation = [N(1:5) * 100, N(6:10) * 1000];                                 %���۴���
    problem = 'DMOPs4_OS';                              %����
    algorithmsandoperators = {
        'KnEA','EAreal',...
        'NSGAII','EAreal',...
        'NSGAIII','EAreal',...
        'IBEA','EAreal',...
        'ARMOEA','FEP',...
        'RVEA','EAreal',...
        'MOEAD','EAreal',...
        'GrEA','EAreal'};                                  %�Ա��㷨�����������(��һ���㷨������������ܵ��㷨)
    base_path = 'F:\��Ŀ���Ż�����\ʵ�顪������\PlatEMO v1.5 (2017-12)\Data\DMOPs4_OS\';
    
    %% ����ʵ��
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
    	fp = fopen(fname,'a');      %��׷�ӵķ�ʽ���ļ������ļ������ڣ��򴴽�
     	fprintf(fp,'%s\n',err.message);            
        exp_4();                        %����ִ�г���
    end
end