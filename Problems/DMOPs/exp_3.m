%% ʵ��3������ʦ��
function [] = exp_3(runs)
    %% ʵ���������
    N = 100;                                           %��Ⱥ��ģ
    evaluation = 100*N;                          %���۴���
    problem = 'DMOPs3_OS';                  %����
    algorithmsandoperators = {
        'NSGAII','DE',...
        'NSGAIII','DE',...
        'IBEA','EAreal',...
        'ARMOEA','FEP',...
        'RVEA','DE',...
        'MOEAD','DE',...
        'KnEA','DE',...
        'GrEA','DE'};                                  %�Ա��㷨�����������(��һ���㷨������������ܵ��㷨)
    
    %% ����ʵ��
    do_exp(N,evaluation,problem,algorithmsandoperators,runs);
    show_save_pf(problem,false,false);                   %�������������Ƿ���ʾ��ϸ����ͼ(ֻ����ڵ�������)��pf��ͼ
end

%% ��ʾ������PF���
function [] = show_save_pf(problem,flag1,flag2)
    title = {'�ܵ���ϳɱ�','�޵׻�ϳɱ�','���͹��л�����','���͹޸���','�ܺĳɱ�'};
    filepath = sprintf('Data/%s/%s.xls',problem,problem);
    PF=load(sprintf('%s.mat',problem));
    
    %�������ݵ�excel
    save2excel(PF.PF.objs,title,filepath);
    
    %��ʾ��ϸ���ȸ���ͼ
    if flag1
        [~, scheduleplans] = fat_3(PF.PF.decs);
        for i = 1:length(scheduleplans)
            figure;
            printgante(scheduleplans{i});
            pause(1);%��ͣ1�������
        end
    end
    
    %��ʾ��ͼ
    if flag2
        boxplot(PF.PF.objs);
    end
end