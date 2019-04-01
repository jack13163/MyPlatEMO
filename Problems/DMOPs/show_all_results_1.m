%% ��ʾ������PF���
function [] = show_all_results_1(  )
    problem = 'DMOPs1_OS';                  %����
    title = {'�ܵ���ϳɱ�','���͹��л��ɱ�','�޵׻�ϳɱ�','���͹޸���'};
    filepath = sprintf('F:\\��Ŀ���Ż�����\\ʵ�顪������\\PlatEMO v1.5 (2017-12)\\Data/%s/%s.xls',problem,problem);
    PF=load(sprintf('F:\\��Ŀ���Ż�����\\ʵ�顪������\\PlatEMO v1.5 (2017-12)\\Data/%s/%s.mat',problem,problem));
    
    %�������ݵ�excel
    if exist(filepath,'file')
        delete(filepath);
    end
    save2excel(PF.PF.objs,title,filepath);
        
    %��������ͼ
    Draw_bubble(PF.PF.objs);
end

function [  ] = Draw_bubble( data )
%% ��������ͼ
    x = data(:,2)';
    y = data(:,3)';
    z = data(:,4)';
    
    %������ɫ
 	gnums = unique(data(:,1));
	c = parula(length(gnums));    
    c = fliplr(c')';
    
    %����ʹ�õĹ޵ĸ�������ͼ��
    figure;
    g_num = data(:,1);
    legendtext = cell(1,length(gnums));
    for i=1:length(gnums)
        ind = find(g_num== gnums(i));  %���͹޸���Ϊtmp(i)�ļ�¼�±�
        scatter3(x(ind),y(ind),z(ind),80,c(i,:),'fill');%���ݴ�С���л�����
        legendtext{i} = num2str(gnums(i));
        hold on;
    end
    legend(legendtext);
    title('Pareto���Ž⼯');
    range = [min(x), max(x), min(y), max(y), min(z), max(z)];
    axis(range);% ���������ʾ��Χ
    xlabel('���͹��л�����');% x�����ơ����弰��С
    ylabel('�޵׻�ϳɱ�');% y�����ơ����弰��С
    zlabel('�ܵ���ϳɱ�');% z�����ơ����弰��С
end