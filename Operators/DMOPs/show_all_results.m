%% ��ʾ������PF���
function [] = show_all_results(  )
    problem = 'DMOPs4_OS';                  %����
    title = {'�ܵ���ϳɱ�','�޵׻�ϳɱ�','���͹��л�����','���͹޸���','�ܺĳɱ�'};
    filepath = sprintf('Data/%s/%s.xls',problem,problem);
    PF=load(sprintf('%s.mat',problem));
    
    %�������ݵ�excel
    if exist(filepath,'file')
        delete(filepath);
    end
    save2excel(PF.PF.objs,title,filepath);
        
    %��������ͼ
    Draw_bubble(PF.PF.objs);
end