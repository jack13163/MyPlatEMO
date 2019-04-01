function [ Gbest ] = UpdateGbest( Archive, N )
%UPDATEGBEST 
    %% Find the non-dominated solutions
    Gbest_crowd_val = CrowdingDistance(Archive);
    num = size(Archive,2);
    %�������е����ӣ������ӵ�ȫ��������
    for i = 1:N                                                    %���Ӹ���Ϊpopsize
        a1 = randi(num,1);                           %���ⲿ�����������ѡ����������
        a2 = randi(num,1);
        %�Ƚ��������ӵ�ӵ������ֵ��ѡ��ӵ������ֵ����ⲿ�������е�������Ϊ�����ӵ�ȫ��������
        if Gbest_crowd_val(a1) >= Gbest_crowd_val(a2)
            Gbest(i) = Archive(a1);
        else
            Gbest(i) = Archive(a2);
        end
    end
end

