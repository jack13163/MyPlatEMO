function [MaxObjV,MaxChrom]=EliteInduvidual(Chrom,ObjV,MaxObjV,MaxChrom)
%% �˹�ѡ������
MP=length(Chrom);  %��Ⱥ��
for i=1:MP
    [MaxO,maxI]=max(ObjV{i});   %�ҳ���i��Ⱥ�����Ÿ���
    if MaxO>MaxObjV(i)
        MaxObjV(i)=MaxO;         %��¼����Ⱥ�ľ�������
        MaxChrom(i,:)=Chrom{i}(maxI,:);  %��¼����Ⱥ��������ı���
    end
end
