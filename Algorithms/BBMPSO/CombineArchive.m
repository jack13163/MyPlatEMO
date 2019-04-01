function [ Archive ] = CombineArchive( Population,OldArchive )
%   �ϲ�������
%   Detailed explanation goes here
    %% ����ȺΪ�գ���ֱ�ӷ���ԭ������
    if isempty(Population)
        Archive = OldArchive;
        return;
    end
        
    %% ԭ��������Ϊ��
    if isempty(OldArchive)
        PopulationObj = Population.objs;%Ŀ�����
        for i = 1:size(Population,2)
            %�ҵ���һ�����еĽ⣬��Ϊ�������ĳ�ʼ�⣬��������ȺΪʣ���
            if PopulationObj(i,1) ~= inf
                Archive = Population(i);
                Population = Population(i+1:end);
                break;
            end
        end
    else
        Archive = OldArchive;
    end
    
	M = size(Population.objs,2);%Ŀ�����
    
    %% ԭ�������ǿ�
    for i = 1:size(Population,2)
        %��ÿ����ѡ�����ⲿ�������е�Ԫ�ؽ��бȽϣ����´�����
        Archive = up_vac0(Population(i),Archive,M);
    end
end

%% �����ⲿ������
function [ new_AC ] = up_vac0( par_eff,old_AC,M )
%par_eff����ѡ����
%old_AC���ⲿ������
%new_AC���µ��ⲿ������

    %��ȡ�ⲿ�����������ӵĸ���
    ss1 = size(old_AC,2);
    old_ACobjs = old_AC.objs;
    par_effobj = par_eff.objs;
    %�ж��������ⲿ�����������ӵ�֧���ϵ
    for i = 1:ss1
        bb1 = 0;
        bb2 = 0;
        for j = 1:M     %Ŀ�꺯���ĸ���M
            aa1 = old_ACobjs(i,j);
            aa2 = par_effobj(1,j);
            if aa2 < aa1
                bb1 = bb1 + 1;
            elseif aa2 == aa1
                bb2 = bb2 + 1;
            end
        end
        %�ж�֧���ϵ
        if bb1 == M         %��ѡ��֧��old_AC(i,:)������ѡ��֧���ⲿ�������еĵ�ǰ����
            old_ACobjs(i,:) = inf;          %ɾ���ⲿ�������б���ѡ��֧�������
        elseif bb2 > 0 && bb1 == M - bb2       %��ѡ��֧��old_AC(i,:)������ѡ��֧���ⲿ�������еĵ�ǰ����
            old_ACobjs(i,:) = inf;
        elseif bb1 + bb2 == 0       %��ѡ�߲�֧��old_AC(i,:)
            par_effobj(1,:) = inf;         %ɾ����ѡ�����ӣ���ѡ�����ӱ���̭�����ܽ����ⲿ������
            break;      %��ǰ��ѡ��֧�����ⲿ�������еĵ�ǰ���ӣ��˳�
        elseif bb2 ~= 0 && bb1 == 0         %��ѡ�߲�֧��old_AC(i,:)
            par_effobj(1,:) = inf;
            break;
        end
    end
    %����ѡ�߼��뵽��ʱ�ⲿ�������У���ʱ�ⲿ�������б�֧����������ڵ��б�Ϊinf���������
    tmp_AC = [par_eff,old_AC]; 
    tmp_ACobjs = [par_effobj;old_ACobjs];
    new_AC = [];
    %�����ϲ������ʱ�ⲿ������
    for i = 1:ss1+1
        %�жϵ�ǰ�����Ƿ�����������֧�䣬�������������µ��ⲿ��������
        if tmp_ACobjs(i,1) ~= inf
            new_AC = [new_AC,tmp_AC(i)];
        end
    end
end
