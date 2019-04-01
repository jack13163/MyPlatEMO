function [ crowd_value ] = CrowdingDistance( Population )
%������Ⱥ��Ԫ�ص�ӵ������ֵ
%Population����Ⱥ

    num = size(Population,2);   %��Ⱥ��ģ
    objs = Population.objs;
    M = size(objs,2); %Ŀ�����
    crowd_value = zeros(1,num);
    lim_f = zeros(M,2);

    %���ӵ��������ڵ���2���ܼ���ӵ������ֵ������ӵ������ֵȡ1
    if num >= 2
        for i = 1:M
            lim_f(i,1) = min(objs(:,i));         %Ŀ�꺯��i�ļ�Сֵ
            lim_f(i,2) = max(objs(:,i));        %Ŀ�꺯��i�ļ���ֵ
        end
        DD = zeros(num,M);
        %����Ŀ�꺯��
        for i = 1:M
            [~,ind] = sort(objs(:,i));        %��Ŀ�꺯��i��ֵ����ind��������val��ֵ
            %�����������е�����
            for j = 1:num
                if j == 1
                    DD(ind(j),i) = inf;     %����Ŀ�꺯��i�ļ�Сֵ��ӵ������ֵΪ�����
                elseif j == num
                    DD(ind(j),i) = inf;     %����Ŀ�꺯��i�ļ���ֵ��ӵ������ֵΪ�����
                else
                    DD(ind(j),i) = objs(ind(j+1),i) - objs(ind(j-1),i)/2*(lim_f(i,2) - lim_f(i,1));     %����Ŀ�꺯��i�Ǳ߽���ӵ������ֵ
                end
            end
        end
        for jj = 1:num
            crowd_value(jj) = sum(DD(jj,:));        %����������ӵ�ӵ������ֵ
        end
    else
        crowd_value(1) = 1;
    end
end