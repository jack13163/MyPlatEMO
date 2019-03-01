function [  ] = Draw_bubble( data )
%% ��������ͼ
    x = data(:,1)';
    y = data(:,2)';
    z = data(:,3)';
    
    %�����ܺĵȼ�
 	energe_level = unique(data(:,5));
    level = zeros(1,size(data,1));
    for i=1:size(data,1)
        for j=1:length(energe_level)
            if data(i,5) == energe_level(j)
                k=j;
                break;
            end
        end
        level(i) = k;
    end
    sz_base = 90;
    u = 50 + sz_base * level;
    
    %������ɫ
 	gnums = unique(data(:,4));
	c = parula(length(gnums));    
    c = fliplr(c')';
    
    %����ʹ�õĹ޵ĸ�������ͼ��
    figure;
    g_num = data(:,4);
    legendtext = cell(1,length(gnums));
    for i=1:length(gnums)
        ind = find(g_num== gnums(i));  %���͹޸���Ϊtmp(i)�ļ�¼�±�
        scatter3(x(ind),y(ind),z(ind),u(ind),c(i,:),'fill');%��ɫ���ù޸��������ݴ�С���л�����
        legendtext{i} = num2str(gnums(i));
        hold on;
    end
    legend(legendtext);
    title('Visible of Prato front');    
    range = [min(x), max(x), min(y), max(y), min(z), max(z)];
    axis(range);% ���������ʾ��Χ
    xlabel('Pipeline mixing cost');% x�����ơ����弰��С
    ylabel('Tank bottom mixing cost');% y�����ơ����弰��С
    zlabel('Number of switchs');% z�����ơ����弰��С
%    text(max(x), min(y), max(z),'���ݵ����Խ�󣬴����ܺ�Խ��','fontsize',12);
%     colorbar;%��ʾcolorbar
%     caxis([6.5,10.5]);%����colorbar�ķ�Χ
%     colormap(c);%����colorbar����ɫ
end