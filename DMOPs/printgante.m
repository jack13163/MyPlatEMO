%% ���Ƹ���ͼ
function [] = printgante( a,title_str )
    %a��ÿһ�еĸ�ʽ��������/�ܵ� �͹� ��ʼʱ�� ����ʱ�� ԭ������
    w=0.5;       %�������
    color=['r','g','b','c','m','y','w'];    %��ɫ

    dset = [];   %������/�ܵ��ļ���
    %ͳ��������/�ܵ��ĸ���
    for i=1:size(a,1)
        if ~ismember(a(i,1), dset)
            dset = [dset, a(i,1)];
        end
    end

    xlabel('time'); 
    ylabel('distillating tower/pipe line'); 
    title(title_str);
    axis([0 270 0 10]); 
    set(gca,'Box','on'); 
    set(gca,'XTick',0:30:240);      %���ÿ̶�
    set(gca,'YTick',0:10/(size(dset,2)+1):10);      %���ÿ̶�
    set(gca,'YTickLabel',{'';num2str((1:size(dset,2)-1)','DS%d');'Pipe';''});       %���ñ��

    for ii=1:size(a,1) 
        %�����������ı߽��
        x=a(ii,[3 3 4 4]);
        y=a(ii,1)*10/(size(dset,2)+1)+[-w/2 w/2 w/2 -w/2];
        %ʹ�ò�ͬ����ɫ����յľ�������
        if a(ii,5)==1
            p=patch('xdata',x,'ydata',y,'facecolor',color(1),'edgecolor','k'); 
        elseif a(ii,5)==2
            p=patch('xdata',x,'ydata',y,'facecolor',color(2),'edgecolor','k'); 
        elseif a(ii,5)==3
            p=patch('xdata',x,'ydata',y,'facecolor',color(3),'edgecolor','k'); 
        elseif a(ii,5)==4
            p=patch('xdata',x,'ydata',y,'facecolor',color(4),'edgecolor','k'); 
        elseif a(ii,5)==5
            p=patch('xdata',x,'ydata',y,'facecolor',color(5),'edgecolor','k'); 
        elseif a(ii,5)==6
            p=patch('xdata',x,'ydata',y,'facecolor',color(6),'edgecolor','k'); 
        elseif a(ii,5)==0
            p=patch('xdata',x,'ydata',y,'facecolor',color(7),'edgecolor','k'); 
        end
        %���þ������������
        text(a(ii,3)+0.5,a(ii,1)*10/(size(dset,2)+1),[num2str(a(ii,2))]); 
    end 
end

