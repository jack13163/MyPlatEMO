function [ A1 ] = StoreArchive( Global,Archive )
%(�ϲ�������)�����������浽mat�ļ���

    fileName = sprintf('%s\\Problems\\DMOPs\\PF.mat',pwd);%PF
    tmpfileName = sprintf('%s\\Problems\\DMOPs\\PF_Archive.mat',pwd);%������
    %�ж�mat�ļ��Ƿ����
    if exist(tmpfileName,'file')
        A1 = [];
        load(tmpfileName,'A1');
        if exist('A1','var')
            %���Ѿ����ڴ�����������Ҫִ�кϲ�����
            A1 = UpdateArchive(Archive,A1,Global.N,Global.D);
        else
            A1 = Archive;
        end
    else        
        A1 = Archive;
    end
    save(tmpfileName,'A1');%���洢����
    %Ŀ��ı�׼��
	Parameter = {A1.objs};
	if exist(fileName,'file') == 2
        save(fileName,'Parameter');%����PF
	end
end

