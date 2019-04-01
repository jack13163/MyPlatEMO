function [ A1 ] = StoreArchive( Global,Archive )
%(合并操作后)将储备集保存到mat文件中

    fileName = sprintf('%s\\Problems\\DMOPs\\PF.mat',pwd);%PF
    tmpfileName = sprintf('%s\\Problems\\DMOPs\\PF_Archive.mat',pwd);%储备集
    %判断mat文件是否存在
    if exist(tmpfileName,'file')
        A1 = [];
        load(tmpfileName,'A1');
        if exist('A1','var')
            %若已经存在储备集，则需要执行合并操作
            A1 = UpdateArchive(Archive,A1,Global.N,Global.D);
        else
            A1 = Archive;
        end
    else        
        A1 = Archive;
    end
    save(tmpfileName,'A1');%保存储备集
    %目标的标准化
	Parameter = {A1.objs};
	if exist(fileName,'file') == 2
        save(fileName,'Parameter');%保存PF
	end
end

