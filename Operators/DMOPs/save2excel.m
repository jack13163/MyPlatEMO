%% �������ݵ�excel
function [] = save2excel(data,title,filename)
    [m, n] = size(data);
    data_cell = mat2cell(data, ones(m,1), ones(n,1));    % ��data�и��m*n��cell����
    result = [title; data_cell];                                            % ���������ƺ���ֵ�鼯��result
    xlswrite(filename, result);                                      % ��resultд�뵽wind.xls�ļ���
end