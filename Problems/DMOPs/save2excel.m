%% 保存数据到excel
function [] = save2excel(data,title,filename)
    [m, n] = size(data);
    data_cell = mat2cell(data, ones(m,1), ones(n,1));    % 将data切割成m*n的cell矩阵
    result = [title; data_cell];                                            % 将变量名称和数值组集到result
    xlswrite(filename, result);                                      % 将result写入到wind.xls文件中
end