function [ PopDiv ] = GetPopDiv( Global,Population )
%% 计算种群的熵值
	result = tabulate(Population.objs*rand(Global.M,1));
	PopDiv = -sum((result(:,3)*0.01).*log(result(:,3)*0.01));
end