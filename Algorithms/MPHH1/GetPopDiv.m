function [ PopDiv ] = GetPopDiv( Global,Population )
%% ������Ⱥ����ֵ
	result = tabulate(Population.objs*rand(Global.M,1));
	PopDiv = -sum((result(:,3)*0.01).*log(result(:,3)*0.01));
end