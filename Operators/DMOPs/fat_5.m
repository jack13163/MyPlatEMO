%% 计算解的适应度
% 编码长度为：3*50=150，并行化
function [ eff,scheduleplans ] = fat_5( pop )
    [popsize, w] = size(pop);   %计算种群的规模和决策变量的维数
    eff = pop(:,1:w);
    scheduleplans = cell(1,popsize);    %调度计划
    global iteration;                                                                               %记录函数执行次数
    
    %% 遍历粒子种群
    p = 1;
    while p <= popsize
        %p
       	iteration = 0;                                                                              %重新回溯计数
        x = pop(p,1:w);                                                                          %获取个体
        [RT,DSFR,PIPE1FR,PIPE2FR,FPORDER,FPPackage,c1,c2,c3,TKS,FP,VPIPE,TKHOT,H_DS] = parameter_setting();             %参数设置
        if TKS(TKHOT,3) < VPIPE                                                          %判断加热管道的原油是否足够
            disp('不能完全加热管道，因此不可调度/n');
        end
        
        UD = [1,2,3,4];                                                                          %未实现炼油计划的蒸馏塔集合
        
        %初始调度（使用原油初始库存）
        [~,scheduleplan,FPPackage] = schedule_oil(FPORDER,FPPackage,TKS,DSFR,TKHOT,VPIPE,PIPE2FR);
        
        %首先判断当前状态是否可以调度
        if schedulable(scheduleplan,H_DS, DSFR, RT, UD)
            %回溯调度，并修正不可行解
            [ Flag,x,step,FP,FPPackage,scheduleplan ] = schedule( x,1,TKS,FP,scheduleplan,DSFR,PIPE1FR,PIPE2FR,FPORDER,FPPackage,RT,H_DS );
            if ~Flag
                pop(p,1:w) = rand(1,length(x));                                            %重新初始化个体
                continue;
            else                
                %% 输出方案
                scheduleplans{p} = scheduleplan;
                scheduleplan=sortrows(scheduleplan,1);                
               	%计算适应度函数
               	f1 = gdmix(scheduleplan, c1);                                                 %管道混合成本
              	f2 = gdimix(scheduleplan, c2);                                                %罐底混合成本
              	f3 = gchange(scheduleplan);                                                   %蒸馏塔的供油罐切换次数
               	f4 = gnum(scheduleplan);                                                       %供油罐个数
             	f5 = roundn(energecost(scheduleplan, c3),0);                           %能耗成本
                %更新
                eff(p,1:w) = x(1:w);   
                eff(p,w+1) = f1;
                eff(p,w+2) = f2;
                eff(p,w+3) = f3;
                eff(p,w+4) = f4;
                eff(p,w+5) = f5;
                %接着计算下一个个体
                p = p + 1;
            end
        else
            disp('系统不可调度！\n');
            break;
        end
    end
end

%% 回溯调度(只需要找到一个可行解即可)
function [ Flag,X,step,FP,FPPackage,SchduleTable ] = schedule( X,step,TKS,FP,SchduleTable,DSFR,PIPE1FR,PIPE2FR,FPORDER,FPPackage,RT,H_DS )
% 只需要依次遍历下去，若出现不可调度的地方就回溯
% 然后接着遍历，直到FP为空，代表找到一个可行解
    Flag = false;%是否找到可行解的标志
    %保存旧值
    Old_step = step;
    Old_FP = FP;
    Old_FPPackage = FPPackage;
    Old_SchduleTable = SchduleTable;
    times = 100;                                                                                   %函数最大执行次数
    global iteration;                                                                               %记录函数执行次数
   	iteration = iteration + 1;
    if iteration > times
        return;
    end
    
    %结束条件
    if round(sum(FP)) == 0
        Flag = true;
    elseif step <= 50   %步长不得超过50
    	%计算ET和UD(随后只需要检查冲突即可)
        T = max(SchduleTable(SchduleTable(:,1) > 4,4));                     %插入决策
        if isempty(T)
            T = 0;
        end
      	[ET,UD] = get_et_ud(TKS,SchduleTable,FPORDER,FPPackage,T);
        %记录遍历的足迹
        footprint = zeros(length(ET) + 1, length(UD));
        
    	%深度优先遍历
        while ~all(all(footprint)) && ~Flag && iteration <= times
           	TK_NO=getInt(X(3*step-2), size(ET, 2));                                                                                         %供油罐序号
            DS_NO = getInt(X(3*step-1), size(UD, 2)-1)+1;                                                                               %蒸馏塔
          	DS = UD(DS_NO);
           	R = getInt(X(3*step), 3-1)+1;                                                                                                         %管道转运速率
            COT = FPORDER(DS,find(FPPackage(DS,:), 1 ));                                                                         %确定原油类型   
          	PET1 = max(SchduleTable(SchduleTable(:,1) == 5, 4));                                                                 %计算管道1结束转运时间
            if isempty(PET1)
                PET1 = 0;
            end
         	PET2 = max(SchduleTable(SchduleTable(:,1) == 6, 4));                                                                 %计算管道2结束转运时间
          	DSFET = max(SchduleTable(SchduleTable(:,1) == DS, 4));                                                            %计算DS塔炼油结束时间
            
           	f1 = false;                                                                                                                                           %当前转运的安全油量过低时，不进行转运
            
            %为了避免回溯过深造成的计算代价过高，因此，采用先调度高熔点原油再调度低熔点原油的策略
            if ismember(H_DS,UD) && DS ~= H_DS
                f1 = false;
            elseif DS == H_DS
                %高熔点转运及其炼油
                if TK_NO ~= 0                                                                                                                                %当供油罐选择不为0，代表需要进行转运操作，否则，不进行任何操作
                    TK = ET(TK_NO) ;                                                                                                                      %确定供油罐
                    %计算当前PET1时刻原油库存总油量
                    total = gettotal(DSFR,SchduleTable,PET2);
                    %安全的原油体积
                    V = getVolume(DS,TK,UD,H_DS,RT,total,DSFR,PIPE1FR(R),PIPE2FR(2),SchduleTable);
                    
                    %体积不可以低于特定的阈值
                    if V >= 15000
                        if TKS(TK, 1) < V
                            V = TKS(TK, 1);
                        end
                        if FPPackage(DS,find(FPPackage(DS,:), 1 )) < V
                            V = FPPackage(DS,find(FPPackage(DS,:), 1 ));
                        end

                        %进料包
                        FP(COT) = FP(COT) - V;
                        FPPackage(DS,find(FPPackage(DS,:), 1 )) = FPPackage(DS,find(FPPackage(DS,:), 1 )) - V;
                     	SchduleTable = [SchduleTable; 6, TK, PET2, PET2+V/PIPE2FR(2), COT,R,V];                               %管道2转运记录(最低速转运)
                       	SchduleTable = [SchduleTable; DS, TK, DSFET, DSFET + V/DSFR(DS), COT,2,V];                        %管道2炼油记录
                        f1 = true;
                    end
                end
            else
                if TK_NO==0
                    PipeStoptime  = roundn(getPipeStoptime(DSFR, RT, H_DS,SchduleTable,TKS,PET1),-1);%四舍五入，取6位小数
                    if PipeStoptime > 0
                        %停运(尽量不停运)
                        SchduleTable = stop(SchduleTable,PipeStoptime);
                        f1 = true;
                    end
                else
                    TK = ET(TK_NO) ;                                                           %确定供油罐
                    %计算当前PET1时刻原油库存总油量
                    total = gettotal(DSFR,SchduleTable,PET1);
                    %安全的原油体积
                    V = getVolume(DS,TK,UD,H_DS,RT,total,DSFR,PIPE1FR(R),PIPE2FR(2),SchduleTable);
                    if V >= 15000
                        if TKS(TK, 1) < V
                            V = TKS(TK, 1);
                        end
                        if FPPackage(DS,find(FPPackage(DS,:), 1 )) < V
                            V = FPPackage(DS,find(FPPackage(DS,:), 1 ));
                        end

                        %进料包
                        FP(COT) = FP(COT) - V;
                        FPPackage(DS,find(FPPackage(DS,:), 1 )) = FPPackage(DS,find(FPPackage(DS,:), 1 )) - V;
                       	%更新转运记录和炼油记录
                     	SchduleTable = [SchduleTable;  length(DSFR)+1, TK, PET1, PET1 + V / PIPE1FR(R), COT,R,V];     %只记录使用的速率序号
                       	SchduleTable = [SchduleTable;  DS, TK, DSFET, DSFET + V / DSFR(DS), COT,R,V];
                       	f1 = true;
                    end
                end
            end
            
                    
         	%判断试调度是否成功，下一状态是否可以调度
          	if f1 && conflic_detect(SchduleTable) && schedulable(SchduleTable, H_DS, DSFR, RT, UD)
                f2 = true;
            else
                f2 = false;
            end
            
          	%可调度，进行下一步
          	if f2
              	step = step + 1;
                [ Flag,X,step,FP,FPPackage,SchduleTable ] = schedule( X,step,TKS,FP,SchduleTable,DSFR,PIPE1FR,PIPE2FR,FPORDER,FPPackage,RT,H_DS );
            end
                        
            %若未找到可行解，标记并回滚
            if ~Flag && iteration <= times
%                 clc;
%                 step
%                 FPPackage

                footprint(TK_NO + 1,DS_NO) = 1;                                                 %先尝试不同的塔，再尝试不同的罐
                %数据回滚
                [step,FP,FPPackage,SchduleTable] = rollback(Old_step,Old_FP,Old_FPPackage,Old_SchduleTable);
            
               	%更改策略（先试试是否可以停运）
                f_j = false;
                if ~all(all(footprint))                                                                         %还存在未尝试过的路径
                    for i = 1:size(footprint,1)
                        if i == 1 && any(footprint(i,:))                                                  %尝试过停运后不再继续尝试
                            footprint(i,:) = 1;
                        end
                        if ~all(footprint(i,:))                                                                 %当前供油罐还存在未尝试过的蒸馏塔
                            while TK_NO + 1 ~= i
                                X(3*step-2) = rand;
                                TK_NO = getInt(X(3*step-2),size(ET, 2));
                            end
                            for j = 1:size(footprint,2)
                                if footprint(i,j) == 0                                                              %未尝试过的路径
                                    while DS_NO ~= j
                                        X(3*step-1) = rand;
                                        DS_NO = getInt(X(3*step-1),size(UD, 2) - 1) + 1;
                                    end
                                    f_j = true;
                                    break;
                                end
                            end
                            if f_j
                                break;
                            end
                        end
                    end
                end
            end
        end
    end
end

%判断是否存在矛盾的调度(flag==true，代表不冲突)
function [ flag,TK ] = conflic_detect(SchduleTable)
    flag = true;
    TKs = unique(SchduleTable(:,2));
    for i=1:length(TKs)
        TK = TKs(i);
        if TK == 0
            continue;
        end
        Tmp_SchduleTable = SchduleTable(SchduleTable(:,2) == TK,:);
        s = 1;
        k = 1;
        t = [];                                                                                 %区间
        while s <= size(Tmp_SchduleTable,1)
            if s <= size(Tmp_SchduleTable,1) - 1
                if Tmp_SchduleTable(s,1) > 4 && Tmp_SchduleTable(s+1,1) <= 4
                    t(k,1) = Tmp_SchduleTable(s,3);
                    t(k,2) = Tmp_SchduleTable(s+1,4);
                    s = s + 2;
                else                
                    t(k,1) = Tmp_SchduleTable(s,3);
                    t(k,2) = Tmp_SchduleTable(s,4);
                    s = s + 1;
                end
            else
                t(k,1) = Tmp_SchduleTable(s,3);
                t(k,2) = Tmp_SchduleTable(s,4);
               	s = s + 1;
            end
            k = k + 1;
        end
        
        %检查区间是否冲突
        if size(t,1) > 1
            for m=1:size(t,1) - 1
                for n=m+1:size(t,1)
                    t1 = t(m,:);
                    t2 = t(n,:);
                    if max(t1(1),t2(1)) < min(t1(2),t2(2))
                        flag = false;           %检查到冲突，退出
                        return;
                    end
                end
            end
        end
    end
end

%% 判断当前系统状态是否可调度(独立检查低熔点蒸馏塔油量是否满足下限，同时判断是否冲突)
function [ flag ] = schedulable(SchduleTable, H_DS, DSFR, RT, UD)
    flag = true;
    %分离出高熔点塔和低熔点塔
    if isempty(UD)        
        disp('ud为空');
    end
    LT_UD = UD(UD~=H_DS);
    HT_UD =H_DS;
    
    %判断低熔点塔是否满足最低油量下限
    if ~isempty(LT_UD)
        PET1 = max(SchduleTable(SchduleTable(:,1) == 5,4));                                   %管道1的转运结束时间
        if isempty(PET1)
            PET1 = 0;
        end
        total = gettotal(DSFR,SchduleTable,PET1);                                                   %库存总油量
        
        for i=1:length(LT_UD)
            DSN = LT_UD(i);
            arfai = RT * DSFR(DSN);

            if round(total(DSN)) < length(LT_UD) * arfai                                                %判断是否满足最低油量约束
                flag=false;
                %disp('低熔点原油油量不足！');
                return;
            end
        end
    end
    
    %判断高熔点塔是否满足最低油量下限(能够使得输入的油和输出的油相等即可)
    if ~isempty(HT_UD)
        PET2 = max(SchduleTable(SchduleTable(:,1) == 6,4));                                   %管道2的转运结束时间
        if isempty(PET2)
            PET2 = 0;
        end
        total = gettotal(DSFR,SchduleTable,PET2);                                                   %库存总油量
    
        if length(HT_UD) > 1                                                                                    %高熔点塔多于1个的情况暂时不考虑
            flag = false;
            disp('高熔点塔多于1个的情况暂时不考虑！');
            return;
        end
        if total(HT_UD) < RT * DSFR(HT_UD);                %低速转运
            %disp('高熔点原油油量不足');
            flag = false;
            return;
        end
    end
end

%% 安全油量
function [V] = getVolume(DS,TK,UD,H_DS,RT,total,DSFR,PIPE1FR,PIPE2FR,SchduleTable)
    %高熔点原油和低熔点原油都要考虑，其中，PIPE1FR,PIPE2FR为一个数字，而非数组
    if isempty(UD)        
        disp('ud为空');
    end
    LT_UD = UD(UD~=H_DS);
    HT_UD =H_DS;
    
    if DS ~= H_DS
        if ~isempty(LT_UD)
            %计算满足供油罐占用时间不冲突的最大体积
            Tmp_SchduleTable = SchduleTable(SchduleTable(:,2) == TK,:);
            if size(Tmp_SchduleTable,1) > 0
                s = 1;
                k = 1;
                interval = [];                                                                                 %区间
                while s <= size(Tmp_SchduleTable,1)
                    if s <= size(Tmp_SchduleTable,1) - 1
                        if Tmp_SchduleTable(s,1) > 4 && Tmp_SchduleTable(s+1,1) <= 4
                            interval(k,1) = Tmp_SchduleTable(s,3);
                            interval(k,2) = Tmp_SchduleTable(s+1,4);
                            s = s + 2;
                        else                
                            interval(k,1) = Tmp_SchduleTable(s,3);
                            interval(k,2) = Tmp_SchduleTable(s,4);
                            s = s + 1;
                        end
                    else
                        interval(k,1) = Tmp_SchduleTable(s,3);
                        interval(k,2) = Tmp_SchduleTable(s,4);
                        s = s + 1;
                    end
                    k = k + 1;
                end
                PET1 = max(SchduleTable(SchduleTable(:,1) == 5,4));                                   %管道1的转运结束时间
                if isempty(PET1)
                    PET1 = 0;
                end

                for i=size(interval,1)
                    if interval(i,1) <= PET1 && interval(i,2) > PET1                                            %检查是否冲突，则安全体积为0
                        V = 0;
                        return;
                    end
                end
            end
            
            %判断低熔点塔是否满足最低油量下限
            for i=1:length(LT_UD)
                DSL = LT_UD(i);                                                                                             %低熔点蒸馏塔
                arfai = RT * DSFR(DSL);                                                                                 %低熔点蒸馏塔的原油体积下限
              	Vr = round(total(DSL)) - length(LT_UD) * arfai;
                t = Vr/DSFR(DSL);
            end
            
            V = round(t * PIPE1FR);
        end
    else
        %判断高熔点塔是否满足最低油量下限(能够使得输入的油和输出的油相等即可)
        if ~isempty(HT_UD)
            if length(HT_UD) > 1                                                                                    %高熔点塔多于1个的情况暂时不考虑
                disp('高熔点塔多于1个的情况暂时不考虑！');
                return;
            end
            
            %计算满足供油罐占用时间不冲突的最大体积
        	Tmp_SchduleTable = SchduleTable(SchduleTable(:,2) == TK,:);
            if size(Tmp_SchduleTable,1) > 0
                s = 1;
                k = 1;
                interval = [];                                                                                 %区间
                while s <= size(Tmp_SchduleTable,1)
                    if s <= size(Tmp_SchduleTable,1) - 1
                        if Tmp_SchduleTable(s,1) > 4 && Tmp_SchduleTable(s+1,1) <= 4
                            interval(k,1) = Tmp_SchduleTable(s,3);
                            interval(k,2) = Tmp_SchduleTable(s+1,4);
                            s = s + 2;
                        else                
                            interval(k,1) = Tmp_SchduleTable(s,3);
                            interval(k,2) = Tmp_SchduleTable(s,4);
                            s = s + 1;
                        end
                    else
                        interval(k,1) = Tmp_SchduleTable(s,3);
                        interval(k,2) = Tmp_SchduleTable(s,4);
                        s = s + 1;
                    end
                    k = k + 1;
                end
                PET2 = max(SchduleTable(SchduleTable(:,1) == 6,4));                                   %管道2的转运结束时间
                if isempty(PET2)
                    PET2 = 0;
                end
                
                for i=size(interval,1)
                    if interval(i,1) <= PET2 && interval(i,2) > PET2
                        V = 0;
                        return;
                    end
                end
            end
    
            %计算库存原油体积下限
          	DSH = HT_UD(1);
           	Vr = total(DSH) - RT * DSFR(DSH);                %低速转运
         	t = Vr/DSFR(DSH);
                        
            %计算最大转运体积
            V = round(t * PIPE2FR);
        end
    end
end

%% 停运
function [SchduleTable] = stop(SchduleTable,PipeStoptime)
	%停运期间，各个蒸馏塔是否炼油结束
	PETOLD = max(SchduleTable(SchduleTable(:,1) == 5, 4));                                                         %计算低熔点转运结束时间
    if isempty(PETOLD)
        PETOLD = 0;
    end
    tmp = ones(1,size(SchduleTable,1))*inf;
    for i=1:size(SchduleTable,1)
        if SchduleTable(i,1) <= 4 && SchduleTable(i,4) > PETOLD && SchduleTable(i,4) <= PETOLD + PipeStoptime               %停运
            tmp(i) = SchduleTable(i,4);
        end
    end
    [m,n] = min(tmp);
    if m ~= inf
        tk = SchduleTable(n,2);
    else
        tk = 0;
    end
    
	%有油罐释放
	if tk ~= 0
        %更新转运结束时间为供油罐可用的时间
        PET = m;
    else
        %正常停运
        PET = PETOLD + PipeStoptime;
	end
	%更新转运记录表，不需要挑选供油罐供油
	SchduleTable = [SchduleTable; 5, 0, PETOLD, PET, 0,0,0];                        %这里5代表管道1停运
end

%% 根据调度记录，获取ET和UD
function [ET,UD] = get_et_ud(TKS,scheduletable,FPORDER,FPPackage,t)
    ET = [];
    UD = [];
    %查找可用罐
    for i=1:size(TKS,1)
        tmp = scheduletable(scheduletable(:,2)==i,:);
        
        %尚未使用的供油罐
        if isempty(tmp)
            ET = [ET,i];
            continue;
        end
        
        %原油初始库存还没用完
        if tmp(1,6) == 0 && tmp(1,4) > t
            continue;
        end
        
        %判断时刻t是否在使用
        flag = true;
        for j=1:size(tmp,1)
            %转运记录后面是炼油记录，则两者之间也不容许使用。
            if j + 1 <= size(tmp,1)
                if tmp(j,1) > 4 && tmp(j+1) <= 4            %其中4代表蒸馏塔的个数
                    if tmp(j,3) <= t && tmp(j+1,4) >= t
                        flag = false;
                        break;
                    end
                end
            end
            
            %单独的转运记录或者供油记录，则只需要检查t是否处于其开始时间和结束时间即可。
            if tmp(j,3) <= t && tmp(j,4) >= t
                flag = false;
                break;
            end
        end
        
        if flag
            ET = [ET,i];
        end
    end
    
    %检测未完成炼油任务的蒸馏塔
    for i=1:size(FPORDER,1)
        if round(sum(FPPackage(i,:))) ~= 0
            UD = [UD,i];            
        end
    end
end

%% 调度最初库存中存在的原油
function [total,scheduleplan,FPPackage,PET1,PET2] = schedule_oil(FPORDER,FPPackage,TKS,DSFR,TKHOT,VPIPE,PIPE2FR)
    scheduleplan = [];                                        %调度计划
    DSFET = zeros(1,size(FPORDER,1));           %各个塔最后一次供油结束时间
    total = zeros(1,size(FPORDER,1));               %各个塔总的油量
    PET1 = 0;
    PET2 = 0;
    
    %计算各个塔的原油总量，最后一次供油结束时间，以及初始调度计划
    for i=1:size(TKS,1)  
        if TKS(i,2)~=inf                                                              %过滤出空罐
            if i == TKHOT                                                            %盛放加热管道所需要的原油的供油罐
                PET2 = TKS(i,3) / PIPE2FR(1);                               %高速加热管道2
                scheduleplan = [scheduleplan; 6, i, 0, PET2, TKS(i,2),0,TKS(i,3)];                                          %记录管道2加热管道记录(高速)
                scheduleplan = [scheduleplan; 6, i, PET2, PET2 + VPIPE/PIPE2FR(1), TKS(i,2),0,VPIPE];     %记录管道2将低熔点原油挤出管道记录（）
                PET2 = PET2 + VPIPE/PIPE2FR(1);                       %挤出低熔点原油
                TKS(i,3) = VPIPE;                                                  %继续供油
            end
           	OT = TKS(i,2);                                    %原油类型
           	[m,n] = find(FPORDER == OT);

          	%选择蒸馏塔供油
          	if length(m) == 1 && length(n) == 1
            	ds = m;
            else
                [~,ind] = min(DSFET(m));                %为最需要油的塔供油（贪心规则：在两个塔同时需要一种油时，先满足最需要的塔的需求）
               	ds = m(ind);
            end

            %获取ds塔的炼油顺序，并计算FPPackage
            DS_FPORDER = FPORDER(ds,:);
            for j=1:length(DS_FPORDER)
                if OT == DS_FPORDER(j)
                    FPPackage(ds,j) = FPPackage(ds,j) - TKS(i,3);
                    break;
                end
            end
            
         	total(ds) = total(ds) + TKS(i,3);           %计算各个蒸馏塔的总油量
           	Temp = DSFET(ds);                           %保存前一次的供油结束时间
          	%计算供油罐的供油结束时间，并添加到供油罐对应的供油记录表中
           	DSFET(ds) = DSFET(ds) + TKS(i,3) / DSFR(ds);      %DSFET(i)代表蒸馏塔i的最后一次供油结束时间，TKS(k,3)代表供油罐k的当前油量，DSFR(i)代表蒸馏塔i的炼油速率
          	scheduleplan = [scheduleplan; ds, i, Temp, DSFET(ds), TKS(i,2),0,TKS(i,3)];     %记录初始调度计划
        end
   	end
end

%% 参数设置
function [RT,DSFR,PIPE1FR,PIPE2FR,FPORDER,FPPackage,c1,c2,c3,TKS,FP,VPIPE,TKHOT,H_DS] = parameter_setting()
 	%% 参数设置(不变参数)
   	RT = 4;                                                                          %驻留时间
	DSFR = [250,334,250,625];                                             %蒸馏塔炼油速率(可计算蒸馏塔个数)
    PIPE1FR = [833.3 1250 1375];                                        %管道1的输油速率
    PIPE2FR = [625,420];                                                     %管道2的输油速率
    FPORDER = [3,4,8,8;1,2,11,11;5,4,6,10;7,8,9,9];              %蒸馏塔的炼油顺序
    Plan = [25000,24000,44055,0;16000,26000,79285,0;15000,17000,18000,43055;22000,103000,107638,0];            %蒸馏塔的蒸馏包
    VPIPE = 18000;                                                              %管道的体积
    TKHOT = 8;                                                                    %装有加热原油的供油罐（假设该供油罐中的原油体积已经足够加热管道）
    H_DS = 2;                                                                      %标识高熔点塔，即代表2号塔为高熔点塔
    
	c1 = [0	7	25	27	26	28	23	20	21	24	10;
        7	0	29	26	30	26	24	27	28	25	12;
        25	29	0	11	12	13	10	15	16	7	31;
        27	26	11	0	11	12	13	10	15	16	28;
        26	30	12	11	0	10	12	13	10	15	25;
        28	26	13	12	10	0	11	12	13	17	27;
        23	24	10	13	12	11	0	11	12	15	32;
        20	27	15	10	13	12	11	0	9	18	35;
        21	28	16	15	10	13	12	9	0	13	21;
        24	25	7	16	15	17	15	18	13	0	24;
        10	12	31	28	25	27	32	35	21	24	0];                                                    %管道混合成本
  	c2 = [0	7	25	27	26	28	23	20	21	24	10;
        7	0	29	26	30	26	24	27	28	25	12;
        25	29	0	11	12	13	10	15	16	7	31;
        27	26	11	0	11	12	13	10	15	16	28;
        26	30	12	11	0	10	12	13	10	15	25;
        28	26	13	12	10	0	11	12	13	17	27;
        23	24	10	13	12	11	0	11	12	15	32;
        20	27	15	10	13	12	11	0	9	18	35;
        21	28	16	15	10	13	12	9	0	13	21;
        24	25	7	16	15	17	15	18	13	0	24;
        10	12	31	28	25	27	32	35	21	24	0];                                                 %罐底混合成本
	c3 = [1 2 3];                                                                   %能耗成本
    
	%% 参数设置(可变参数)
	TKS = [
            20000,inf,0,0,0,0;
            34000,1,16000,0,0,0;
            34000,2,26000,0,0,0;
            25000,inf,0,0,0,0;
            30000,3,25000,0,0,0;
            34000,5,15000,0,0,0;
            20000,4,17000,0,0,0;
            30000,6,25000,0,0,0;
            34000,7,22000,0,0,0
            34000,8,19000,0,0,0
            20000,inf,0,0,0,0
            ];                                        %油罐初始库存，后面会自动添加蒸馏塔和供油开始时间和结束时间
	EO = TKS(TKS(:,2)~=inf,:);             %所有的非空油罐
	TKO = zeros(1,length(unique(FPORDER)));
	for i=1:size(EO,1)
        k = EO(i,2);
        TKO(k) = TKO(k) + EO(i,3);
	end
	%Total = [16000 26000 25000 41000 15000 18000 22000 147055 107638 43055 79285];
    Total = zeros(size(TKO));
    for i=1:size(FPORDER,1)
        for j=1:size(FPORDER,2)
            OT = FPORDER(i,j);
            Total(OT) = Total(OT) + Plan(i,j);
        end
    end
	FP = Total - TKO;               %进料包(可计算原油种类)
    FP(FP<0)=0;                      %去除库存中多余的油（即负值）
    FPPackage = Plan;            %暂时令FPPackage==Plan
end

%% 数据回滚  [step,FP,SchduleTable] = rollback(Old_step,Old_FP,Old_SchduleTable);
function [step,FP,FPPackage,SchduleTable] = rollback(Old_step,Old_FP,Old_FPPackage,Old_SchduleTable)
        step = Old_step;
        FP = Old_FP;
        FPPackage = Old_FPPackage;
        SchduleTable = Old_SchduleTable;
end

%% 计算可用原油总量
function [total] = gettotal(DSFR,ScheduleTable,T)
   	total = zeros(1,length(DSFR));
    %根据调度计划，计算现有库存
    for i=1:size(ScheduleTable,1)
        if ScheduleTable(i,1) <= length(total)
            ds = ScheduleTable(i,1);
            if T <  ScheduleTable(i,4)
                if T <= ScheduleTable(i,3);
                    total(ds) = total(ds) + ScheduleTable(i,7);
                else
                    total(ds) = total(ds) + DSFR(ds) * (ScheduleTable(i,4) - T);
                end
            end
        end
    end
end

%% 供油罐切换次数
function [ x ] = gchange( a )
    K = 0;
    
    %获取转运记录
    for i=1:size(a,1)
        if a(i,1) >= 5
            K = i;
            break;
        end
    end
    x = K - 1;
end

%% 罐底混合次数
function [ x ] = gdimix( a, c2 )
%c2：罐底混合成本
    m2 = zeros(size(c2));
    K = 0;
    
    %获取转运记录
    for i=1:size(a,1)
        if a(i,1) >= 5
            K = i;
            break;
        end
    end
    K = K - 1;
    record = a(1:K,:);
    record = sortrows(record,2);
    
    %统计各个油罐的使用次数
    A = record(:,2)';
    B = unique(A);
    tks = [];
    for i=1:size(B,2)
        tks(i) = sum(A(:)==B(i));
    end
    
	%根据炼油记录表分析罐底混合次数
    for i=1:size(tks,2)     %i代表油罐
        tmp = 0;
        for j=1:(size(record,1) - 1)      %j代表炼油记录
            if record(j,2) == i
                tmp =  tmp + 1;
                if tmp >= tks(i)     %tks(i)：油罐i中的使用次数
                    break;
                end
                if record(j,5) ~= record(j+1,5)
                    m2(record(j,5), record(j+1,5)) = m2(record(j,5), record(j+1,5)) + 1;
                end
            end
        end
    end
    x = sum(sum(m2.*c2));          %罐底混合成本
end

%% 获取管道混合成本
function [ x ] = gdmix( a, c1 )
%c1：管道混合成本
    m1 = zeros(size(c1));
    K = 0;
    
    %获取转运记录
    for i=1:size(a,1)
        if a(i,1) >= 5
            K = i;
            break;
        end
    end
    PIPE = a(K:size(a,1),:);
    PIPE = PIPE(PIPE(:,5)~=0,:);   %删除其中的停运记录
    PIPE1 = PIPE(PIPE(:,1)==5,:);   %管道1的停运记录
    PIPE2 = PIPE(PIPE(:,1)==6,:);   %管道2的停运记录
    %分别计算管道1和管道2的原油混合成本
    for i=1:size(PIPE1,1)-1
        if PIPE1(i,5) ~= PIPE1(i+1,5)
            m1(PIPE1(i,5),PIPE1(i+1,5)) = m1(PIPE1(i,5),PIPE1(i+1,5)) +1;
        end
    end
	for i=1:size(PIPE2,1)-1
        if PIPE2(i,5) ~= PIPE2(i+1,5)
            m1(PIPE2(i,5),PIPE2(i+1,5)) = m1(PIPE2(i,5),PIPE2(i+1,5)) +1;
        end
	end
    x = sum(sum(m1.*c1));       %管道混合次数
end

%% 供油罐个数
function [ x ] = gnum( a )
    TKSet = [];
    for i=1:size(a,1)
        if a(i,2) ~= 0 && ~ismember(a(i,2), TKSet)
            TKSet = [TKSet, a(i,2)];
        end
    end
    x = size(TKSet,2);      %供油罐的个数
end

%% 能耗成本
function [ x ] = energecost( a,c )
    pipplan = a(a(:,1) >= 5 & a(:,6) ~= 0,:);
    x = c(pipplan(:,6))*(pipplan(:,4)-pipplan(:,3));       %时间*单位时间能耗成本
end

%% 根据[0~1]之间的随机数，返回一个[0~length]整数
function [ r ] = getInt( a, length )
    %a：[0~1]
    if a == 1
        r = length;
    else
        r = fix((length + 1) * a);
    end
end

%% 获取管道停运安全时间
function [ t ] = getPipeStoptime(DSFR, RT, H_DS,ScheduleTable,TKS,T)
    t=inf;
    [total] = gettotal(DSFR,ScheduleTable,T);
    
    %油量约束(高熔点原油和低熔点原油分开考虑)
    DSFR_L = DSFR;
    DSFR_L(H_DS) = [];                                                                %低熔点蒸馏塔的炼油速率
    total(H_DS) = [];
    
    %检查是否有停运的必要
    if length(unique(ScheduleTable(:,2))) < size(TKS,1)                     %存在没有使用的供油罐，因此，没有停运的必要
        t = 0;
        return;
    else
    	TKs = unique(ScheduleTable(:,2));
        Temp_ScheduleTable = ScheduleTable(ScheduleTable(:,1)<=4,:);
        for i=1:length(TKs)
            TK = TKs(i);
            if max(Temp_ScheduleTable(Temp_ScheduleTable(:,2)==TK,4)) < T            %存在供油罐，即不需要转运
                t = 0;
                return;
            end
        end
    end
    
    %低熔点原油油量约束
    K = length(DSFR_L);
    for i=1:K
        arfai = RT * DSFR_L(i);

        %预留一个调度周期的原油
      	total(i) = total(i) - K * arfai;

        tmp = total(i) / DSFR_L(i);
        if tmp < t
            t = tmp;
        end
    end
end