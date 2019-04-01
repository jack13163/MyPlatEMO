%% ��������Ӧ��
% ���볤��Ϊ��3*50=150�����л�
function [ eff,scheduleplans ] = fat_5( pop )
    [popsize, w] = size(pop);   %������Ⱥ�Ĺ�ģ�;��߱�����ά��
    eff = pop(:,1:w);
    scheduleplans = cell(1,popsize);    %���ȼƻ�
    global iteration;                                                                               %��¼����ִ�д���
    
    %% ����������Ⱥ
    p = 1;
    while p <= popsize
        %p
       	iteration = 0;                                                                              %���»��ݼ���
        x = pop(p,1:w);                                                                          %��ȡ����
        [RT,DSFR,PIPE1FR,PIPE2FR,FPORDER,FPPackage,c1,c2,c3,TKS,FP,VPIPE,TKHOT,H_DS] = parameter_setting();             %��������
        if TKS(TKHOT,3) < VPIPE                                                          %�жϼ��ȹܵ���ԭ���Ƿ��㹻
            disp('������ȫ���ȹܵ�����˲��ɵ���/n');
        end
        
        UD = [1,2,3,4];                                                                          %δʵ�����ͼƻ�������������
        
        %��ʼ���ȣ�ʹ��ԭ�ͳ�ʼ��棩
        [~,scheduleplan,FPPackage] = schedule_oil(FPORDER,FPPackage,TKS,DSFR,TKHOT,VPIPE,PIPE2FR);
        
        %�����жϵ�ǰ״̬�Ƿ���Ե���
        if schedulable(scheduleplan,H_DS, DSFR, RT, UD)
            %���ݵ��ȣ������������н�
            [ Flag,x,step,FP,FPPackage,scheduleplan ] = schedule( x,1,TKS,FP,scheduleplan,DSFR,PIPE1FR,PIPE2FR,FPORDER,FPPackage,RT,H_DS );
            if ~Flag
                pop(p,1:w) = rand(1,length(x));                                            %���³�ʼ������
                continue;
            else                
                %% �������
                scheduleplans{p} = scheduleplan;
                scheduleplan=sortrows(scheduleplan,1);                
               	%������Ӧ�Ⱥ���
               	f1 = gdmix(scheduleplan, c1);                                                 %�ܵ���ϳɱ�
              	f2 = gdimix(scheduleplan, c2);                                                %�޵׻�ϳɱ�
              	f3 = gchange(scheduleplan);                                                   %�������Ĺ��͹��л�����
               	f4 = gnum(scheduleplan);                                                       %���͹޸���
             	f5 = roundn(energecost(scheduleplan, c3),0);                           %�ܺĳɱ�
                %����
                eff(p,1:w) = x(1:w);   
                eff(p,w+1) = f1;
                eff(p,w+2) = f2;
                eff(p,w+3) = f3;
                eff(p,w+4) = f4;
                eff(p,w+5) = f5;
                %���ż�����һ������
                p = p + 1;
            end
        else
            disp('ϵͳ���ɵ��ȣ�\n');
            break;
        end
    end
end

%% ���ݵ���(ֻ��Ҫ�ҵ�һ�����н⼴��)
function [ Flag,X,step,FP,FPPackage,SchduleTable ] = schedule( X,step,TKS,FP,SchduleTable,DSFR,PIPE1FR,PIPE2FR,FPORDER,FPPackage,RT,H_DS )
% ֻ��Ҫ���α�����ȥ�������ֲ��ɵ��ȵĵط��ͻ���
% Ȼ����ű�����ֱ��FPΪ�գ������ҵ�һ�����н�
    Flag = false;%�Ƿ��ҵ����н�ı�־
    %�����ֵ
    Old_step = step;
    Old_FP = FP;
    Old_FPPackage = FPPackage;
    Old_SchduleTable = SchduleTable;
    times = 100;                                                                                   %�������ִ�д���
    global iteration;                                                                               %��¼����ִ�д���
   	iteration = iteration + 1;
    if iteration > times
        return;
    end
    
    %��������
    if round(sum(FP)) == 0
        Flag = true;
    elseif step <= 50   %�������ó���50
    	%����ET��UD(���ֻ��Ҫ����ͻ����)
        T = max(SchduleTable(SchduleTable(:,1) > 4,4));                     %�������
        if isempty(T)
            T = 0;
        end
      	[ET,UD] = get_et_ud(TKS,SchduleTable,FPORDER,FPPackage,T);
        %��¼�������㼣
        footprint = zeros(length(ET) + 1, length(UD));
        
    	%������ȱ���
        while ~all(all(footprint)) && ~Flag && iteration <= times
           	TK_NO=getInt(X(3*step-2), size(ET, 2));                                                                                         %���͹����
            DS_NO = getInt(X(3*step-1), size(UD, 2)-1)+1;                                                                               %������
          	DS = UD(DS_NO);
           	R = getInt(X(3*step), 3-1)+1;                                                                                                         %�ܵ�ת������
            COT = FPORDER(DS,find(FPPackage(DS,:), 1 ));                                                                         %ȷ��ԭ������   
          	PET1 = max(SchduleTable(SchduleTable(:,1) == 5, 4));                                                                 %����ܵ�1����ת��ʱ��
            if isempty(PET1)
                PET1 = 0;
            end
         	PET2 = max(SchduleTable(SchduleTable(:,1) == 6, 4));                                                                 %����ܵ�2����ת��ʱ��
          	DSFET = max(SchduleTable(SchduleTable(:,1) == DS, 4));                                                            %����DS�����ͽ���ʱ��
            
           	f1 = false;                                                                                                                                           %��ǰת�˵İ�ȫ��������ʱ��������ת��
            
            %Ϊ�˱�����ݹ�����ɵļ�����۹��ߣ���ˣ������ȵ��ȸ��۵�ԭ���ٵ��ȵ��۵�ԭ�͵Ĳ���
            if ismember(H_DS,UD) && DS ~= H_DS
                f1 = false;
            elseif DS == H_DS
                %���۵�ת�˼�������
                if TK_NO ~= 0                                                                                                                                %�����͹�ѡ��Ϊ0��������Ҫ����ת�˲��������򣬲������κβ���
                    TK = ET(TK_NO) ;                                                                                                                      %ȷ�����͹�
                    %���㵱ǰPET1ʱ��ԭ�Ϳ��������
                    total = gettotal(DSFR,SchduleTable,PET2);
                    %��ȫ��ԭ�����
                    V = getVolume(DS,TK,UD,H_DS,RT,total,DSFR,PIPE1FR(R),PIPE2FR(2),SchduleTable);
                    
                    %��������Ե����ض�����ֵ
                    if V >= 15000
                        if TKS(TK, 1) < V
                            V = TKS(TK, 1);
                        end
                        if FPPackage(DS,find(FPPackage(DS,:), 1 )) < V
                            V = FPPackage(DS,find(FPPackage(DS,:), 1 ));
                        end

                        %���ϰ�
                        FP(COT) = FP(COT) - V;
                        FPPackage(DS,find(FPPackage(DS,:), 1 )) = FPPackage(DS,find(FPPackage(DS,:), 1 )) - V;
                     	SchduleTable = [SchduleTable; 6, TK, PET2, PET2+V/PIPE2FR(2), COT,R,V];                               %�ܵ�2ת�˼�¼(�����ת��)
                       	SchduleTable = [SchduleTable; DS, TK, DSFET, DSFET + V/DSFR(DS), COT,2,V];                        %�ܵ�2���ͼ�¼
                        f1 = true;
                    end
                end
            else
                if TK_NO==0
                    PipeStoptime  = roundn(getPipeStoptime(DSFR, RT, H_DS,SchduleTable,TKS,PET1),-1);%�������룬ȡ6λС��
                    if PipeStoptime > 0
                        %ͣ��(������ͣ��)
                        SchduleTable = stop(SchduleTable,PipeStoptime);
                        f1 = true;
                    end
                else
                    TK = ET(TK_NO) ;                                                           %ȷ�����͹�
                    %���㵱ǰPET1ʱ��ԭ�Ϳ��������
                    total = gettotal(DSFR,SchduleTable,PET1);
                    %��ȫ��ԭ�����
                    V = getVolume(DS,TK,UD,H_DS,RT,total,DSFR,PIPE1FR(R),PIPE2FR(2),SchduleTable);
                    if V >= 15000
                        if TKS(TK, 1) < V
                            V = TKS(TK, 1);
                        end
                        if FPPackage(DS,find(FPPackage(DS,:), 1 )) < V
                            V = FPPackage(DS,find(FPPackage(DS,:), 1 ));
                        end

                        %���ϰ�
                        FP(COT) = FP(COT) - V;
                        FPPackage(DS,find(FPPackage(DS,:), 1 )) = FPPackage(DS,find(FPPackage(DS,:), 1 )) - V;
                       	%����ת�˼�¼�����ͼ�¼
                     	SchduleTable = [SchduleTable;  length(DSFR)+1, TK, PET1, PET1 + V / PIPE1FR(R), COT,R,V];     %ֻ��¼ʹ�õ��������
                       	SchduleTable = [SchduleTable;  DS, TK, DSFET, DSFET + V / DSFR(DS), COT,R,V];
                       	f1 = true;
                    end
                end
            end
            
                    
         	%�ж��Ե����Ƿ�ɹ�����һ״̬�Ƿ���Ե���
          	if f1 && conflic_detect(SchduleTable) && schedulable(SchduleTable, H_DS, DSFR, RT, UD)
                f2 = true;
            else
                f2 = false;
            end
            
          	%�ɵ��ȣ�������һ��
          	if f2
              	step = step + 1;
                [ Flag,X,step,FP,FPPackage,SchduleTable ] = schedule( X,step,TKS,FP,SchduleTable,DSFR,PIPE1FR,PIPE2FR,FPORDER,FPPackage,RT,H_DS );
            end
                        
            %��δ�ҵ����н⣬��ǲ��ع�
            if ~Flag && iteration <= times
%                 clc;
%                 step
%                 FPPackage

                footprint(TK_NO + 1,DS_NO) = 1;                                                 %�ȳ��Բ�ͬ�������ٳ��Բ�ͬ�Ĺ�
                %���ݻع�
                [step,FP,FPPackage,SchduleTable] = rollback(Old_step,Old_FP,Old_FPPackage,Old_SchduleTable);
            
               	%���Ĳ��ԣ��������Ƿ����ͣ�ˣ�
                f_j = false;
                if ~all(all(footprint))                                                                         %������δ���Թ���·��
                    for i = 1:size(footprint,1)
                        if i == 1 && any(footprint(i,:))                                                  %���Թ�ͣ�˺��ټ�������
                            footprint(i,:) = 1;
                        end
                        if ~all(footprint(i,:))                                                                 %��ǰ���͹޻�����δ���Թ���������
                            while TK_NO + 1 ~= i
                                X(3*step-2) = rand;
                                TK_NO = getInt(X(3*step-2),size(ET, 2));
                            end
                            for j = 1:size(footprint,2)
                                if footprint(i,j) == 0                                                              %δ���Թ���·��
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

%�ж��Ƿ����ì�ܵĵ���(flag==true��������ͻ)
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
        t = [];                                                                                 %����
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
        
        %��������Ƿ��ͻ
        if size(t,1) > 1
            for m=1:size(t,1) - 1
                for n=m+1:size(t,1)
                    t1 = t(m,:);
                    t2 = t(n,:);
                    if max(t1(1),t2(1)) < min(t1(2),t2(2))
                        flag = false;           %��鵽��ͻ���˳�
                        return;
                    end
                end
            end
        end
    end
end

%% �жϵ�ǰϵͳ״̬�Ƿ�ɵ���(���������۵������������Ƿ��������ޣ�ͬʱ�ж��Ƿ��ͻ)
function [ flag ] = schedulable(SchduleTable, H_DS, DSFR, RT, UD)
    flag = true;
    %��������۵����͵��۵���
    if isempty(UD)        
        disp('udΪ��');
    end
    LT_UD = UD(UD~=H_DS);
    HT_UD =H_DS;
    
    %�жϵ��۵����Ƿ����������������
    if ~isempty(LT_UD)
        PET1 = max(SchduleTable(SchduleTable(:,1) == 5,4));                                   %�ܵ�1��ת�˽���ʱ��
        if isempty(PET1)
            PET1 = 0;
        end
        total = gettotal(DSFR,SchduleTable,PET1);                                                   %���������
        
        for i=1:length(LT_UD)
            DSN = LT_UD(i);
            arfai = RT * DSFR(DSN);

            if round(total(DSN)) < length(LT_UD) * arfai                                                %�ж��Ƿ������������Լ��
                flag=false;
                %disp('���۵�ԭ���������㣡');
                return;
            end
        end
    end
    
    %�жϸ��۵����Ƿ����������������(�ܹ�ʹ��������ͺ����������ȼ���)
    if ~isempty(HT_UD)
        PET2 = max(SchduleTable(SchduleTable(:,1) == 6,4));                                   %�ܵ�2��ת�˽���ʱ��
        if isempty(PET2)
            PET2 = 0;
        end
        total = gettotal(DSFR,SchduleTable,PET2);                                                   %���������
    
        if length(HT_UD) > 1                                                                                    %���۵�������1���������ʱ������
            flag = false;
            disp('���۵�������1���������ʱ�����ǣ�');
            return;
        end
        if total(HT_UD) < RT * DSFR(HT_UD);                %����ת��
            %disp('���۵�ԭ����������');
            flag = false;
            return;
        end
    end
end

%% ��ȫ����
function [V] = getVolume(DS,TK,UD,H_DS,RT,total,DSFR,PIPE1FR,PIPE2FR,SchduleTable)
    %���۵�ԭ�ͺ͵��۵�ԭ�Ͷ�Ҫ���ǣ����У�PIPE1FR,PIPE2FRΪһ�����֣���������
    if isempty(UD)        
        disp('udΪ��');
    end
    LT_UD = UD(UD~=H_DS);
    HT_UD =H_DS;
    
    if DS ~= H_DS
        if ~isempty(LT_UD)
            %�������㹩�͹�ռ��ʱ�䲻��ͻ��������
            Tmp_SchduleTable = SchduleTable(SchduleTable(:,2) == TK,:);
            if size(Tmp_SchduleTable,1) > 0
                s = 1;
                k = 1;
                interval = [];                                                                                 %����
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
                PET1 = max(SchduleTable(SchduleTable(:,1) == 5,4));                                   %�ܵ�1��ת�˽���ʱ��
                if isempty(PET1)
                    PET1 = 0;
                end

                for i=size(interval,1)
                    if interval(i,1) <= PET1 && interval(i,2) > PET1                                            %����Ƿ��ͻ����ȫ���Ϊ0
                        V = 0;
                        return;
                    end
                end
            end
            
            %�жϵ��۵����Ƿ����������������
            for i=1:length(LT_UD)
                DSL = LT_UD(i);                                                                                             %���۵�������
                arfai = RT * DSFR(DSL);                                                                                 %���۵���������ԭ���������
              	Vr = round(total(DSL)) - length(LT_UD) * arfai;
                t = Vr/DSFR(DSL);
            end
            
            V = round(t * PIPE1FR);
        end
    else
        %�жϸ��۵����Ƿ����������������(�ܹ�ʹ��������ͺ����������ȼ���)
        if ~isempty(HT_UD)
            if length(HT_UD) > 1                                                                                    %���۵�������1���������ʱ������
                disp('���۵�������1���������ʱ�����ǣ�');
                return;
            end
            
            %�������㹩�͹�ռ��ʱ�䲻��ͻ��������
        	Tmp_SchduleTable = SchduleTable(SchduleTable(:,2) == TK,:);
            if size(Tmp_SchduleTable,1) > 0
                s = 1;
                k = 1;
                interval = [];                                                                                 %����
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
                PET2 = max(SchduleTable(SchduleTable(:,1) == 6,4));                                   %�ܵ�2��ת�˽���ʱ��
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
    
            %������ԭ���������
          	DSH = HT_UD(1);
           	Vr = total(DSH) - RT * DSFR(DSH);                %����ת��
         	t = Vr/DSFR(DSH);
                        
            %�������ת�����
            V = round(t * PIPE2FR);
        end
    end
end

%% ͣ��
function [SchduleTable] = stop(SchduleTable,PipeStoptime)
	%ͣ���ڼ䣬�����������Ƿ����ͽ���
	PETOLD = max(SchduleTable(SchduleTable(:,1) == 5, 4));                                                         %������۵�ת�˽���ʱ��
    if isempty(PETOLD)
        PETOLD = 0;
    end
    tmp = ones(1,size(SchduleTable,1))*inf;
    for i=1:size(SchduleTable,1)
        if SchduleTable(i,1) <= 4 && SchduleTable(i,4) > PETOLD && SchduleTable(i,4) <= PETOLD + PipeStoptime               %ͣ��
            tmp(i) = SchduleTable(i,4);
        end
    end
    [m,n] = min(tmp);
    if m ~= inf
        tk = SchduleTable(n,2);
    else
        tk = 0;
    end
    
	%���͹��ͷ�
	if tk ~= 0
        %����ת�˽���ʱ��Ϊ���͹޿��õ�ʱ��
        PET = m;
    else
        %����ͣ��
        PET = PETOLD + PipeStoptime;
	end
	%����ת�˼�¼������Ҫ��ѡ���͹޹���
	SchduleTable = [SchduleTable; 5, 0, PETOLD, PET, 0,0,0];                        %����5����ܵ�1ͣ��
end

%% ���ݵ��ȼ�¼����ȡET��UD
function [ET,UD] = get_et_ud(TKS,scheduletable,FPORDER,FPPackage,t)
    ET = [];
    UD = [];
    %���ҿ��ù�
    for i=1:size(TKS,1)
        tmp = scheduletable(scheduletable(:,2)==i,:);
        
        %��δʹ�õĹ��͹�
        if isempty(tmp)
            ET = [ET,i];
            continue;
        end
        
        %ԭ�ͳ�ʼ��滹û����
        if tmp(1,6) == 0 && tmp(1,4) > t
            continue;
        end
        
        %�ж�ʱ��t�Ƿ���ʹ��
        flag = true;
        for j=1:size(tmp,1)
            %ת�˼�¼���������ͼ�¼��������֮��Ҳ������ʹ�á�
            if j + 1 <= size(tmp,1)
                if tmp(j,1) > 4 && tmp(j+1) <= 4            %����4�����������ĸ���
                    if tmp(j,3) <= t && tmp(j+1,4) >= t
                        flag = false;
                        break;
                    end
                end
            end
            
            %������ת�˼�¼���߹��ͼ�¼����ֻ��Ҫ���t�Ƿ����俪ʼʱ��ͽ���ʱ�伴�ɡ�
            if tmp(j,3) <= t && tmp(j,4) >= t
                flag = false;
                break;
            end
        end
        
        if flag
            ET = [ET,i];
        end
    end
    
    %���δ������������������
    for i=1:size(FPORDER,1)
        if round(sum(FPPackage(i,:))) ~= 0
            UD = [UD,i];            
        end
    end
end

%% �����������д��ڵ�ԭ��
function [total,scheduleplan,FPPackage,PET1,PET2] = schedule_oil(FPORDER,FPPackage,TKS,DSFR,TKHOT,VPIPE,PIPE2FR)
    scheduleplan = [];                                        %���ȼƻ�
    DSFET = zeros(1,size(FPORDER,1));           %���������һ�ι��ͽ���ʱ��
    total = zeros(1,size(FPORDER,1));               %�������ܵ�����
    PET1 = 0;
    PET2 = 0;
    
    %�����������ԭ�����������һ�ι��ͽ���ʱ�䣬�Լ���ʼ���ȼƻ�
    for i=1:size(TKS,1)  
        if TKS(i,2)~=inf                                                              %���˳��չ�
            if i == TKHOT                                                            %ʢ�ż��ȹܵ�����Ҫ��ԭ�͵Ĺ��͹�
                PET2 = TKS(i,3) / PIPE2FR(1);                               %���ټ��ȹܵ�2
                scheduleplan = [scheduleplan; 6, i, 0, PET2, TKS(i,2),0,TKS(i,3)];                                          %��¼�ܵ�2���ȹܵ���¼(����)
                scheduleplan = [scheduleplan; 6, i, PET2, PET2 + VPIPE/PIPE2FR(1), TKS(i,2),0,VPIPE];     %��¼�ܵ�2�����۵�ԭ�ͼ����ܵ���¼����
                PET2 = PET2 + VPIPE/PIPE2FR(1);                       %�������۵�ԭ��
                TKS(i,3) = VPIPE;                                                  %��������
            end
           	OT = TKS(i,2);                                    %ԭ������
           	[m,n] = find(FPORDER == OT);

          	%ѡ������������
          	if length(m) == 1 && length(n) == 1
            	ds = m;
            else
                [~,ind] = min(DSFET(m));                %Ϊ����Ҫ�͵������ͣ�̰�Ĺ�����������ͬʱ��Ҫһ����ʱ������������Ҫ����������
               	ds = m(ind);
            end

            %��ȡds��������˳�򣬲�����FPPackage
            DS_FPORDER = FPORDER(ds,:);
            for j=1:length(DS_FPORDER)
                if OT == DS_FPORDER(j)
                    FPPackage(ds,j) = FPPackage(ds,j) - TKS(i,3);
                    break;
                end
            end
            
         	total(ds) = total(ds) + TKS(i,3);           %���������������������
           	Temp = DSFET(ds);                           %����ǰһ�εĹ��ͽ���ʱ��
          	%���㹩�͹޵Ĺ��ͽ���ʱ�䣬����ӵ����͹޶�Ӧ�Ĺ��ͼ�¼����
           	DSFET(ds) = DSFET(ds) + TKS(i,3) / DSFR(ds);      %DSFET(i)����������i�����һ�ι��ͽ���ʱ�䣬TKS(k,3)�����͹�k�ĵ�ǰ������DSFR(i)����������i����������
          	scheduleplan = [scheduleplan; ds, i, Temp, DSFET(ds), TKS(i,2),0,TKS(i,3)];     %��¼��ʼ���ȼƻ�
        end
   	end
end

%% ��������
function [RT,DSFR,PIPE1FR,PIPE2FR,FPORDER,FPPackage,c1,c2,c3,TKS,FP,VPIPE,TKHOT,H_DS] = parameter_setting()
 	%% ��������(�������)
   	RT = 4;                                                                          %פ��ʱ��
	DSFR = [250,334,250,625];                                             %��������������(�ɼ�������������)
    PIPE1FR = [833.3 1250 1375];                                        %�ܵ�1����������
    PIPE2FR = [625,420];                                                     %�ܵ�2����������
    FPORDER = [3,4,8,8;1,2,11,11;5,4,6,10;7,8,9,9];              %������������˳��
    Plan = [25000,24000,44055,0;16000,26000,79285,0;15000,17000,18000,43055;22000,103000,107638,0];            %�������������
    VPIPE = 18000;                                                              %�ܵ������
    TKHOT = 8;                                                                    %װ�м���ԭ�͵Ĺ��͹ޣ�����ù��͹��е�ԭ������Ѿ��㹻���ȹܵ���
    H_DS = 2;                                                                      %��ʶ���۵�����������2����Ϊ���۵���
    
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
        10	12	31	28	25	27	32	35	21	24	0];                                                    %�ܵ���ϳɱ�
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
        10	12	31	28	25	27	32	35	21	24	0];                                                 %�޵׻�ϳɱ�
	c3 = [1 2 3];                                                                   %�ܺĳɱ�
    
	%% ��������(�ɱ����)
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
            ];                                        %�͹޳�ʼ��棬������Զ�����������͹��Ϳ�ʼʱ��ͽ���ʱ��
	EO = TKS(TKS(:,2)~=inf,:);             %���еķǿ��͹�
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
	FP = Total - TKO;               %���ϰ�(�ɼ���ԭ������)
    FP(FP<0)=0;                      %ȥ������ж�����ͣ�����ֵ��
    FPPackage = Plan;            %��ʱ��FPPackage==Plan
end

%% ���ݻع�  [step,FP,SchduleTable] = rollback(Old_step,Old_FP,Old_SchduleTable);
function [step,FP,FPPackage,SchduleTable] = rollback(Old_step,Old_FP,Old_FPPackage,Old_SchduleTable)
        step = Old_step;
        FP = Old_FP;
        FPPackage = Old_FPPackage;
        SchduleTable = Old_SchduleTable;
end

%% �������ԭ������
function [total] = gettotal(DSFR,ScheduleTable,T)
   	total = zeros(1,length(DSFR));
    %���ݵ��ȼƻ����������п��
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

%% ���͹��л�����
function [ x ] = gchange( a )
    K = 0;
    
    %��ȡת�˼�¼
    for i=1:size(a,1)
        if a(i,1) >= 5
            K = i;
            break;
        end
    end
    x = K - 1;
end

%% �޵׻�ϴ���
function [ x ] = gdimix( a, c2 )
%c2���޵׻�ϳɱ�
    m2 = zeros(size(c2));
    K = 0;
    
    %��ȡת�˼�¼
    for i=1:size(a,1)
        if a(i,1) >= 5
            K = i;
            break;
        end
    end
    K = K - 1;
    record = a(1:K,:);
    record = sortrows(record,2);
    
    %ͳ�Ƹ����͹޵�ʹ�ô���
    A = record(:,2)';
    B = unique(A);
    tks = [];
    for i=1:size(B,2)
        tks(i) = sum(A(:)==B(i));
    end
    
	%�������ͼ�¼������޵׻�ϴ���
    for i=1:size(tks,2)     %i�����͹�
        tmp = 0;
        for j=1:(size(record,1) - 1)      %j�������ͼ�¼
            if record(j,2) == i
                tmp =  tmp + 1;
                if tmp >= tks(i)     %tks(i)���͹�i�е�ʹ�ô���
                    break;
                end
                if record(j,5) ~= record(j+1,5)
                    m2(record(j,5), record(j+1,5)) = m2(record(j,5), record(j+1,5)) + 1;
                end
            end
        end
    end
    x = sum(sum(m2.*c2));          %�޵׻�ϳɱ�
end

%% ��ȡ�ܵ���ϳɱ�
function [ x ] = gdmix( a, c1 )
%c1���ܵ���ϳɱ�
    m1 = zeros(size(c1));
    K = 0;
    
    %��ȡת�˼�¼
    for i=1:size(a,1)
        if a(i,1) >= 5
            K = i;
            break;
        end
    end
    PIPE = a(K:size(a,1),:);
    PIPE = PIPE(PIPE(:,5)~=0,:);   %ɾ�����е�ͣ�˼�¼
    PIPE1 = PIPE(PIPE(:,1)==5,:);   %�ܵ�1��ͣ�˼�¼
    PIPE2 = PIPE(PIPE(:,1)==6,:);   %�ܵ�2��ͣ�˼�¼
    %�ֱ����ܵ�1�͹ܵ�2��ԭ�ͻ�ϳɱ�
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
    x = sum(sum(m1.*c1));       %�ܵ���ϴ���
end

%% ���͹޸���
function [ x ] = gnum( a )
    TKSet = [];
    for i=1:size(a,1)
        if a(i,2) ~= 0 && ~ismember(a(i,2), TKSet)
            TKSet = [TKSet, a(i,2)];
        end
    end
    x = size(TKSet,2);      %���͹޵ĸ���
end

%% �ܺĳɱ�
function [ x ] = energecost( a,c )
    pipplan = a(a(:,1) >= 5 & a(:,6) ~= 0,:);
    x = c(pipplan(:,6))*(pipplan(:,4)-pipplan(:,3));       %ʱ��*��λʱ���ܺĳɱ�
end

%% ����[0~1]֮��������������һ��[0~length]����
function [ r ] = getInt( a, length )
    %a��[0~1]
    if a == 1
        r = length;
    else
        r = fix((length + 1) * a);
    end
end

%% ��ȡ�ܵ�ͣ�˰�ȫʱ��
function [ t ] = getPipeStoptime(DSFR, RT, H_DS,ScheduleTable,TKS,T)
    t=inf;
    [total] = gettotal(DSFR,ScheduleTable,T);
    
    %����Լ��(���۵�ԭ�ͺ͵��۵�ԭ�ͷֿ�����)
    DSFR_L = DSFR;
    DSFR_L(H_DS) = [];                                                                %���۵�����������������
    total(H_DS) = [];
    
    %����Ƿ���ͣ�˵ı�Ҫ
    if length(unique(ScheduleTable(:,2))) < size(TKS,1)                     %����û��ʹ�õĹ��͹ޣ���ˣ�û��ͣ�˵ı�Ҫ
        t = 0;
        return;
    else
    	TKs = unique(ScheduleTable(:,2));
        Temp_ScheduleTable = ScheduleTable(ScheduleTable(:,1)<=4,:);
        for i=1:length(TKs)
            TK = TKs(i);
            if max(Temp_ScheduleTable(Temp_ScheduleTable(:,2)==TK,4)) < T            %���ڹ��͹ޣ�������Ҫת��
                t = 0;
                return;
            end
        end
    end
    
    %���۵�ԭ������Լ��
    K = length(DSFR_L);
    for i=1:K
        arfai = RT * DSFR_L(i);

        %Ԥ��һ���������ڵ�ԭ��
      	total(i) = total(i) - K * arfai;

        tmp = total(i) / DSFR_L(i);
        if tmp < t
            t = tmp;
        end
    end
end