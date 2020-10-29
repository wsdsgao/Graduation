function out = GMSK_demod(wav, c0, c1, oversamp, indi, iter)
% input: 
%        wav        待解调GMSK基带信号波形
%        c0         匹配滤波器1
%        c1         匹配滤波器2
%        oversamp   过采样倍数
%        indi       =0，最终维特比译码长度为iter  
%                   =1，最终维特比译码长度为iter-1
%        iter       维特比译码长度       
%
% output：
%        out   解调结果

    Nc0=length(c0);  % 匹配滤波器1长度
    Nc1=length(c1);  % 匹配滤波器2长度

    r0=conv(c0, wav);  % 将待解调的波形通过匹配滤波器1
    r0n=r0(Nc0:oversamp:end);  % 8倍抽取
    r1=conv(c1, wav);  % 将待解调的波形通过匹配滤波器2
    r1n=r1(Nc1:oversamp:end);  % 8倍抽取
    length1=1e03;  % 路径1 度量值
    length2=1e03;  % 路径2 度量值
    length3=1e03;  % 路径3 度量值
    length4=1e03;  % 路径4 度量值
    l1=[1 1];   % 路径1 路径
    l2=[-1 1];  % 路径2 路径
    l3=[1 -1];  % 路径3 路径
    l4=[-1 -1]; % 路径4 路径
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SOVA = [0 0];  % 每一个soft output值为log(P(a=0|r)/P(a=1|r)) 
                   % >0 表示为0概率大
                   % <0 表示为1概率大

    for n=1:(iter)/2-1
        
        % 奇数位
        g1=l1;  % 继承上一次计算留下的路径1的路径
        g2=l2;  % 继承上一次计算留下的路径2的路径
        g3=l3;  % 继承上一次计算留下的路径3的路径
        g4=l4;  % 继承上一次计算留下的路径4的路径
        length11=imag(r0n(2*n-1))-real(r1n(2*n-1))+length1;  %  1  1 ->  1
        length12=imag(r0n(2*n-1))+real(r1n(2*n-1))+length2;  % -1  1 ->  1
        length21=imag(r0n(2*n-1))+real(r1n(2*n-1))+length3;  %  1 -1 ->  1
        length22=imag(r0n(2*n-1))-real(r1n(2*n-1))+length4;  % -1 -1 ->  1
        length31=-imag(r0n(2*n-1))+real(r1n(2*n-1))+length1; %  1  1 -> -1
        length32=-imag(r0n(2*n-1))-real(r1n(2*n-1))+length2; % -1  1 -> -1
        length41=-imag(r0n(2*n-1))-real(r1n(2*n-1))+length3; %  1 -1 -> -1
        length42=-imag(r0n(2*n-1))+real(r1n(2*n-1))+length4; % -1 -1 -> -1

        % 根据计算所得路径度量值进行取舍
        if length11>length12
            length1=length11;
            l1=[g1,1];
        else
            length1=length12;
            l1=[g2,1];
        end

        if length21>length22
            length2=length21;
            l2=[g3,1];
        else
            length2=length22;
            l2=[g4,1];
        end

        if length31>length32
            length3=length31;
            l3=[g1,-1];
        else
            length3=length32;
            l3=[g2,-1];
        end

        if length41>length42
            length4=length41;
            l4=[g3,-1];
        else
            length4=length42;
            l4=[g4,-1];
        end
        
        length_min = min(min(length1, length2), min(length3, length4))-1e-03;
        SOVA = [SOVA log((length3+length4-length_min*2)/(length1+length2-length_min*2))];
        

        % 通过indi判断 最后1次偶数位的计算是否需要跳过
        if(n == (iter)/2-1)
            if(indi)
                continue;
            end
        end
        
        % 偶数位
        g1=l1;  % 继承上一次计算留下的路径1的路径
        g2=l2;  % 继承上一次计算留下的路径2的路径
        g3=l3;  % 继承上一次计算留下的路径3的路径
        g4=l4;  % 继承上一次计算留下的路径4的路径
        length11=real(r0n(2*n))-imag(r1n(2*n))+length1;  %  1  1 ->  1
        length12=real(r0n(2*n))+imag(r1n(2*n))+length2;  % -1  1 ->  1
        length21=real(r0n(2*n))+imag(r1n(2*n))+length3;  %  1 -1 ->  1
        length22=real(r0n(2*n))-imag(r1n(2*n))+length4;  % -1 -1 ->  1
        length31=-real(r0n(2*n))+imag(r1n(2*n))+length1; %  1  1 -> -1
        length32=-real(r0n(2*n))-imag(r1n(2*n))+length2; % -1  1 -> -1
        length41=-real(r0n(2*n))-imag(r1n(2*n))+length3; %  1 -1 -> -1
        length42=-real(r0n(2*n))+imag(r1n(2*n))+length4; % -1 -1 -> -1

        % 根据计算所得路径度量值进行取舍
        if length11>length12
            length1=length11;
            l1=[g1,1];
        else
            length1=length12;
            l1=[g2,1];
        end

        if length21>length22
            length2=length21;
            l2=[g3,1];
        else
            length2=length22;
            l2=[g4,1];
        end

        if length31>length32
            length3=length31;
            l3=[g1,-1];
        else
            length3=length32;
            l3=[g2,-1];
        end

        if length41>length42
            length4=length41;
            l4=[g3,-1];
        else
            length4=length42;
            l4=[g4,-1];
        end
        
        length_min = min(min(length1, length2), min(length3, length4));
        SOVA = [SOVA log((length3+length4-length_min*2)/(length1+length2-length_min*2))];
        
    end

    % 计算完全部波形，选择四条路径中度量值最大的一条
    path_length=max(max(length1,length2),max(length3,length4));
    if length1==path_length
        out=l1;
    elseif length2==path_length
        out=l2;
    elseif length3==path_length
        out=l3;
    elseif length4==path_length
        out=l4;
    end 

    a1 = out(2:1:end); 
    SOVA_1 = SOVA(2:1:end);
    a2 = out(1:1:end-1);
    SOVA_2 = SOVA(1:1:end-1);
    out = a1 .* a2;  % 将译码结果后一位与前一位相乘
    SOVA_out = SOVA_1 .* SOVA_2;
    out(1:2:end)=-out(1:2:end); % 相乘结果奇数位取反，即为最终解调结果; 
                                % 解调序列与波形对应的原始数据序列的对应关系为:                                 
                                % Viterbi(3～end) = bit(4~end)
    SOVA_out(1:2:end) = -SOVA_out(1:2:end);
    aaaaa=1;  % 无用