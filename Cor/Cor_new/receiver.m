function result = receiver(rx, fh_pat_lib, th_pat_lib, pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, frame_num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          参数定义
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode_sel = 0;  % 用于记录最终判断当前接收波形属于哪种速率模式

num_bits_pn = 24;  % 同步头S1\S2长度
num_bits_pn_2 = 21;  % 同步头S3\S4长度

bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间
fs_IF = 1024e6;  % 射频、中频信号采样速率
fs_BB = 64e6;  % 基带信号采样速率
oversamp_BB = T * fs_BB;  % 基带信号过采样倍数
oversamp_IF = T * fs_IF;  % 中频信号过采样倍数

num_bits_pulse = 304;  % 2Mbps A\500Kbps\250Kbps 一个脉冲包的长度
                       % 2Mbps B 每个脉冲去掉前后跳数据部分后的长度
frame_counter = 1;  % 记录已检测到多少帧数据（仿真用）
% flag_frame = 0;  % 解调完一帧，该标志位置1，等待一帧时间后再捕获下一帧（仿真用）
% time_counter_mode1 = 6720*oversamp_IF - 0.5*oversamp_IF;  % 2Mbps A模式 等待时长
% time_counter_mode2 = 7956*oversamp_IF - 0.5*oversamp_IF;  % 2Mbps B模式 等待时长
% time_counter_mode3 = 26880*oversamp_IF + 10*oversamp_IF - 0.5*oversamp_IF;  % 500Kbps模式 等待时长
% time_counter_mode4 = 53760*oversamp_IF - 0.5*oversamp_IF;  % 250Kbps模式 等待时长
counter = 0;

flag_Capture_C = 0;  % 指示是否捕获到一帧

% 各模式一帧总长度
time_frame_mode1 = (304*12+512*6+103*12);  %统一都算了跳时和固定的间隔
time_frame_mode2 = (304*12+512*6+103*12);
time_frame_mode3 = (304*12+512*6)*4+103*48;
time_frame_mode4 = (304*12+512*6)*8+103*96;

% 各模式一帧有效数据长度
length_frame_mode1 = 3072;
length_frame_mode2 = 6144;
length_frame_mode3 = 10272;
length_frame_mode4 = 20544;

% 各模式一帧脉冲数量
num_pulses_mode1 = 12;
num_pulses_mode2 = 12;
num_pulses_mode3 = 48;
num_pulses_mode4 = 96;

% 各模式跳频图案（只有一种图案，下同）
fh_pat_mode1 = fh_pat_lib(1:num_pulses_mode1);
fh_pat_mode2 = fh_pat_lib(1:num_pulses_mode2);
fh_pat_mode3 = fh_pat_lib(1:num_pulses_mode3);
fh_pat_mode4 = fh_pat_lib(1:num_pulses_mode4);

% 各模式跳时图案
th_pat_mode1 = th_pat_lib(1:num_pulses_mode1);
th_pat_mode2 = th_pat_lib(1:num_pulses_mode2);
th_pat_mode3 = th_pat_lib(1:num_pulses_mode3);
th_pat_mode4 = th_pat_lib(1:num_pulses_mode4);


% 载入已存数据
load('lib/f_trans.mat');  % 21个频点
% 接收端滤波器
load('lib/filter/BPF_CHAN1_UD.mat');  % 带通 通道1
load('lib/filter/BPF_CHAN2_UD.mat');  % 带通 通道2
load('lib/filter/BPF_CHAN3_UD.mat');  % 带通 通道3
load('lib/filter/BPF_CHAN4_UD.mat');  % 带通 通道4
load('lib/filter/BPF_CHAN5_UD.mat');  % 带通 通道5
load('lib/filter/BPF_CHAN12_2_UD.mat');  % 带通 通道1、2 第二级滤波
load('lib/filter/BPF_CHAN34_2_UD.mat');  % 带通 通道3、4 第二级滤波
load('lib/filter/BPF_CHAN5_2_UD.mat');  % 带通 通道5 第二级滤波
load('lib/filter/c0_f_128.mat');  % 解调匹配滤波器1
load('lib/filter/c1_f_128.mat');  % 解调匹配滤波器2
load('filter/LPF_1.mat'); % 采样率1024MHz，通带20MHz
load('filter/LPF_2.mat'); % 采样率256MHz，通带10MHz
S_lpf1 = 32;
S_lpf2 = 32;
S_bpf = 253;

%生成本地波形（仅mode4用于同步捕获和精细校准）
[wav_S1_mode1, wav_S2_mode1] = wav_gen(pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 1);

[wav_S1_mode2, wav_S2_mode2] = wav_gen(pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 2);
    
[wav_S1_S3_mode3, wav_S4_S2_mode3] = wav_gen(pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 3);

[wav_S1_S3_mode4, wav_S4_S2_mode4] = wav_gen(pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 4);

%考虑一下前后是否还需补0，以及要是一帧被连续捕获了好几次怎么办？！

for i = 20:5:200   % 等待一帧的时间长度
    
    %   % 解调完一帧，等待一帧时间长度再做捕获（仿真用）
    %   if (flag_frame)  
    %       if (counter < time_counter - 1)
    %           counter = counter + 1;
    %           continue;
    %       else
    %           counter = 0;
    %           flag_frame = 0;
    %       end
    %   end
    
    
     % 未捕获成功
     if(~flag_Capture_C)
             
     % for mode1 
     
        % 截取对应一帧长度的波形
        temp_rx_mode1 = rx(i:i+time_frame_mode1*oversamp_IF-1);
        
        % 根据跳频图案对应的频点找到对应通道的波形
        wav_temp_mode1 = zeros(num_pulses_mode1, num_bits_pulse * oversamp_IF);
        % 第1个跳时跳频图案
        for pulse_idx = 1:num_pulses_mode1
            f_idx = fh_pat_mode1(pulse_idx); % 当前脉冲对应频点
            
            % 计算当前脉冲在一帧中的起始位置
            if pulse_idx == 1
                pos_pulse_mode1 = th_pat_mode1(1)/2;
            else
                pos_pulse_mode1 = sum(th_pat_mode1(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103) + th_pat_mode1(pulse_idx)/2;
            end

            % 射频 -> 中频 （只取出数据位）
            if (f_idx >= 1) && (f_idx <= 5)
                wav_temp1_mode1 = downConv_IF(temp_rx_mode1(pos_pulse_mode1*oversamp_IF+1:pos_pulse_mode1*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % 通道1
            elseif (f_idx >= 6) && (f_idx <= 10)
                wav_temp1_mode1 = downConv_IF(temp_rx_mode1(pos_pulse_mode1*oversamp_IF+1:pos_pulse_mode1*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % 通道2
            elseif (f_idx >= 11) && (f_idx <= 14)
                wav_temp1_mode1 = downConv_IF(temp_rx_mode1(pos_pulse_mode1*oversamp_IF+1:pos_pulse_mode1*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % 通道3
            elseif (f_idx >= 15) && (f_idx <= 18)
                wav_temp1_mode1 = downConv_IF(temp_rx_mode1(pos_pulse_mode1*oversamp_IF+1:pos_pulse_mode1*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % 通道4 
            elseif (f_idx >= 19) && (f_idx <= 21)
                wav_temp1_mode1 = downConv_IF(temp_rx_mode1(pos_pulse_mode1*oversamp_IF+1:pos_pulse_mode1*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % 通道5 
            end

            wav_temp_mode1(pulse_idx,:) = wav_temp1_mode1;

        end


        % 根据跳频图案对应的频点将对应通道的波形滤波至零中频
        % wav_temp(采样率1024) -> rx_pulse_mat(采样率64)
        rx_pulse_mat_mode1 = downConv(wav_temp_mode1(:,:), num_pulses_mode1, fh_pat_mode1);  % 对应前10ms
        % rx_pulse_mat_2_mode1 = downConv(wav_temp_mode1(:,:,2), num_pulses_mode1, fh_pat_2_mode1);  % 对应后10ms

        % plot(real(rx_pulse_mat_mode1(1, :)));
        % value = sum(abs(rx_pulse_mat_mode1(1, :)))


        % % 预取前后24bit位置的波形（即同步头） 等待后续处理
        % 前10ms
        D_S1_mode1 = rx_pulse_mat_mode1(:,1:24*oversamp_BB);  %一行是一个脉冲
        D_S1_one_mode1 = D_S1_mode1(:,4:oversamp_BB:end);  %和数据速率相同，一个点对应一个bit
        D_S2_mode1 = rx_pulse_mat_mode1(:,1120+1:1120+24*oversamp_BB);  %2240=280*oversamp_BB
        D_S2_one_mode1 = D_S2_mode1(:,4:oversamp_BB:end);

        % % 后10ms
        % D_S1_2_mode1 = rx_pulse_mat_2_mode1(:,1:24*oversamp_BB);
        % D_S1_2_one_mode1 = D_S1_2_mode1(:,8:oversamp_BB:end);
        % D_S2_2_mode1 = rx_pulse_mat_2_mode1(:,2240+1:2240+24*oversamp_BB);
        % D_S2_2_one_mode1 = D_S2_2_mode1(:,8:oversamp_BB:end);

        % 同步头捕获过程                 
        % 1倍数据速率波形数据做捕获
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% 
        %
        %                           前10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%             

        %原算法
        rx_corr_S1_pat_mode1 = zeros(1,num_pulses_mode1);
        rx_corr_S2_pat_mode1 = zeros(1,num_pulses_mode1);
        for j = 1:num_pulses_mode1

            rx_wav_S1_pat_mode1 = D_S1_one_mode1(j,1:22) .* conj(wav_S1_mode1(j,1:22));
            rx_wav_S2_pat_mode1 = D_S2_one_mode1(j,1:22) .* conj(wav_S2_mode1(j,1:22));   

            rx_corr_S1_pat_mode1(j) = abs(sum(rx_wav_S1_pat_mode1));
            rx_corr_S2_pat_mode1(j) = abs(sum(rx_wav_S2_pat_mode1));

            rx_corr_pat_mode1(j) =  rx_corr_S1_pat_mode1(j) + rx_corr_S2_pat_mode1(j);

        end

        corr_value_mode1 = sum(rx_corr_pat_mode1);

        corr_value_mode1_plot1((i-20)/5+1) = corr_value_mode1;

        % 求前后同步头相关峰(改进算法)
        corr_value_mode1_plot2((i-20)/5+1) = corr(rx_pulse_mat_mode1, wav_S1_S3_mode4, wav_S4_S2_mode4, 1);
        % rx_corr_S1_pat_mode1 = zeros(1,num_pulses_mode1);
        % rx_corr_S2_pat_mode1 = zeros(1,num_pulses_mode1);
        % rx_corr_pat_mode1 = zeros(1,num_pulses_mode1);
        % for j = 1:num_pulses_mode1

        %     rx_wav_S1_pat_mode1_1 = D_S1_one_mode1(j,1:8) .* (wav_S1_S3_mode4(j,1:8));
        %     rx_wav_S1_pat_mode1_2 = D_S1_one_mode1(j,9:16) .* conj(wav_S1_S3_mode4(j,9:16));
        %     rx_wav_S1_pat_mode1_3 = D_S1_one_mode1(j,17:24) .* conj(wav_S1_S3_mode4(j,17:24)); %考虑一下是否为1-22比较合理
        %     rx_wav_S2_pat_mode1_1 = D_S2_one_mode1(j,1:8) .* (wav_S4_S2_mode4(j,22:29)); %考虑一下是否为24-45比较合理
        %     rx_wav_S2_pat_mode1_2 = D_S2_one_mode1(j,9:16) .* conj(wav_S4_S2_mode4(j,30:37));
        %     rx_wav_S2_pat_mode1_3 = D_S2_one_mode1(j,17:24) .* conj(wav_S4_S2_mode4(j,38:45));   

        %     rx_corr_S1_pat_mode1(j) = abs(sum(rx_wav_S1_pat_mode1_1)+sum(rx_wav_S1_pat_mode1_2)+sum(rx_wav_S1_pat_mode1_3));
        %     rx_corr_S2_pat_mode1(j) = abs(sum(rx_wav_S2_pat_mode1_1)+sum(rx_wav_S2_pat_mode1_2)+sum(rx_wav_S2_pat_mode1_3));

        %     rx_corr_pat_mode1(j) =  rx_corr_S1_pat_mode1(j) + rx_corr_S2_pat_mode1(j);

        % end

        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                           后10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % 求前后同步头相关峰
        % rx_corr_S1_pat_2_mode1 = zeros(1,num_pulses_mode1);  
        % rx_corr_S2_pat_2_mode1 = zeros(1,num_pulses_mode1);
        % for j = 1:num_pulses_mode1

        %     rx_wav_S1_pat_2_mode1 = D_S1_2_one_mode1(j,1:22) .* conj(wav_S1_2_mode1(j,1:22));
        %     rx_wav_S2_pat_2_mode1 = D_S2_2_one_mode1(j,1:22) .* conj(wav_S2_2_mode1(j,1:22));   

        %     rx_corr_S1_pat_2_mode1(j) = abs(sum(rx_wav_S1_pat_2_mode1));  %求一行的和
        %     rx_corr_S2_pat_2_mode1(j) = abs(sum(rx_wav_S2_pat_2_mode1));

        % end
        
        % 同步捕获判决条件
        % 前(或后)同步头相关峰值大于阈值
        if (corr_value_mode1 > 44*num_pulses_mode1) %判决条件要和信号能量去比，还没加
            mode_sel = 1;
            flag_Capture_C = 1;
            % if (sum(rx_corr_S1_pat_1_mode1) > sum(rx_corr_S2_pat_1_mode1))
            %     tag = 1;
            % else
            %     tag = 2;
            % end

        % elseif (sum(rx_corr_S1_pat_2_mode1)>20*num_pulses_mode1) || (sum(rx_corr_S2_pat_2_mode1)>20*num_pulses_mode1)
        %     mode_sel = 1;
        %     flag_Capture_C = 1;
        %     if (sum(rx_corr_S1_pat_2_mode1) > sum(rx_corr_S2_pat_2_mode1))
        %         tag = 3;
        %     else
        %         tag = 4;
        %     end
        end        

        
    % for mode2
    
        % 截取对应一帧长度的波形
        temp_rx_mode2 = rx(i:i+time_frame_mode2*oversamp_IF-1);
        % 根据跳频图案对应的频点找到对应通道的波形
        % 两个连续10ms的图案
        wav_temp_mode2 = zeros(num_pulses_mode2, num_bits_pulse * oversamp_IF);  %没包含跳时的数据
        % 第1个跳时跳频图案
        for pulse_idx = 1:num_pulses_mode2
            f_idx = fh_pat_mode2(pulse_idx);  % 当前脉冲对应的频点

            % 计算当前脉冲在一帧中的起始位置
            if pulse_idx == 1
                pos_pulse_mode2 = 0;
            else
                pos_pulse_mode2 = sum(th_pat_mode2(1:pulse_idx-1)) + 103 * (pulse_idx-1) + (pulse_idx-1)*num_bits_pulse;
            end  

            % 计算当前脉冲的长度
            num_bits_pulse_rx = th_pat_mode2(pulse_idx) + num_bits_pulse; %包含了跳时

            % 射频 -> 中频
            if (f_idx >= 1) && (f_idx <= 5)    
                wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % 通道1
            elseif (f_idx >= 6) && (f_idx <= 10)
                wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % 通道2
            elseif (f_idx >= 11) && (f_idx <= 14)
                wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % 通道3
            elseif (f_idx >= 15) && (f_idx <= 18)
                wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % 通道4 
            elseif (f_idx >= 19) && (f_idx <= 21)
                wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % 通道5 
            end

            wav_temp_mode2(pulse_idx,:) = wav_temp1_mode2((th_pat_mode2(pulse_idx)/2)*oversamp_IF+1:(th_pat_mode2(pulse_idx)/2+num_bits_pulse)*oversamp_IF);  %没包含跳时的数据
        end


        % % 第2个跳时跳频图案
        % for pulse_idx = 1:num_pulses_mode2
        %     f_idx = fh_pat_2_mode2(pulse_idx);  % 当前脉冲对应的频点

        %     % 计算当前脉冲在一帧中的起始位置
        %     if pulse_idx == 1
        %         pos_pulse_mode2 = 100;
        %     else
        %         pos_pulse_mode2 = sum(th_pat_2_mode2(1:pulse_idx-1)) + 3 * (pulse_idx-1) + (pulse_idx-1)*num_bits_pulse + 100 * pulse_idx;
        %     end
            
        %     % 计算当前脉冲长度
        %     num_bits_pulse_rx = th_pat_2_mode2(pulse_idx) + 3 + num_bits_pulse;

        %     % 射频 -> 中频
        %     if (f_idx >= 1) && (f_idx <= 5)    
        %         wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % 通道1
        %     elseif (f_idx >= 6) && (f_idx <= 10)
        %         wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % 通道2
        %     elseif (f_idx >= 11) && (f_idx <= 14)
        %         wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % 通道3
        %     elseif (f_idx >= 15) && (f_idx <= 18)
        %         wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % 通道4 
        %     elseif (f_idx >= 19) && (f_idx <= 21)
        %         wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % 通道5 
        %     end

        %     wav_temp_mode2(pulse_idx,:,2) = wav_temp1_mode2((th_pat_2_mode2(pulse_idx)/2+3)*oversamp_IF+1:(th_pat_2_mode2(pulse_idx)/2+3+num_bits_pulse)*oversamp_IF);
        % end


        % 根据跳频图案对应的频点将对应通道的波形滤波至零中频
        % wav_temp(采样率1024) -> rx_pulse_mat(采样率64)
        rx_pulse_mat_mode2 = downConv_2MLL(wav_temp_mode2(:,:), num_pulses_mode2, fh_pat_mode2, th_pat_mode2);  % 对应前10ms
        % rx_pulse_mat_2_mode2 = downConv_2MLL(wav_temp_mode2(:,:,2), num_pulses_mode2, fh_pat_2_mode2, th_pat_2_mode2);  % 对应后10ms


        % % 预取前后24bit位置的波形 等待后续处理
        % % 前10ms
        % D_S1_mode2 = rx_pulse_mat_mode2(:,1:num_bits_pn*oversamp_BB);
        % D_S1_one_mode2 = D_S1_mode2(:,4:oversamp_BB:end);
        % D_S2_mode2 = rx_pulse_mat_mode2(:,2240+1:2240+num_bits_pn*oversamp_BB);
        % D_S2_one_mode2 = D_S2_mode2(:,4:oversamp_BB:end);
        % % 后10ms
        % D_S1_2_mode2 = rx_pulse_mat_2_mode2(:,1:num_bits_pn*oversamp_BB);
        % D_S1_2_one_mode2 = D_S1_2_mode2(:,8:oversamp_BB:end);
        % D_S2_2_mode2 = rx_pulse_mat_2_mode2(:,2240+1:2240+num_bits_pn*oversamp_BB);
        % D_S2_2_one_mode2 = D_S2_2_mode2(:,8:oversamp_BB:end);


        % 同步头捕获过程                 
        % 1倍数据速率波形数据做捕获
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                   前10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % 求前后同步头相关峰(带频偏、相偏)
        corr_value_mode2 = corr(rx_pulse_mat_mode2, wav_S1_S3_mode4, wav_S4_S2_mode4, 2);
        % rx_corr_S1_pat_1_mode2 = zeros(1,num_pulses_mode2);
        % rx_corr_S2_pat_1_mode2 = zeros(1,num_pulses_mode2);
        % for j = 1:num_pulses_mode2

        %     rx_wav_S1_pat_1 = D_S1_1_one_mode2(j,1:22) .* conj(wav_S1_1_mode2(j,1:22));
        %     rx_wav_S2_pat_1 = D_S2_1_one_mode2(j,1:22) .* conj(wav_S2_1_mode2(j,1:22));   

        %     rx_corr_S1_pat_1_mode2(j) = abs(sum(rx_wav_S1_pat_1));
        %     rx_corr_S2_pat_1_mode2(j) = abs(sum(rx_wav_S2_pat_1));

        % end


        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                       后10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % % 求前后同步头相关峰(带频偏、相偏)
        % rx_corr_S1_pat_2_mode2 = zeros(1,num_pulses_mode2);  
        % rx_corr_S2_pat_2_mode2 = zeros(1,num_pulses_mode2);
        % for j = 1:num_pulses_mode2

        %     rx_wav_S1_pat_2 = D_S1_2_one_mode2(j,1:22) .* conj(wav_S1_2_mode2(j,1:22));
        %     rx_wav_S2_pat_2 = D_S2_2_one_mode2(j,1:22) .* conj(wav_S2_2_mode2(j,1:22));   

        %     rx_corr_S1_pat_2_mode2(j) = abs(sum(rx_wav_S1_pat_2));
        %     rx_corr_S2_pat_2_mode2(j) = abs(sum(rx_wav_S2_pat_2));

        % end


         % 同步捕获判决条件
         % 前(或后)同步头相关峰值大于阈值       
        if (corr_value_mode2 > 44*num_pulses_mode2) 
            mode_sel = 2;
            flag_Capture_C = 1; 
            % if (sum(rx_corr_S1_pat_1_mode2) > sum(rx_corr_S2_pat_1_mode2))
            %     tag = 1;
            % else
            %     tag = 2;
            % end
        %  elseif (sum(rx_corr_S1_pat_2_mode2)>20*num_pulses_mode2) || (sum(rx_corr_S2_pat_2_mode2)>20*num_pulses_mode2)
        %     mode_sel = 2;
        %     flag_Capture_C = 1;
        %     if (sum(rx_corr_S1_pat_2_mode2) > sum(rx_corr_S2_pat_2_mode2))
        %         tag = 3;
        %     else
        %         tag = 4;
        %     end
        end 

    % for mode3

        if (i+time_frame_mode3*oversamp_IF-1) > length(rx)
            continue;
        end
    
        % 截取对应48个脉冲长度的波形
        temp_rx_mode3 = rx(i:i+time_frame_mode3*oversamp_IF-1);
        
        % 根据跳频图案对应的频点找到对应通道的波形
        % 两个连续10ms的图案
        wav_temp_mode3 = zeros(num_pulses_mode3, num_bits_pulse * oversamp_IF);
        % 第1个跳时跳频图案
        for pulse_idx = 1:num_pulses_mode3
            f_idx = fh_pat_mode3(pulse_idx);  % 当前脉冲对应的频点
            
            % 计算当前脉冲在一帧中的起始位置
            if pulse_idx == 1
                pos_pulse_mode3 = th_pat_mode3(1)/2;
            else
                pos_pulse_mode3 = sum(th_pat_mode3(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103) + th_pat_mode3(pulse_idx)/2;
            end

            % 射频 -> 中频
            if (f_idx >= 1) && (f_idx <= 5)    
                wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % 通道1
            elseif (f_idx >= 6) && (f_idx <= 10)
                wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % 通道2
            elseif (f_idx >= 11) && (f_idx <= 14)
                wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % 通道3
            elseif (f_idx >= 15) && (f_idx <= 18)
                wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % 通道4 
            elseif (f_idx >= 19) && (f_idx <= 21)
                wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % 通道5 
            end

            wav_temp_mode3(pulse_idx,:) = wav_temp1_mode3;

        end


        % 第2个跳时跳频图案
        % for pulse_idx = 1:num_pulses_mode3
        %     f_idx = fh_pat_2_mode3(pulse_idx);  % 当前脉冲对应的频点
            
        %     % 计算当前脉冲在一帧中的起始位置
        %     if pulse_idx == 1
        %         pos_pulse_mode3 = th_pat_2_mode3(1)/2;
        %     else
        %         pos_pulse_mode3 = sum(th_pat_2_mode3(1:pulse_idx-1)) + (pulse_idx-1)*num_bits_pulse + th_pat_2_mode3(pulse_idx)/2;
        %     end

        %     % 射频 -> 中频
        %     if (f_idx >= 1) && (f_idx <= 5)    
        %         wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % 通道1
        %     elseif (f_idx >= 6) && (f_idx <= 10)
        %         wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % 通道2
        %     elseif (f_idx >= 11) && (f_idx <= 14)
        %         wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % 通道3
        %     elseif (f_idx >= 15) && (f_idx <= 18)
        %         wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % 通道4 
        %     elseif (f_idx >= 19) && (f_idx <= 21)
        %         wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % 通道5 
        %     end

        %     wav_temp_mode3(pulse_idx,:,2) = wav_temp1_mode3;

        % end



        % 根据跳频图案对应的频点将对应通道的波形滤波至零中频
        % wav_temp(采样率1024) -> rx_pulse_mat(采样率64)
        rx_pulse_mat_mode3 = downConv(wav_temp_mode3(:,:), num_pulses_mode3, fh_pat_mode3);  % 对应前10ms
        % rx_pulse_mat_2_mode3 = downConv_500K(wav_temp_mode3(:,:,2), num_pulses_mode3, fh_pat_2_mode3);  % 对应后10ms


        % % 预取前后24bit位置的波形 等待后续处理
        % % 前10ms
        % D_S1_mode3 = rx_pulse_mat_mode3(:,1:24*oversamp_BB);
        % D_S1_one_mode3 = D_S1_mode3(:,4:oversamp_BB:end);
        % D_S2_mode3 = rx_pulse_mat_mode3(:,280*oversamp_BB+1:304*oversamp_BB);
        % D_S2_one_mode3 = D_S2_mode3(:,4:oversamp_BB:end);
        % D_S3_mode3 = rx_pulse_mat_mode3(:,24*oversamp_BB+1:45*oversamp_BB);
        % D_S3_one_mode3 = D_S3_mode3(:,4:oversamp_BB:end);
        % D_S4_mode3 = rx_pulse_mat_mode3(:,259*oversamp_BB+1:280*oversamp_BB);
        % D_S4_one_mode3 = D_S4_mode3(:,4:oversamp_BB:end);

        % % 后10ms
        % D_S1_2_mode3 = rx_pulse_mat_2_mode3(:,1:24*oversamp_BB);
        % D_S1_2_one_mode3 = D_S1_2_mode3(:,8:oversamp_BB:end);
        % D_S2_2_mode3 = rx_pulse_mat_2_mode3(:,280*oversamp_BB+1:304*oversamp_BB);
        % D_S2_2_one_mode3 = D_S2_2_mode3(:,8:oversamp_BB:end);
        % D_S3_2_mode3 = rx_pulse_mat_2_mode3(:,24*oversamp_BB+1:45*oversamp_BB);
        % D_S3_2_one_mode3 = D_S3_2_mode3(:,8:oversamp_BB:end);
        % D_S4_2_mode3 = rx_pulse_mat_2_mode3(:,259*oversamp_BB+1:280*oversamp_BB);
        % D_S4_2_one_mode3 = D_S4_2_mode3(:,8:oversamp_BB:end);            


        % 同步头捕获过程    (S1_S3____S4_S2)             
        % 1倍数据速率波形数据做捕获
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                   前10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % 求前后同步头相关峰
        corr_value_mode3 = corr(rx_pulse_mat_mode3, wav_S1_S3_mode4, wav_S4_S2_mode4, 3);
        % rx_corr_S1_pat_1_mode3 = zeros(1,num_pulses_mode3);
        % rx_corr_S2_pat_1_mode3 = zeros(1,num_pulses_mode3);
        % rx_corr_S3_pat_1_mode3 = zeros(1,num_pulses_mode3);
        % rx_corr_S4_pat_1_mode3 = zeros(1,num_pulses_mode3);
        % for j = 1:num_pulses_mode3

        %     rx_wav_S1_pat_1 = D_S1_1_one_mode3(j,1:24) .* conj(wav_S1_S3_1_mode3(j,1:24));     % S1 24bit
        %     rx_wav_S2_pat_1 = D_S2_1_one_mode3(j,1:22) .* conj(wav_S4_S2_1_mode3(j,21+1:43));  % S2 22bit  
        %     rx_wav_S3_pat_1 = D_S3_1_one_mode3(j,1:19) .* conj(wav_S1_S3_1_mode3(j,24+1:43));  % S3 19bit
        %     rx_wav_S4_pat_1 = D_S4_1_one_mode3(j,1:21) .* conj(wav_S4_S2_1_mode3(j,1:21));     % S4 21bit

        %     rx_corr_S1_pat_1_mode3(j) = abs(sum(rx_wav_S1_pat_1));
        %     rx_corr_S2_pat_1_mode3(j) = abs(sum(rx_wav_S2_pat_1));
        %     rx_corr_S3_pat_1_mode3(j) = abs(sum(rx_wav_S3_pat_1));
        %     rx_corr_S4_pat_1_mode3(j) = abs(sum(rx_wav_S4_pat_1));

        % end

        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                       后10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % rx_corr_S1_pat_2_mode3 = zeros(1,num_pulses_mode3);  
        % rx_corr_S2_pat_2_mode3 = zeros(1,num_pulses_mode3);
        % rx_corr_S3_pat_2_mode3 = zeros(1,num_pulses_mode3);
        % rx_corr_S4_pat_2_mode3 = zeros(1,num_pulses_mode3);
        % for j = 1:num_pulses_mode3

        %     rx_wav_S1_pat_2 = D_S1_2_one_mode3(j,1:24) .* conj(wav_S1_S3_2_mode3(j,1:24));     % S1 24bit
        %     rx_wav_S2_pat_2 = D_S2_2_one_mode3(j,1:22) .* conj(wav_S4_S2_2_mode3(j,21+1:43));  % S2 22bit  
        %     rx_wav_S3_pat_2 = D_S3_2_one_mode3(j,1:19) .* conj(wav_S1_S3_2_mode3(j,24+1:43));  % S3 19bit
        %     rx_wav_S4_pat_2 = D_S4_2_one_mode3(j,1:21) .* conj(wav_S4_S2_2_mode3(j,1:21));     % S4 21bit
            
        %     rx_corr_S1_pat_2_mode3(j) = abs(sum(rx_wav_S1_pat_2));
        %     rx_corr_S2_pat_2_mode3(j) = abs(sum(rx_wav_S2_pat_2));
        %     rx_corr_S3_pat_2_mode3(j) = abs(sum(rx_wav_S3_pat_2));
        %     rx_corr_S4_pat_2_mode3(j) = abs(sum(rx_wav_S4_pat_2));

        % end


         % 同步捕获判决条件
         % 前(或后)同步头解调误码数的12个脉冲平均值不超过4        
         if (corr_value_mode3 > 82*num_pulses_mode3)
             flag_Capture_C = 1;
             mode_sel = 3;
             
        %      if (sum(rx_corr_S1_pat_1_mode3) + sum(rx_corr_S3_pat_1_mode3) > sum(rx_corr_S2_pat_1_mode3) + sum(rx_corr_S4_pat_1_mode3))
        %          tag = 1
        %      else
        %          tag = 2
        %      end

        %  elseif ((sum(rx_corr_S1_pat_2_mode3)>22*num_pulses_mode3) && (sum(rx_corr_S3_pat_2_mode3)>17*num_pulses_mode3)) ...
        %          || ((sum(rx_corr_S2_pat_2_mode3)>20*num_pulses_mode3) && (sum(rx_corr_S4_pat_2_mode3)>19*num_pulses_mode3))
        %      flag_Capture_C = 1;
        %      mode_sel = 3;
        %      if (sum(rx_corr_S1_pat_2_mode3) + sum(rx_corr_S3_pat_2_mode3) > sum(rx_corr_S2_pat_2_mode3) + sum(rx_corr_S4_pat_2_mode3))
        %         tag = 3         
        %      else
        %         tag = 4
        %      end

         end         
    

    % for mode4
    
        if (i+time_frame_mode4*oversamp_IF-1) > length(rx)
            continue;
        end

        % 截取对应96个脉冲长度的波形
        temp_rx_mode4 = rx(i:i+time_frame_mode4*oversamp_IF-1);
        
        % 根据跳频图案对应的频点找到对应通道的波形
        % 两个连续10ms的图案
        wav_temp_mode4 = zeros(num_pulses_mode4, num_bits_pulse * oversamp_IF);
        % 第1个跳时跳频图案
        for pulse_idx = 1:num_pulses_mode4
            f_idx = fh_pat_mode4(pulse_idx);  % 当前脉冲对应的频点
            
            % 计算当前脉冲在一帧中的起始位置
            if pulse_idx == 1
                pos_pulse_mode4 = th_pat_mode4(1)/2;
            else
                pos_pulse_mode4 = sum(th_pat_mode4(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103) + th_pat_mode4(pulse_idx)/2;
            end
            
            % 射频 -> 中频
            if (f_idx >= 1) && (f_idx <= 5)    
                wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % 通道1
            elseif (f_idx >= 6) && (f_idx <= 10)
                wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % 通道2
            elseif (f_idx >= 11) && (f_idx <= 14)
                wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % 通道3
            elseif (f_idx >= 15) && (f_idx <= 18)
                wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % 通道4 
            elseif (f_idx >= 19) && (f_idx <= 21)
                wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % 通道5 
            end

            wav_temp_mode4(pulse_idx,:) = wav_temp1_mode4;

        end


        % % 第2个跳时跳频图案
        % for pulse_idx = 1:num_pulses_mode4
        %     f_idx = fh_pat_2_mode4(pulse_idx);  % 当前脉冲的对应频点
            
        %     % 计算当前脉冲在一帧中的起始位置
        %     if pulse_idx == 1
        %         pos_pulse_mode4 = th_pat_2_mode4(1)/2;
        %     else
        %         pos_pulse_mode4 = sum(th_pat_2_mode4(1:pulse_idx-1)) + (pulse_idx-1)*num_bits_pulse + th_pat_2_mode4(pulse_idx)/2;
        %     end
            
        %     % 射频 -> 中频
        %     if (f_idx >= 1) && (f_idx <= 5)    
        %         wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % 通道1
        %     elseif (f_idx >= 6) && (f_idx <= 10)
        %         wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % 通道2
        %     elseif (f_idx >= 11) && (f_idx <= 14)
        %         wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % 通道3
        %     elseif (f_idx >= 15) && (f_idx <= 18)
        %         wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % 通道4 
        %     elseif (f_idx >= 19) && (f_idx <= 21)
        %         wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % 通道5 
        %     end

        %     wav_temp_mode4(pulse_idx,:,2) = wav_temp1_mode4;

        % end



        % 根据跳频图案对应的频点将对应通道的波形滤波至零中频
        % wav_temp(采样率1024) -> rx_pulse_mat(采样率64)
        rx_pulse_mat_mode4 = downConv(wav_temp_mode4(:,:), num_pulses_mode4, fh_pat_mode4);  % 对应前10ms
        % rx_pulse_mat_2_mode4 = downConv_500K(wav_temp_mode4(:,:,2), num_pulses_mode4, fh_pat_2_mode4);  % 对应后10ms


        % % 预取前后24bit位置的波形 等待后续处理
        % % 前10ms
        % D_S1_mode4 = rx_pulse_mat_mode4(:,1:24*oversamp_BB);
        % D_S1_one_mode4 = D_S1_mode4(:,4:oversamp_BB:end);
        % D_S2_mode4 = rx_pulse_mat_mode4(:,280*oversamp_BB+1:304*oversamp_BB);
        % D_S2_one_mode4 = D_S2_mode4(:,4:oversamp_BB:end);
        % D_S3_mode4 = rx_pulse_mat_mode4(:,24*oversamp_BB+1:45*oversamp_BB);
        % D_S3_one_mode4 = D_S3_mode4(:,4:oversamp_BB:end);
        % D_S4_mode4 = rx_pulse_mat_mode4(:,259*oversamp_BB+1:280*oversamp_BB);
        % D_S4_one_mode4 = D_S4_mode4(:,4:oversamp_BB:end);

        % % 后10ms
        % D_S1_2_mode4 = rx_pulse_mat_2_mode4(:,1:24*oversamp_BB);
        % D_S1_2_one_mode4 = D_S1_2_mode4(:,8:oversamp_BB:end);
        % D_S2_2_mode4 = rx_pulse_mat_2_mode4(:,280*oversamp_BB+1:304*oversamp_BB);
        % D_S2_2_one_mode4 = D_S2_2_mode4(:,8:oversamp_BB:end);
        % D_S3_2_mode4 = rx_pulse_mat_2_mode4(:,24*oversamp_BB+1:45*oversamp_BB);
        % D_S3_2_one_mode4 = D_S3_2_mode4(:,8:oversamp_BB:end);
        % D_S4_2_mode4 = rx_pulse_mat_2_mode4(:,259*oversamp_BB+1:280*oversamp_BB);
        % D_S4_2_one_mode4 = D_S4_2_mode4(:,8:oversamp_BB:end);            


        % 同步头捕获过程    (S1_S3____S4_S2)             
        % 1倍数据速率波形数据做捕获
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                   前10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % 求前后同步头相关峰
        corr_value_mode4 = corr(rx_pulse_mat_mode4, wav_S1_S3_mode4, wav_S4_S2_mode4, 4);
        % rx_corr_S1_pat_1_mode4 = zeros(1,num_pulses_mode4);
        % rx_corr_S2_pat_1_mode4 = zeros(1,num_pulses_mode4);
        % rx_corr_S3_pat_1_mode4 = zeros(1,num_pulses_mode4);
        % rx_corr_S4_pat_1_mode4 = zeros(1,num_pulses_mode4);
        % for j = 1:num_pulses_mode4

        %     rx_wav_S1_pat_1 = D_S1_1_one_mode4(j,1:24) .* conj(wav_S1_S3_1_mode4(j,1:24));     % S1 24bit
        %     rx_wav_S2_pat_1 = D_S2_1_one_mode4(j,1:22) .* conj(wav_S4_S2_1_mode4(j,21+1:43));  % S2 22bit  
        %     rx_wav_S3_pat_1 = D_S3_1_one_mode4(j,1:19) .* conj(wav_S1_S3_1_mode4(j,24+1:43));  % S3 19bit
        %     rx_wav_S4_pat_1 = D_S4_1_one_mode4(j,1:21) .* conj(wav_S4_S2_1_mode4(j,1:21));     % S4 21bit

        %     rx_corr_S1_pat_1_mode4(j) = abs(sum(rx_wav_S1_pat_1));
        %     rx_corr_S2_pat_1_mode4(j) = abs(sum(rx_wav_S2_pat_1));
        %     rx_corr_S3_pat_1_mode4(j) = abs(sum(rx_wav_S3_pat_1));
        %     rx_corr_S4_pat_1_mode4(j) = abs(sum(rx_wav_S4_pat_1));

        % end

        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                       后10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % % 求前后同步头相关峰
        % rx_corr_S1_pat_2_mode4 = zeros(1,num_pulses_mode4);  
        % rx_corr_S2_pat_2_mode4 = zeros(1,num_pulses_mode4);
        % rx_corr_S3_pat_2_mode4 = zeros(1,num_pulses_mode4);
        % rx_corr_S4_pat_2_mode4 = zeros(1,num_pulses_mode4);
        % for j = 1:num_pulses_mode4

        %     rx_wav_S1_pat_2 = D_S1_2_one_mode4(j,1:24) .* conj(wav_S1_S3_2_mode4(j,1:24));     % S1 24bit
        %     rx_wav_S2_pat_2 = D_S2_2_one_mode4(j,1:22) .* conj(wav_S4_S2_2_mode4(j,21+1:43));  % S2 22bit  
        %     rx_wav_S3_pat_2 = D_S3_2_one_mode4(j,1:19) .* conj(wav_S1_S3_2_mode4(j,24+1:43));  % S3 19bit
        %     rx_wav_S4_pat_2 = D_S4_2_one_mode4(j,1:21) .* conj(wav_S4_S2_2_mode4(j,1:21));     % S4 21bit
            
        %     rx_corr_S1_pat_2_mode4(j) = abs(sum(rx_wav_S1_pat_2));
        %     rx_corr_S2_pat_2_mode4(j) = abs(sum(rx_wav_S2_pat_2));
        %     rx_corr_S3_pat_2_mode4(j) = abs(sum(rx_wav_S3_pat_2));
        %     rx_corr_S4_pat_2_mode4(j) = abs(sum(rx_wav_S4_pat_2));

        % end


         % 同步捕获判决条件
        
         if (corr_value_mode4 > 82*num_pulses_mode4) 
             flag_Capture_C = 1;
             mode_sel = 4;
             
        %      if (sum(rx_corr_S1_pat_1_mode4) + sum(rx_corr_S3_pat_1_mode4) > sum(rx_corr_S2_pat_1_mode4) + sum(rx_corr_S4_pat_1_mode4))
        %          tag = 1
        %      else
        %          tag = 2
        %      end
          

        %  elseif ((sum(rx_corr_S1_pat_2_mode4)>22*num_pulses_mode4) && (sum(rx_corr_S3_pat_2_mode4)>17*num_pulses_mode4)) ...
        %          || ((sum(rx_corr_S2_pat_2_mode4)>20*num_pulses_mode4) && (sum(rx_corr_S4_pat_2_mode4)>19*num_pulses_mode4))
        %      flag_Capture_C = 1;
        %      mode_sel = 4;
             
        %      if (sum(rx_corr_S1_pat_2_mode4) + sum(rx_corr_S3_pat_2_mode4) > sum(rx_corr_S2_pat_2_mode4) + sum(rx_corr_S4_pat_2_mode4))
        %         tag = 3         
        %      else
        %         tag = 4
        %      end

        end
   
        disp(i);
        % figure(1);
        % plot(real(D_S1_one_mode1(1,:)));
        % figure(2);
        % plot(real(wav_S1_mode1(1,:)));

    end % (~flag_Capture_C)

    % 捕获成功
    % 根据判定的速率模式送入对应的解调模块
    if (flag_Capture_C == 1)
        if mode_sel == 1
            num_pulses_rec = num_pulses_mode1;  % 一帧的脉冲总数
            length_frame = length_frame_mode1;  % 一帧的有效数据长度
            % time_counter = time_counter_mode1;  % 解调完一帧的等待时长
            % if tag == 1 || tag == 2             % 确定跳频\跳时及同步头图案
                % fh_pat_rec = fh_pat_mode1;
                % th_pat_rec = th_pat_mode1;
                wav_S1_rec = wav_S1_mode1;    
                wav_S2_rec = wav_S2_mode1;    
            % elseif tag == 3 || tag == 4
            %     fh_pat_rec = fh_pat_2_mode1;
            %     th_pat_rec = th_pat_2_mode1;
            %     wav_S1_rec = wav_S1_2_mode1;    
            %     wav_S2_rec = wav_S2_2_mode1; 
            % end
            time_frame = time_frame_mode1;  % 一帧的总长度(含th和103)
            rx_Cap = rx(i:i+time_frame*oversamp_IF+oversamp_IF*2);  %多取了2个bit 
            result_frame = receiver_Cap_mode1(rx_Cap, fh_pat_mode1, th_pat_mode1, wav_S1_S3_mode4, wav_S4_S2_mode4, wav_S1_rec, wav_S2_rec, time_frame, num_pulses_rec);  % 送入2Mbps A模式对应的解调模块解调
            result(ceil((i+(time_frame-0.5)*oversamp_IF)/(time_frame*oversamp_IF)), 1:length_frame) = result_frame;  % 记录解调结果
            
        elseif mode_sel == 2
            num_pulses_rec = num_pulses_mode2;  % 一帧的脉冲总数
            length_frame = length_frame_mode2;  % 一帧的有效数据长度
            % time_counter = time_counter_mode2;  % 解调完一帧的等待时长
            % if tag == 1 || tag == 2             % 确定跳频\跳时及同步头图案
            %     fh_pat_rec = fh_pat_1_mode2;
            %     th_pat_rec = th_pat_1_mode2;
                wav_S1_rec = wav_S1_mode2;    
                wav_S2_rec = wav_S2_mode2;    
            % elseif tag == 3 || tag == 4
            %     fh_pat_rec = fh_pat_2_mode2;
            %     th_pat_rec = th_pat_2_mode2;
            %     wav_S1_rec = wav_S1_2_mode2;    
            %     wav_S2_rec = wav_S2_2_mode2; 
            % end 
            time_frame = time_frame_mode2;      % 一帧的总长度(含跳时)
            rx_Cap = rx(i:i+time_frame*oversamp_IF+oversamp_IF*2);
            result_frame = receiver_Cap_mode2(rx_Cap, fh_pat_mode2, th_pat_mode2, wav_S1_S3_mode4, wav_S4_S2_mode4, wav_S1_rec, wav_S2_rec, time_frame, num_pulses_rec);  % 送入2Mbps B模式对应的解调模块解调
            result(ceil((i+(time_frame-0.5)*oversamp_IF)/(time_frame*oversamp_IF)), 1:length_frame) = result_frame;  % 记录解调结果（取整）
            
        elseif mode_sel == 3
            num_pulses_rec = num_pulses_mode3;  % 一帧的脉冲总数
            length_frame = length_frame_mode3;  % 一帧的有效数据长度
            % time_counter = time_counter_mode3;  % 解调完一帧的等待时长
            % if tag == 1 || tag == 2             % 确定跳频\跳时及同步头图案
            %     fh_pat_rec = fh_pat_1_mode3;
            %     th_pat_rec = th_pat_1_mode3;
                wav_S1_rec = wav_S1_S3_mode3;    
                wav_S2_rec = wav_S4_S2_mode3;    
            % elseif tag == 3 || tag == 4
            %     fh_pat_rec = fh_pat_2_mode3;
            %     th_pat_rec = th_pat_2_mode3;
            %     wav_S1_rec = wav_S1_S3_2_mode3;    
            %     wav_S2_rec = wav_S4_S2_2_mode3; 
            % end 
            time_frame = time_frame_mode3;      % 一帧的总长度(含跳时)
            rx_Cap = rx(i:i+time_frame*oversamp_IF+oversamp_IF*2);
            result_frame = receiver_Cap_mode3(rx_Cap, fh_pat_mode3, th_pat_mode3, wav_S1_S3_mode4, wav_S4_S2_mode4, wav_S1_rec, wav_S2_rec, time_frame, num_pulses_rec);   % 送入500Kbps模式对应的解调模块解调
            result(ceil((i+(time_frame-0.5)*oversamp_IF)/((time_frame+10)*oversamp_IF)), 1:length_frame) = result_frame;  % 记录解调结果
            
        else
            num_pulses_rec = num_pulses_mode4;   % 一帧的脉冲总数
            length_frame = length_frame_mode4;   % 一帧的有效数据长度
            % time_counter = time_counter_mode4;   % 解调完一帧的等待时长
            % if tag == 1 || tag == 2              % 确定跳频\跳时及同步头图案
            %     fh_pat_rec = fh_pat_1_mode4;
            %     th_pat_rec = th_pat_1_mode4;
                % wav_S1_rec = wav_S1_S3_mode4;    
                % wav_S2_rec = wav_S4_S2_mode4;    
            % elseif tag == 3 || tag == 4
            %     fh_pat_rec = fh_pat_2_mode4;
            %     th_pat_rec = th_pat_2_mode4;
            %     wav_S1_rec = wav_S1_S3_2_mode4;    
            %     wav_S2_rec = wav_S4_S2_2_mode4; 
            % end 
            time_frame = time_frame_mode4;      % 一帧的总长度(含跳时)
            rx_Cap = rx(i:i+time_frame*oversamp_IF+oversamp_IF*2);
            result_frame = receiver_Cap_mode4(rx_Cap, fh_pat_mode4, th_pat_mode4, wav_S1_S3_mode4, wav_S4_S2_mode4, time_frame, num_pulses_rec);   % 送入250Kbps模式对应的解调模块解调
            result(ceil((i+(time_frame-0.5)*oversamp_IF)/(time_frame*oversamp_IF)), 1:length_frame) = result_frame;  % 记录解调结果
            
        end
        
        flag_Capture_C = 0;
        frame_counter = frame_counter + 1;  % 捕获到的帧总数加1
        % flag_frame = 1; 
        % tag = 0;
       
        if frame_counter == frame_num+1  %4帧数据已经全部解调完毕
            return 
        end
   
    end

end  % (for i = ...)

figure(3);
k=1:37;
plot(k, corr_value_mode1_plot1);
figure(4);
k=1:37;
plot(k, corr_value_mode1_plot2);