function rx_pulse_mat = DDC(temp_rx, num_pulses, fh_pat, th_pat, mode)

    num_bits_pulse = 304;
    oversamp_IF = 64;

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
    load('filter/LPF_1.mat'); % 采样率1024MHz，通带20MHz
    load('filter/LPF_2.mat'); % 采样率256MHz，通带10MHz
    S_lpf1 = 32;
    S_lpf2 = 32;
    S_bpf = 253;

    if mode == 1 || mode ==3 || mode ==4
        % 根据跳频图案对应的频点找到对应通道的波形
        wav_temp = zeros(num_pulses, num_bits_pulse * oversamp_IF);

        for pulse_idx = 1:num_pulses
            f_idx = fh_pat(pulse_idx); % 当前脉冲对应频点
            
            % 计算当前脉冲在一帧中的起始位置
            if pulse_idx == 1
                pos_pulse = th_pat(1)/2;
            else
                pos_pulse = sum(th_pat(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103) + th_pat(pulse_idx)/2;
            end

            % 射频 -> 中频 （只取出数据位）
            if (f_idx >= 1) && (f_idx <= 5)
                wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % 通道1
            elseif (f_idx >= 6) && (f_idx <= 10)
                wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % 通道2
            elseif (f_idx >= 11) && (f_idx <= 14)
                wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % 通道3
            elseif (f_idx >= 15) && (f_idx <= 18)
                wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % 通道4 
            elseif (f_idx >= 19) && (f_idx <= 21)
                wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % 通道5 
            end

            wav_temp(pulse_idx,:) = wav_temp1;

        end


        % 根据跳频图案对应的频点将对应通道的波形滤波至零中频
        % wav_temp(采样率1024) -> rx_pulse_mat(采样率64)
        rx_pulse_mat = downConv(wav_temp(:,:), num_pulses, fh_pat);
    
    else % if mode == 2
        % 根据跳频图案对应的频点找到对应通道的波形
        wav_temp = zeros(num_pulses, num_bits_pulse * oversamp_IF);

        for pulse_idx = 1:num_pulses
            f_idx = fh_pat(pulse_idx); % 当前脉冲对应频点
            
            % 计算当前脉冲在一帧中的起始位置
            if pulse_idx == 1
                pos_pulse = 0;
            else
                pos_pulse = sum(th_pat(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103);
            end

            % 计算当前脉冲的长度
            num_bits_pulse_rx = th_pat(pulse_idx) + num_bits_pulse; %包含了跳时

            % 射频 -> 中频 （只取出数据位）
            if (f_idx >= 1) && (f_idx <= 5)
                wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % 通道1
            elseif (f_idx >= 6) && (f_idx <= 10)
                wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % 通道2
            elseif (f_idx >= 11) && (f_idx <= 14)
                wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % 通道3
            elseif (f_idx >= 15) && (f_idx <= 18)
                wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % 通道4 
            elseif (f_idx >= 19) && (f_idx <= 21)
                wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % 通道5 
            end

            wav_temp(pulse_idx,:) = wav_temp1((th_pat(pulse_idx)/2)*oversamp_IF+1:(th_pat(pulse_idx)/2+num_bits_pulse)*oversamp_IF);  %没包含跳时的数据

        end


        % 根据跳频图案对应的频点将对应通道的波形滤波至零中频
        % wav_temp(采样率1024) -> rx_pulse_mat(采样率64)
        rx_pulse_mat = downConv_2MLL(wav_temp(:,:), num_pulses, fh_pat);

    end