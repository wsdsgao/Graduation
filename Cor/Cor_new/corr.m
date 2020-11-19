function corr_value = corr(rx, wav_S1_S3_mode4, wav_S4_S2_mode4, mode)

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

    % 各模式一帧脉冲数量
    num_pulses_mode1 = 12;
    num_pulses_mode2 = 12;
    num_pulses_mode3 = 48;
    num_pulses_mode4 = 96;

    if mode==1
        % 预取前后24bit位置的波形（即同步头） 等待后续处理
        rx_pulse_mat_mode1 = rx;
        D_S1_mode1 = rx_pulse_mat_mode1(:,1:24*oversamp_BB);  %一行是一个脉冲
        D_S1_one_mode1 = D_S1_mode1(:,4:oversamp_BB:end);  %和数据速率相同，一个点对应一个bit
        D_S2_mode1 = rx_pulse_mat_mode1(:,1120+1:1120+24*oversamp_BB);  %1120=280*oversamp_BB
        D_S2_one_mode1 = D_S2_mode1(:,4:oversamp_BB:end);

        % 求前后同步头相关峰并相加
        rx_corr_S1_pat_mode1 = zeros(1,num_pulses_mode1);
        rx_corr_S2_pat_mode1 = zeros(1,num_pulses_mode1);
        rx_corr_pat_mode1 = zeros(1,num_pulses_mode1);
        for j = 1:num_pulses_mode1

            rx_wav_S1_pat_mode1_1 = D_S1_one_mode1(j,1:8) .* (wav_S1_S3_mode4(j,1:8));
            rx_wav_S1_pat_mode1_2 = D_S1_one_mode1(j,9:16) .* conj(wav_S1_S3_mode4(j,9:16));
            rx_wav_S1_pat_mode1_3 = D_S1_one_mode1(j,17:24) .* conj(wav_S1_S3_mode4(j,17:24)); %考虑一下是否为1-22比较合理
            rx_wav_S2_pat_mode1_1 = D_S2_one_mode1(j,1:8) .* (wav_S4_S2_mode4(j,22:29)); %考虑一下是否为24-45比较合理
            rx_wav_S2_pat_mode1_2 = D_S2_one_mode1(j,9:16) .* conj(wav_S4_S2_mode4(j,30:37));
            rx_wav_S2_pat_mode1_3 = D_S2_one_mode1(j,17:24) .* conj(wav_S4_S2_mode4(j,38:45));   

            rx_corr_S1_pat_mode1(j) = abs(sum(rx_wav_S1_pat_mode1_1)+sum(rx_wav_S1_pat_mode1_2)+sum(rx_wav_S1_pat_mode1_3));
            rx_corr_S2_pat_mode1(j) = abs(sum(rx_wav_S2_pat_mode1_1)+sum(rx_wav_S2_pat_mode1_2)+sum(rx_wav_S2_pat_mode1_3));

            rx_corr_pat_mode1(j) =  rx_corr_S1_pat_mode1(j) + rx_corr_S2_pat_mode1(j);

        end

        corr_value = sum(rx_corr_pat_mode1);
    end

    if mode==2
        % 预取前后24bit位置的波形（即同步头） 等待后续处理
        rx_pulse_mat_mode2 = rx;
        D_S1_mode2 = rx_pulse_mat_mode2(:,1:24*oversamp_BB);  %一行是一个脉冲
        D_S1_one_mode2 = D_S1_mode2(:,4:oversamp_BB:end);  %和数据速率相同，一个点对应一个bit
        D_S2_mode2 = rx_pulse_mat_mode2(:,1120+1:1120+24*oversamp_BB);  %2240=280*oversamp_BB
        D_S2_one_mode2 = D_S2_mode2(:,4:oversamp_BB:end);

        % 求前后同步头相关峰并相加
        rx_corr_S1_pat_mode2 = zeros(1,num_pulses_mode2);
        rx_corr_S2_pat_mode2 = zeros(1,num_pulses_mode2);
        rx_corr_pat_mode2 = zeros(1,num_pulses_mode2);
        for j = 1:num_pulses_mode2

            rx_wav_S1_pat_mode2_1 = D_S1_one_mode2(j,1:8) .* conj(wav_S1_S3_mode4(j,1:8));
            rx_wav_S1_pat_mode2_2 = D_S1_one_mode2(j,9:16) .* (wav_S1_S3_mode4(j,9:16));
            rx_wav_S1_pat_mode2_3 = D_S1_one_mode2(j,17:24) .* conj(wav_S1_S3_mode4(j,17:24)); %考虑一下是否为1-22比较合理
            rx_wav_S2_pat_mode2_1 = D_S2_one_mode2(j,1:8) .* conj(wav_S4_S2_mode4(j,22:29)); %考虑一下是否为24-45比较合理
            rx_wav_S2_pat_mode2_2 = D_S2_one_mode2(j,9:16) .* (wav_S4_S2_mode4(j,30:37));
            rx_wav_S2_pat_mode2_3 = D_S2_one_mode2(j,17:24) .* conj(wav_S4_S2_mode4(j,38:45));  

            rx_corr_S1_pat_mode2(j) = abs(sum(rx_wav_S1_pat_mode2_1)+sum(rx_wav_S1_pat_mode2_2)+sum(rx_wav_S1_pat_mode2_3));
            rx_corr_S2_pat_mode2(j) = abs(sum(rx_wav_S2_pat_mode2_1)+sum(rx_wav_S2_pat_mode2_2)+sum(rx_wav_S2_pat_mode2_3));

            rx_corr_pat_mode2(j) =  rx_corr_S1_pat_mode2(j) + rx_corr_S2_pat_mode2(j);

        end

        corr_value = sum(rx_corr_pat_mode2);
    end

    if mode==3
        % 预取前后同步头波形 等待后续处理
        rx_pulse_mat_mode3 = rx;
        D_S1_mode3 = rx_pulse_mat_mode3(:,1:24*oversamp_BB);
        D_S1_one_mode3 = D_S1_mode3(:,4:oversamp_BB:end);
        D_S2_mode3 = rx_pulse_mat_mode3(:,280*oversamp_BB+1:304*oversamp_BB);
        D_S2_one_mode3 = D_S2_mode3(:,4:oversamp_BB:end);
        D_S3_mode3 = rx_pulse_mat_mode3(:,24*oversamp_BB+1:45*oversamp_BB);
        D_S3_one_mode3 = D_S3_mode3(:,4:oversamp_BB:end);
        D_S4_mode3 = rx_pulse_mat_mode3(:,259*oversamp_BB+1:280*oversamp_BB);
        D_S4_one_mode3 = D_S4_mode3(:,4:oversamp_BB:end);

        % 求前后同步头相关峰并相加
        rx_corr_S1_pat_mode3 = zeros(1,num_pulses_mode3);
        rx_corr_S2_pat_mode3 = zeros(1,num_pulses_mode3);
        rx_corr_S3_pat_mode3 = zeros(1,num_pulses_mode3);
        rx_corr_S4_pat_mode3 = zeros(1,num_pulses_mode3);
        rx_corr_pat_mode3 = zeros(1,num_pulses_mode3);
        for j = 1:num_pulses_mode3

            rx_wav_S1_pat_mode3 = D_S1_one_mode3(j,1:24) .* conj(wav_S1_S3_mode4(j,1:24));
            rx_wav_S3_pat_mode3 = D_S3_one_mode3(j,1:21) .* (wav_S1_S3_mode4(j,25:45));
            rx_wav_S4_pat_mode3 = D_S4_one_mode3(j,1:21) .* (wav_S4_S2_mode4(j,1:21));
            rx_wav_S2_pat_mode3 = D_S2_one_mode3(j,1:24) .* conj(wav_S4_S2_mode4(j,22:45));

            rx_corr_S1_pat_mode3(j) = abs(sum(rx_wav_S1_pat_mode3));
            rx_corr_S2_pat_mode3(j) = abs(sum(rx_wav_S2_pat_mode3));
            rx_corr_S3_pat_mode3(j) = abs(sum(rx_wav_S3_pat_mode3));
            rx_corr_S4_pat_mode3(j) = abs(sum(rx_wav_S4_pat_mode3));

            rx_corr_pat_mode3(j) =  rx_corr_S1_pat_mode3(j) + rx_corr_S2_pat_mode3(j) + rx_corr_S3_pat_mode3(j) + rx_corr_S4_pat_mode3(j);

        end

        corr_value = sum(rx_corr_pat_mode3);
    end

    if mode==4
        % 预取前后同步头波形 等待后续处理
        rx_pulse_mat_mode4 = rx;
        D_S1_mode4 = rx_pulse_mat_mode4(:,1:24*oversamp_BB);
        D_S1_one_mode4 = D_S1_mode4(:,4:oversamp_BB:end);
        D_S2_mode4 = rx_pulse_mat_mode4(:,280*oversamp_BB+1:304*oversamp_BB);
        D_S2_one_mode4 = D_S2_mode4(:,4:oversamp_BB:end);
        D_S3_mode4 = rx_pulse_mat_mode4(:,24*oversamp_BB+1:45*oversamp_BB);
        D_S3_one_mode4 = D_S3_mode4(:,4:oversamp_BB:end);
        D_S4_mode4 = rx_pulse_mat_mode4(:,259*oversamp_BB+1:280*oversamp_BB);
        D_S4_one_mode4 = D_S4_mode4(:,4:oversamp_BB:end);

        % 求前后同步头相关峰并相加
        rx_corr_S1_pat_mode4 = zeros(1,num_pulses_mode4);
        rx_corr_S2_pat_mode4 = zeros(1,num_pulses_mode4);
        rx_corr_S3_pat_mode4 = zeros(1,num_pulses_mode4);
        rx_corr_S4_pat_mode4 = zeros(1,num_pulses_mode4);
        rx_corr_pat_mode4 = zeros(1,num_pulses_mode4);
        for j = 1:num_pulses_mode4

            rx_wav_S1_pat_mode4 = D_S1_one_mode4(j,1:24) .* conj(wav_S1_S3_mode4(j,1:24));
            rx_wav_S3_pat_mode4 = D_S3_one_mode4(j,1:21) .* conj(wav_S1_S3_mode4(j,25:45));
            rx_wav_S4_pat_mode4 = D_S4_one_mode4(j,1:21) .* conj(wav_S4_S2_mode4(j,1:21));
            rx_wav_S2_pat_mode4 = D_S2_one_mode4(j,1:24) .* conj(wav_S4_S2_mode4(j,22:45));

            rx_corr_S1_pat_mode4(j) = abs(sum(rx_wav_S1_pat_mode4));
            rx_corr_S2_pat_mode4(j) = abs(sum(rx_wav_S2_pat_mode4));
            rx_corr_S3_pat_mode4(j) = abs(sum(rx_wav_S3_pat_mode4));
            rx_corr_S4_pat_mode4(j) = abs(sum(rx_wav_S4_pat_mode4));

            rx_corr_pat_mode4(j) =  rx_corr_S1_pat_mode4(j) + rx_corr_S2_pat_mode4(j) + rx_corr_S3_pat_mode4(j) + rx_corr_S4_pat_mode4(j);

        end

        corr_value = sum(rx_corr_pat_mode4);
    end
