function result = receiver_Cap_mode1(rx, fh_pat, th_pat, wav_S1_S3_mode4, wav_S4_S2_mode4, wav_S1, wav_S2, time_frame, num_pulses)

% 2Mbps A模式 的解调模块

% 参数定义
oversamp_IF = 64;  % 射频、中频信号过采样速率
oversamp_BB = 4;  % 基带信号过采样速率
num_bits_pulse = 304;  % 2Mbps A\500Kbps\250Kbps 一个脉冲包的长度
                       % 2Mbps B 每个脉冲去掉前后跳数据部分后的长度
bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间

num_bits_pn = 24;  % 同步头S1\S2长度

t2 = 0:1/oversamp_IF:num_bits_pulse-1/oversamp_IF;
t2 = t2(1:end) + (T/oversamp_IF/2);  % 中心对称

nn = oversamp_IF / 8; %改不改？

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
% load('lib/filter/LPF.mat');  % 低通滤波器 61阶 通带: 4.8MHz 采样频率128MHz
% load('lib/filter/LPF_2.mat');  % 低通滤波器 507阶 通带: 5MHz  采样频率1024MHz
% S_lpf = 30;
% S_lpf2 = 127;
% S_bpf = 253;
 
Wav_str_Cap_mat = zeros(nn, time_frame*oversamp_IF);  %一行就是一帧的长度
% 间隔8个采样点提取一路波形，共8路
for offset = 1:nn
    Wav_str_Cap_mat(offset,:) = rx(1+(offset-1)*8:time_frame*oversamp_IF+(offset-1)*8);  %offset就是偏移量
end

Wav_str_Cap_F_temp = zeros(num_pulses, num_bits_pulse*oversamp_IF, nn);

Corr_S1 = zeros(1, nn);
Corr_S2 = zeros(1, nn); 
for offset = 1:nn

    temp_rx = Wav_str_Cap_mat(offset,:);

    for pulse_idx = 1:num_pulses

        f_idx = fh_pat(pulse_idx);  % 当前脉冲对应的频点

        % 计算当前脉冲在一帧中的起始位置
        if pulse_idx == 1
            pos_pulse = th_pat(1)/2;
        else
            pos_pulse = sum(th_pat(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103) + th_pat(pulse_idx)/2;
        end

        % 射频 -> 中频
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

        Wav_str_Cap_F_temp(pulse_idx,:,offset) = wav_temp1;

    end

    % 根据跳频图案对应的频点将对应通道的波形滤波至零中频
    % wav_temp(采样率1024) -> rx_pulse_mat(采样率64)
    rx_pulse_mat = downConv(Wav_str_Cap_F_temp(:,:,offset), num_pulses, fh_pat);  

    Corr_mode1(offset) = corr(rx_pulse_mat, wav_S1_S3_mode4, wav_S4_S2_mode4, 1);

    % % 预取前后24bit位置的波形 等待后续处理
    % D_S1 = rx_pulse_mat(:,1:24*oversamp_BB);
    % D_S1_one = D_S1(:,4:oversamp_BB:end);
    % D_S2 = rx_pulse_mat(:,2240+1:2240+24*oversamp_BB);
    % D_S2_one = D_S2(:,4:oversamp_BB:end);


    % % 同步头捕获精细校准  
    % % 利用本地的PN波形做互相关 
    % rx_corr_S1 = zeros(1,num_pulses);
    % rx_corr_S2 = zeros(1,num_pulses);
    % for j = 1:num_pulses

    %     rx_wav_S1 = D_S1_one(j,1:22) .* conj(wav_S1(j,1:22));
    %     rx_wav_S2 = D_S2_one(j,1:22) .* conj(wav_S2(j,1:22));   

    %     rx_corr_S1(j) = abs(sum(rx_wav_S1));
    %     rx_corr_S2(j) = abs(sum(rx_wav_S2));

    % end

    % Corr_S1(offset) = sum(rx_corr_S1);
    % Corr_S2(offset) = sum(rx_corr_S2);

end

pos_Corr = find(Corr_mode1 == max(Corr_mode1));
% pos_Corr_S2 = find(Corr_S2 == max(Corr_S2));

% if (tag == 1) || (tag == 3)
%     pos_Corr = pos_Corr_S1;
% elseif (tag == 2) || (tag == 4)
%     pos_Corr = pos_Corr_S2;
% end


% 取9路波形 间隔1个采样点（再精确）
for offset = 1:9
    Wav_str_Cap_mat(offset,:) = rx((pos_Corr-1)*8+(offset-1)-4:time_frame*oversamp_IF+(pos_Corr-1)*8+(offset-1)-4-1);  %4是8的一半
end            

Wav_str_Cap_F_temp2 = zeros(num_pulses, num_bits_pulse*oversamp_IF, 9);

Corr_S1_Cap_F = zeros(1, 9);
Corr_S2_Cap_F = zeros(1, 9);
for offset = 1:9

    temp_rx = Wav_str_Cap_mat(offset,:);

    for pulse_idx = 1:num_pulses

        f_idx = fh_pat(pulse_idx);  % 当前脉冲的对应的频点

        % 计算当前脉冲在一帧中的起始位置
        if pulse_idx == 1
            pos_pulse = th_pat(1)/2;
        else
            pos_pulse = sum(th_pat(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103) + th_pat(pulse_idx)/2;
        end

        % 射频 -> 中频
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

        Wav_str_Cap_F_temp2(pulse_idx,:,offset) = wav_temp1;

    end

    % 根据跳频图案对应的频点将对应通道的波形滤波至零中频
    % wav_temp(采样率1024) -> rx_pulse_mat(采样率64)
    rx_pulse_mat = downConv(Wav_str_Cap_F_temp2(:,:,offset), num_pulses, fh_pat);

    Corr_Cap_F(offset) = corr(rx_pulse_mat, wav_S1_S3_mode4, wav_S4_S2_mode4, 1);

%     % 预取前后24bit位置的波形 等待后续处理
%     D_S1 = rx_pulse_mat(:,1:24*oversamp_BB);
%     D_S1_one = D_S1(:,8:oversamp_BB:end);
%     D_S2 = rx_pulse_mat(:,2240+1:2240+24*oversamp_BB);
%     D_S2_one = D_S2(:,8:oversamp_BB:end);


%     % 同步头捕获精细校准
%     % 求8路频率对正之后的波形与本地波形的相关值（每路间隔8个采样点）
%     for j = 1:num_pulses

%         rx_wav_S1 = D_S1_one(j,1:22) .* conj(wav_S1(j,1:22));
%         rx_wav_S2 = D_S2_one(j,1:22) .* conj(wav_S2(j,1:22));   

%         rx_corr_S1(j) = abs(sum(rx_wav_S1));
%         rx_corr_S2(j) = abs(sum(rx_wav_S2));

%     end

%     Corr_S1_Cap_F(offset) = sum(rx_corr_S1);
%     Corr_S2_Cap_F(offset) = sum(rx_corr_S2);

end

pos_Corr_F = find(Corr_Cap_F == max(Corr_Cap_F));
% pos_Corr_S2_F = find(Corr_S2_Cap_F == max(Corr_S2_Cap_F));

offset_i = pos_Corr;
offset_j = pos_Corr_F;

% if (tag == 1) || (tag == 3)
%     offset_i = pos_Corr_S1
%     offset_j = pos_Corr_S1_F
%     choose_flag = 1;
% elseif (tag == 2) || (tag == 4)
%     offset_i = pos_Corr_S2
%     offset_j = pos_Corr_S2_F
%     choose_flag = 2;

% end

temp_rx_FNL = rx((offset_i-1)*8+(offset_j-1)-4:time_frame*oversamp_IF+(offset_i-1)*8+(offset_j-1)-4-1);
for pulse_idx = 1:num_pulses

    f_idx = fh_pat(pulse_idx);  % 当前脉冲对应的频点
    
    % 计算脉冲在一帧中的起始位置
    if pulse_idx == 1
        pos_pulse = th_pat(1)/2;
    else
        pos_pulse = sum(th_pat(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103) + th_pat(pulse_idx)/2;
    end

    % 射频 -> 中频
    if (f_idx >= 1) && (f_idx <= 5)    
        wav_temp1 = downConv_IF(temp_rx_FNL(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % 通道1
    elseif (f_idx >= 6) && (f_idx <= 10)
        wav_temp1 = downConv_IF(temp_rx_FNL(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % 通道2
    elseif (f_idx >= 11) && (f_idx <= 14)
        wav_temp1 = downConv_IF(temp_rx_FNL(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % 通道3
    elseif (f_idx >= 15) && (f_idx <= 18)
        wav_temp1 = downConv_IF(temp_rx_FNL(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % 通道4 
    elseif (f_idx >= 19) && (f_idx <= 21)
        wav_temp1 = downConv_IF(temp_rx_FNL(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % 通道5 
    end

    wav_temp_FNL(pulse_idx,:) = wav_temp1;

end

% 根据跳频图案对应的频点将对应通道的波形滤波至零中频
% wav_temp(采样率1024) -> rx_pulse_mat(采样率64)
rx_pulse_mat_FNL = downConv(wav_temp_FNL(:,:), num_pulses, fh_pat); 

% 预取前后24bit位置的波形 等待后续处理
D_S1 = rx_pulse_mat_FNL(:,1:24*oversamp_BB);
D_S1_one = D_S1(:,4:oversamp_BB:end);
D_S2 = rx_pulse_mat_FNL(:,1120+1:1120+24*oversamp_BB);
D_S2_one = D_S2(:,4:oversamp_BB:end);

zk_S1 = zeros(num_pulses, (num_bits_pn-2));
zk_S2 = zeros(num_pulses, (num_bits_pn-2));
for j = 1:num_pulses                  
    zk_S1(j,:) = D_S1_one(j,1:num_bits_pn-2) .* conj(wav_S1(j,1:num_bits_pn-2));
    zk_S2(j,:) = D_S2_one(j,1:num_bits_pn-2) .* conj(wav_S2(j,1:num_bits_pn-2)); 
end

% 脉冲完整
for j = 1:num_pulses 
    
    zk_corr_S12(j) = sum(conj(zk_S1(j,:)) .* zk_S2(j,:));
    % zk_corr_S12 = conj(zk_S1(j,:)) .* zk_S2(j,:);
    % deltat_hat_S12(j) = (atan2(-imag(sum(zk_corr_S12)), real(sum(zk_corr_S12)))); 
    % delta_f_hat_S12(j) = deltat_hat_S12(j)/(2*pi*280*T);

end
DF_hat_S12 = atan2(-imag(mean(zk_corr_S12)), real(mean(zk_corr_S12)))/(2*pi*280*T);

% 同步头消除频偏
counter_f_S1 = repmat(DF_hat_S12,[num_pulses,1]) * 2 * pi * t2(16:16:24*oversamp_IF) * T; %这第一个16怎么定？
counter_f_S2 = repmat(DF_hat_S12,[num_pulses,1]) * 2 * pi * t2(280*oversamp_IF+16:16:280*oversamp_IF+24*oversamp_IF) * T;
D_S1_half = D_S1 .* complex(cos(counter_f_S1), sin(counter_f_S1));  % 对前后同步头位置的波形消除频偏
D_S2_half = D_S2 .* complex(cos(counter_f_S2), sin(counter_f_S2));

% 再次利用本地的PN波形计算相偏
delta_theta = zeros(num_pulses, 1);
% delta_theta_S2 = zeros(num_pulses, 1);
for j = 1:num_pulses
    delta_theta_cplx_mat_S1 = D_S1_half(j,4:oversamp_BB:end) .* conj(wav_S1(j,:));
    delta_theta_cplx_mat_S2 = D_S2_half(j,4:oversamp_BB:end) .* conj(wav_S2(j,:));
    delta_theta_cplx = (sum(delta_theta_cplx_mat_S1(1:22)) + sum(delta_theta_cplx_mat_S2(1:22)))/(22*2); %用22个？
    delta_theta(j) = atan2(-imag(delta_theta_cplx), real(delta_theta_cplx));
    % delta_theta_cplx_S1 = sum(delta_theta_cplx_mat_S1(1:22)) / 22;
    % delta_theta_cplx_S2 = sum(delta_theta_cplx_mat_S2(1:22)) / 22;
    % delta_theta_S1(j) = atan2(-imag(delta_theta_cplx_S1), real(delta_theta_cplx_S1));
    % delta_theta_S2(j) = atan2(-imag(delta_theta_cplx_S2), real(delta_theta_cplx_S2));
end

delta_f = DF_hat_S12;
% if (choose_flag == 1)
%     delta_theta = delta_theta_S1;
%     choose_flag = 0;
% elseif (choose_flag == 2)
%     delta_theta = delta_theta_S2;
%     choose_flag = 0;
% end


% 对当前帧的波形利用之前计算所得的频偏及相偏估计值做波形校正，准备送入解调模块（4倍速率）
counter_f = delta_f * 2 * pi * t2(16:16:end) * T;
rx_pulse_mat_FNL = rx_pulse_mat_FNL .* repmat(complex(cos(counter_f), sin(counter_f)),[num_pulses,1]) .* repmat(complex(cos(delta_theta), sin(delta_theta)), [1,num_bits_pulse * oversamp_BB]); 

%只取出数据，舍掉同步头（4倍速率），用于解调
rx_pulse_mat_FNL1 = zreos(num_pulses, 256*oversamp_BB);
rx_pulse_mat_FNL1 = rx_pulse_mat_FNL(:, 24*oversamp_BB+1:280*oversamp_BB);

% GMSK 解调（4倍速率）
out_temp = zeros(num_pulses, 279);
for pulse_idx = 1:num_pulses
    out_temp(pulse_idx,:) = GMSK_demod(rx_pulse_mat_FNL(pulse_idx, 1:(24+256)*oversamp_BB), c0_f, c1_f, oversamp_BB, 0, 280);
end

% 解调结果组成完整的一帧
for pulse_idx = 1:num_pulses
    result((pulse_idx-1)*256+1:pulse_idx*256) = out_temp(pulse_idx, 24:end);  %只有数据
end
