function rx_pulse_mat = downConv_th(wav, num_pulses, fh_pat, th_pat, indi)

% 根据跳频\跳时图案和一帧中包含的总脉冲个数，将接收的中频信号波形下变频到基带并从64倍符号速率降为8倍符号速率
% 用于2Mbps B模式的前后跳数据部分
% indi  用于指示当前是前跳数据部分还是后跳数据部分

% 参数定义
bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间
fs_IF = 1024e6;  % 中频信号采样速率
fs_BB = 64e6;  % 基带信号采样速率
oversamp_BB = T * fs_BB;  % 采样速率为 8倍符号速率
oversamp_IF = T * fs_IF;
num_bits_pulse = 304;

% 低通滤波器
load('filter/LPF_1.mat'); % 采样率1024MHz，通带20MHz
load('filter/LPF_2.mat'); % 采样率256MHz，通带10MHz
S_lpf1 = 32;
S_lpf2 = 32;


% 21个频点下变频至零中频所需各点频率
f_chan12_1 = 213.34e6;  % 通道1、2 第1个频点对应的频率
f_chan12_2 = 226.67e6;  % 通道1、2 第2个频点对应的频率
f_chan12_3 = 240e6;     % 通道1、2 第3个频点对应的频率
f_chan12_4 = 253.33e6;  % 通道1、2 第4个频点对应的频率
f_chan12_5 = 266.66e6;  % 通道1、2 第5个频点对应的频率
f_chan34_1 = 220.005e6; % 通道3、4 第1个频点对应的频率
f_chan34_2 = 233.335e6; % 通道3、4 第2个频点对应的频率
f_chan34_3 = 246.665e6; % 通道3、4 第3个频点对应的频率
f_chan34_4 = 259.995e6; % 通道3、4 第4个频点对应的频率
f_chan5_1 = f_chan12_2; % 通道5 第1个频点对应的频率
f_chan5_2 = f_chan12_3; % 通道5 第2个频点对应的频率
f_chan5_3 = f_chan12_4; % 通道5 第3个频点对应的频率

rx_pulse_mat = zeros(num_pulses, 256 * oversamp_BB);
for pulse_idx = 1:num_pulses
    
    f_idx = fh_pat(pulse_idx);  % 按脉冲序号和跳频图案，找到当前脉冲对应的频点
    th = th_pat(pulse_idx)/2;  % 按脉冲序号和跳时图案，找到当前脉冲对应的跳数据长度
    
    t_240_chip = -(num_bits_pulse + th_pat(pulse_idx))*T/2:1/fs_IF:(num_bits_pulse + th_pat(pulse_idx))*T/2-1/fs_IF;
    t_240_pulse_temp = t_240_chip(1:end) + (T/oversamp_IF/2);  % 中心对称;
    
    if (indi == 1)  % 首跳时数据部分
        range = 1:th*oversamp_IF;
        range2 = 1:th*oversamp_BB;
        t_240_pulse = t_240_pulse_temp(1:th*oversamp_IF);
    elseif (indi == 2)  % 尾跳时数据部分
        range = th*oversamp_IF+1:(2*th)*oversamp_IF;
        range2 = 1:th*oversamp_BB;
        t_240_pulse = t_240_pulse_temp((th+304)*oversamp_IF+1:end);
    end
    
    switch f_idx
                        
        % 变换至零中频
        case 1  % 频点1
            rx_f1_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan12_1*t_240_pulse)); % 取差频
            rx_f1_256_temp = conv(rx_f1_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f1_256 = downsample(rx_f1_256_temp(S_lpf1+1:S_lpf1+length(rx_f1_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f1_64_temp = conv(rx_f1_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f1_64 = downsample(rx_f1_64_temp(S_lpf2+1:S_lpf2+length(rx_f1_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f1_64;



        case 2  % 频点2
            rx_f2_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan12_2*t_240_pulse)); % 取差频
            rx_f2_256_temp = conv(rx_f2_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f2_256 = downsample(rx_f2_256_temp(S_lpf1+1:S_lpf1+length(rx_f2_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f2_64_temp = conv(rx_f2_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f2_64 = downsample(rx_f2_64_temp(S_lpf2+1:S_lpf2+length(rx_f2_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f2_64;



        case 3  % 频点3
            rx_f3_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan12_3*t_240_pulse)); % 取差频
            rx_f3_256_temp = conv(rx_f3_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f3_256 = downsample(rx_f3_256_temp(S_lpf1+1:S_lpf1+length(rx_f3_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f3_64_temp = conv(rx_f3_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f3_64 = downsample(rx_f3_64_temp(S_lpf2+1:S_lpf2+length(rx_f3_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f3_64;



        case 4  % 频点4
            rx_f4_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan12_4*t_240_pulse)); % 取差频
            rx_f4_256_temp = conv(rx_f4_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f4_256 = downsample(rx_f4_256_temp(S_lpf1+1:S_lpf1+length(rx_f4_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f4_64_temp = conv(rx_f4_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f4_64 = downsample(rx_f4_64_temp(S_lpf2+1:S_lpf2+length(rx_f4_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f4_64;



        case 5  % 频点5
            rx_f5_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan12_5*t_240_pulse)); % 取差频
            rx_f5_256_temp = conv(rx_f5_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f5_256 = downsample(rx_f5_256_temp(S_lpf1+1:S_lpf1+length(rx_f5_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f5_64_temp = conv(rx_f5_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f5_64 = downsample(rx_f5_64_temp(S_lpf2+1:S_lpf2+length(rx_f5_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f5_64;



        case 6  % 频点6                            
            rx_f6_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan12_1*t_240_pulse)); % 取差频
            rx_f6_256_temp = conv(rx_f6_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f6_256 = downsample(rx_f6_256_temp(S_lpf1+1:S_lpf1+length(rx_f6_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f6_64_temp = conv(rx_f6_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f6_64 = downsample(rx_f6_64_temp(S_lpf2+1:S_lpf2+length(rx_f6_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f6_64;



        case 7  % 频点7
            rx_f7_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan12_2*t_240_pulse)); % 取差频
            rx_f7_256_temp = conv(rx_f7_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f7_256 = downsample(rx_f7_256_temp(S_lpf1+1:S_lpf1+length(rx_f7_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f7_64_temp = conv(rx_f7_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f7_64 = downsample(rx_f7_64_temp(S_lpf2+1:S_lpf2+length(rx_f7_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f7_64;



        case 8  % 频点8    
            rx_f8_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan12_3*t_240_pulse)); % 取差频
            rx_f8_256_temp = conv(rx_f8_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f8_256 = downsample(rx_f8_256_temp(S_lpf1+1:S_lpf1+length(rx_f8_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f8_64_temp = conv(rx_f8_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f8_64 = downsample(rx_f8_64_temp(S_lpf2+1:S_lpf2+length(rx_f8_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f8_64;
     


        case 9  % 频点9  
            rx_f9_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan12_4*t_240_pulse)); % 取差频
            rx_f9_256_temp = conv(rx_f9_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f9_256 = downsample(rx_f9_256_temp(S_lpf1+1:S_lpf1+length(rx_f9_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f9_64_temp = conv(rx_f9_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f9_64 = downsample(rx_f9_64_temp(S_lpf2+1:S_lpf2+length(rx_f9_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f9_64;



        case 10 % 频点10   
            rx_f10_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan12_5*t_240_pulse)); % 取差频
            rx_f10_256_temp = conv(rx_f10_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f10_256 = downsample(rx_f10_256_temp(S_lpf1+1:S_lpf1+length(rx_f10_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f10_64_temp = conv(rx_f10_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f10_64 = downsample(rx_f10_64_temp(S_lpf2+1:S_lpf2+length(rx_f10_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f10_64;



        case 11 % 频点11                           
            rx_f11_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan34_1*t_240_pulse)); % 取差频
            rx_f11_256_temp = conv(rx_f11_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f11_256 = downsample(rx_f11_256_temp(S_lpf1+1:S_lpf1+length(rx_f11_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f11_64_temp = conv(rx_f11_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f11_64 = downsample(rx_f11_64_temp(S_lpf2+1:S_lpf2+length(rx_f11_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f11_64;



        case 12 % 频点12
            rx_f12_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan34_2*t_240_pulse)); % 取差频
            rx_f12_256_temp = conv(rx_f12_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f12_256 = downsample(rx_f12_256_temp(S_lpf1+1:S_lpf1+length(rx_f12_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f12_64_temp = conv(rx_f12_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f12_64 = downsample(rx_f12_64_temp(S_lpf2+1:S_lpf2+length(rx_f12_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f12_64;



        case 13 % 频点13
            rx_f13_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan34_3*t_240_pulse)); % 取差频
            rx_f13_256_temp = conv(rx_f13_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f13_256 = downsample(rx_f13_256_temp(S_lpf1+1:S_lpf1+length(rx_f13_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f13_64_temp = conv(rx_f13_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f13_64 = downsample(rx_f13_64_temp(S_lpf2+1:S_lpf2+length(rx_f13_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f13_64;



        case 14 % 频点14
            rx_f14_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan34_4*t_240_pulse)); % 取差频
            rx_f14_256_temp = conv(rx_f14_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f14_256 = downsample(rx_f14_256_temp(S_lpf1+1:S_lpf1+length(rx_f14_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f14_64_temp = conv(rx_f14_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f14_64 = downsample(rx_f14_64_temp(S_lpf2+1:S_lpf2+length(rx_f14_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f14_64;



        case 15 % 频点15
            rx_f15_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan34_1*t_240_pulse)); % 取差频
            rx_f15_256_temp = conv(rx_f15_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f15_256 = downsample(rx_f15_256_temp(S_lpf1+1:S_lpf1+length(rx_f15_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f15_64_temp = conv(rx_f15_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f15_64 = downsample(rx_f15_64_temp(S_lpf2+1:S_lpf2+length(rx_f15_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f15_64;



        case 16 % 频点16
            rx_f16_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan34_2*t_240_pulse)); % 取差频
            rx_f16_256_temp = conv(rx_f16_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f16_256 = downsample(rx_f16_256_temp(S_lpf1+1:S_lpf1+length(rx_f16_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f16_64_temp = conv(rx_f16_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f16_64 = downsample(rx_f16_64_temp(S_lpf2+1:S_lpf2+length(rx_f16_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f16_64;



        case 17 % 频点17
            rx_f17_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan34_3*t_240_pulse)); % 取差频
            rx_f17_256_temp = conv(rx_f17_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f17_256 = downsample(rx_f17_256_temp(S_lpf1+1:S_lpf1+length(rx_f17_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f17_64_temp = conv(rx_f17_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f17_64 = downsample(rx_f17_64_temp(S_lpf2+1:S_lpf2+length(rx_f17_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f17_64;
                            


        case 18 % 频点18
            rx_f18_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan34_4*t_240_pulse)); % 取差频
            rx_f18_256_temp = conv(rx_f18_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f18_256 = downsample(rx_f18_256_temp(S_lpf1+1:S_lpf1+length(rx_f18_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f18_64_temp = conv(rx_f18_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f18_64 = downsample(rx_f18_64_temp(S_lpf2+1:S_lpf2+length(rx_f18_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f18_64;



        case 19 % 频点19
            rx_f19_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan5_1*t_240_pulse)); % 取差频
            rx_f19_256_temp = conv(rx_f19_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f19_256 = downsample(rx_f19_256_temp(S_lpf1+1:S_lpf1+length(rx_f19_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f19_64_temp = conv(rx_f19_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f19_64 = downsample(rx_f19_64_temp(S_lpf2+1:S_lpf2+length(rx_f19_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f19_64;
                          


        case 20 % 频点20
            rx_f20_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan5_2*t_240_pulse)); % 取差频
            rx_f20_256_temp = conv(rx_f20_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f20_256 = downsample(rx_f20_256_temp(S_lpf1+1:S_lpf1+length(rx_f20_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f20_64_temp = conv(rx_f20_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f20_64 = downsample(rx_f20_64_temp(S_lpf2+1:S_lpf2+length(rx_f20_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f20_64;



        case 21 % 频点21
            rx_f21_temp1 = wav(pulse_idx,range) .* exp(1i*(-2*pi*f_chan5_3*t_240_pulse)); % 取差频
            rx_f21_256_temp = conv(rx_f21_temp1, LPF_1);  % 低通滤波(第1次抽取前)
            rx_f21_256 = downsample(rx_f21_256_temp(S_lpf1+1:S_lpf1+length(rx_f21_temp1)), 4, 3);  % 4倍抽取， 采样率降为256（3是偏移量）
            rx_f21_64_temp = conv(rx_f21_256, LPF_2);  % 低通滤波(第2次抽取前)
            rx_f21_64 = downsample(rx_f21_64_temp(S_lpf2+1:S_lpf2+length(rx_f21_256)), 4, 3);  % 4倍抽取， 采样率降为64（3是偏移量）
            rx_pulse_mat(pulse_idx, range2) = rx_f21_64;
  
    end
                
    
end