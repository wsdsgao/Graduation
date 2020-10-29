function rx_pulse_mat = downConv_2MLL(wav, num_pulses, fh_pat, th_pat)

% 根据跳频\跳时图案和一帧中包含的总脉冲个数，将接收的中频信号波形下变频到基带并从64倍符号速率降为8倍符号速率
% 用于2Mbps B模式

% 参数定义
bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间
fs_IF = 1024e6;  % 中频信号采样速率
fs_BB = 128e6;  % 基带信号采样速率
oversamp_BB = T * fs_BB;  % 基带信号过采样速率
oversamp_IF = T * fs_IF;  % 中频信号过采样速率
num_bits_pulse = 304;  % 2Mbps A\500Kbps\250Kbps 一个脉冲包的长度
                       % 2Mbps B 每个脉冲去掉前后跳数据部分后的长度

% 低通滤波器
load('lib/filter/LPF_2.mat'); % 采样率1024MHz，通带8MHz
load('lib/filter/LPF.mat'); % 采样率128MHz，通带5MHz
S_lpf = 30;
S_lpf2 = 127;

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

rx_pulse_mat = zeros(num_pulses, num_bits_pulse * oversamp_BB);
for pulse_idx = 1:num_pulses
    
    f_idx = fh_pat(pulse_idx);  % 按脉冲序号和跳频图案，找到当前脉冲对应的频点
    
    th = th_pat(pulse_idx) + 3; % 按脉冲序号和跳时图案，找到当前脉冲对应的跳数据长度
    
    t_240_chip = -(num_bits_pulse+th)*T/2:1/fs_IF:(num_bits_pulse+th)*T/2-1/fs_IF;
    t_240_pulse_temp = t_240_chip(1:end) + (T/oversamp_IF/2);  % 中心对称;
    t_240_pulse = t_240_pulse_temp((th_pat(pulse_idx)/2+3)*oversamp_IF+1:(th_pat(pulse_idx)/2+3+num_bits_pulse)*oversamp_IF);
    
    switch f_idx % 根据频点序号做对应的滤波
                        
        % 变换至零中频
        case 1  % 频点1
            rx_f1_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_1*t_240_pulse)); % 取差频
            rx_f1_128_temp1 = conv(rx_f1_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f1_128_temp2 = conv(rx_f1_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f1_128_1 = downsample(rx_f1_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128 
            rx_f1_128_2 = downsample(rx_f1_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f1_temp_1 = conv(rx_f1_128_1, LPF);
            rx_f1_temp_2 = conv(rx_f1_128_2, LPF);
            rx_f1 = [rx_f1_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f1_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f1;



        case 2  % 频点2
            rx_f2_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_2*t_240_pulse)); % 取差频
            rx_f2_128_temp1 = conv(rx_f2_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f2_128_temp2 = conv(rx_f2_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f2_128_1 = downsample(rx_f2_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128  
            rx_f2_128_2 = downsample(rx_f2_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f2_temp_1 = conv(rx_f2_128_1, LPF);
            rx_f2_temp_2 = conv(rx_f2_128_2, LPF);
            rx_f2 = [rx_f2_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f2_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f2;



        case 3  % 频点3
            rx_f3_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_3*t_240_pulse)); % 取差频
            rx_f3_128_temp1 = conv(rx_f3_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f3_128_temp2 = conv(rx_f3_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f3_128_1 = downsample(rx_f3_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128  
            rx_f3_128_2 = downsample(rx_f3_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);  
            rx_f3_temp_1 = conv(rx_f3_128_1, LPF);
            rx_f3_temp_2 = conv(rx_f3_128_2, LPF);
            rx_f3 = [rx_f3_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f3_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f3;



        case 4  % 频点4
            rx_f4_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_4*t_240_pulse)); % 取差频
            rx_f4_128_temp1 = conv(rx_f4_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f4_128_temp2 = conv(rx_f4_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f4_128_1 = downsample(rx_f4_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128  
            rx_f4_128_2 = downsample(rx_f4_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);  
            rx_f4_temp_1 = conv(rx_f4_128_1, LPF);
            rx_f4_temp_2 = conv(rx_f4_128_2, LPF);
            rx_f4 = [rx_f4_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f4_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f4;



        case 5  % 频点5
            rx_f5_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_5*t_240_pulse)); % 取差频
            rx_f5_128_temp1 = conv(rx_f5_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f5_128_temp2 = conv(rx_f5_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f5_128_1 = downsample(rx_f5_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128
            rx_f5_128_2 = downsample(rx_f5_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);  
            rx_f5_temp_1 = conv(rx_f5_128_1, LPF);
            rx_f5_temp_2 = conv(rx_f5_128_2, LPF);
            rx_f5 = [rx_f5_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f5_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f5;



        case 6  % 频点6                            
            rx_f6_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_1*t_240_pulse)); % 取差频
            rx_f6_128_temp1 = conv(rx_f6_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f6_128_temp2 = conv(rx_f6_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f6_128_1 = downsample(rx_f6_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128        
            rx_f6_128_2 = downsample(rx_f6_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f6_temp_1 = conv(rx_f6_128_1, LPF);
            rx_f6_temp_2 = conv(rx_f6_128_2, LPF);
            rx_f6 = [rx_f6_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f6_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f6;



        case 7  % 频点7
            rx_f7_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_2*t_240_pulse)); % 取差频
            rx_f7_128_temp1 = conv(rx_f7_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f7_128_temp2 = conv(rx_f7_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f7_128_1 = downsample(rx_f7_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128        
            rx_f7_128_2 = downsample(rx_f7_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f7_temp_1 = conv(rx_f7_128_1, LPF);
            rx_f7_temp_2 = conv(rx_f7_128_2, LPF);
            rx_f7 = [rx_f7_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f7_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f7;



        case 8  % 频点8    
            rx_f8_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_3*t_240_pulse)); % 取差频
            rx_f8_128_temp1 = conv(rx_f8_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f8_128_temp2 = conv(rx_f8_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f8_128_1 = downsample(rx_f8_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128        
            rx_f8_128_2 = downsample(rx_f8_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f8_temp_1 = conv(rx_f8_128_1, LPF);
            rx_f8_temp_2 = conv(rx_f8_128_2, LPF);
            rx_f8 = [rx_f8_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f8_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f8;
     


        case 9  % 频点9  
            rx_f9_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_4*t_240_pulse)); % 取差频
            rx_f9_128_temp1 = conv(rx_f9_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f9_128_temp2 = conv(rx_f9_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f9_128_1 = downsample(rx_f9_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128  
            rx_f9_128_2 = downsample(rx_f9_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f9_temp_1 = conv(rx_f9_128_1, LPF);
            rx_f9_temp_2 = conv(rx_f9_128_2, LPF);
            rx_f9 = [rx_f9_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f9_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f9;



        case 10 % 频点10   
            rx_f10_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_5*t_240_pulse)); % 取差频
            rx_f10_128_temp1 = conv(rx_f10_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f10_128_temp2 = conv(rx_f10_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f10_128_1 = downsample(rx_f10_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128   
            rx_f10_128_2 = downsample(rx_f10_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f10_temp_1 = conv(rx_f10_128_1, LPF);
            rx_f10_temp_2 = conv(rx_f10_128_2, LPF);
            rx_f10 = [rx_f10_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f10_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f10;



        case 11 % 频点11                           
            rx_f11_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_1*t_240_pulse)); % 取差频
            rx_f11_128_temp1 = conv(rx_f11_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f11_128_temp2 = conv(rx_f11_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f11_128_1 = downsample(rx_f11_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128   
            rx_f11_128_2 = downsample(rx_f11_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7); 
            rx_f11_temp_1 = conv(rx_f11_128_1, LPF);
            rx_f11_temp_2 = conv(rx_f11_128_2, LPF);
            rx_f11 = [rx_f11_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f11_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f11;



        case 12 % 频点12
            rx_f12_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_2*t_240_pulse)); % 取差频
            rx_f12_128_temp1 = conv(rx_f12_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f12_128_temp2 = conv(rx_f12_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f12_128_1 = downsample(rx_f12_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128
            rx_f12_128_2 = downsample(rx_f12_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f12_temp_1 = conv(rx_f12_128_1, LPF);
            rx_f12_temp_2 = conv(rx_f12_128_2, LPF);
            rx_f12 = [rx_f12_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f12_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f12;



        case 13 % 频点13
            rx_f13_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_3*t_240_pulse)); % 取差频
            rx_f13_128_temp1 = conv(rx_f13_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f13_128_temp2 = conv(rx_f13_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f13_128_1 = downsample(rx_f13_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128 
            rx_f13_128_2 = downsample(rx_f13_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f13_temp_1 = conv(rx_f13_128_1, LPF);
            rx_f13_temp_2 = conv(rx_f13_128_2, LPF);
            rx_f13 = [rx_f13_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f13_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f13;



        case 14 % 频点14
            rx_f14_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_4*t_240_pulse)); % 取差频
            rx_f14_128_temp1 = conv(rx_f14_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f14_128_temp2 = conv(rx_f14_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f14_128_1 = downsample(rx_f14_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128        
            rx_f14_128_2 = downsample(rx_f14_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7); 
            rx_f14_temp_1 = conv(rx_f14_128_1, LPF);
            rx_f14_temp_2 = conv(rx_f14_128_2, LPF);
            rx_f14 = [rx_f14_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f14_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f14;



        case 15 % 频点15
            rx_f15_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_1*t_240_pulse)); % 取差频
            rx_f15_128_temp1 = conv(rx_f15_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f15_128_temp2 = conv(rx_f15_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f15_128_1 = downsample(rx_f15_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128        
            rx_f15_128_2 = downsample(rx_f15_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f15_temp_1 = conv(rx_f15_128_1, LPF);
            rx_f15_temp_2 = conv(rx_f15_128_2, LPF);
            rx_f15 = [rx_f15_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f15_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f15;



        case 16 % 频点16
            rx_f16_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_2*t_240_pulse)); % 取差频
            rx_f16_128_temp1 = conv(rx_f16_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f16_128_temp2 = conv(rx_f16_temp1(280*oversamp_IF+1:end), LPF_2);    % 低通滤波(后24chip)
            rx_f16_128_1 = downsample(rx_f16_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128
            rx_f16_128_2 = downsample(rx_f16_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f16_temp_1 = conv(rx_f16_128_1, LPF);
            rx_f16_temp_2 = conv(rx_f16_128_2, LPF);
            rx_f16 = [rx_f16_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f16_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f16;



        case 17 % 频点17
            rx_f17_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_3*t_240_pulse)); % 取差频
            rx_f17_128_temp1 = conv(rx_f17_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f17_128_temp2 = conv(rx_f17_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f17_128_1 = downsample(rx_f17_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128        
            rx_f17_128_2 = downsample(rx_f17_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f17_temp_1 = conv(rx_f17_128_1, LPF);
            rx_f17_temp_2 = conv(rx_f17_128_2, LPF);
            rx_f17 = [rx_f17_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f17_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f17;
                            


        case 18 % 频点18
            rx_f18_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_4*t_240_pulse)); % 取差频
            rx_f18_128_temp1 = conv(rx_f18_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f18_128_temp2 = conv(rx_f18_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f18_128_1 = downsample(rx_f18_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128        
            rx_f18_128_2 = downsample(rx_f18_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f18_temp_1 = conv(rx_f18_128_1, LPF);
            rx_f18_temp_2 = conv(rx_f18_128_2, LPF);
            rx_f18 = [rx_f18_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f18_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f18;



        case 19 % 频点19
            rx_f19_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan5_1*t_240_pulse)); % 取差频
            rx_f19_128_temp1 = conv(rx_f19_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f19_128_temp2 = conv(rx_f19_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f19_128_1 = downsample(rx_f19_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128      
            rx_f19_128_2 = downsample(rx_f19_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f19_temp_1 = conv(rx_f19_128_1, LPF);
            rx_f19_temp_2 = conv(rx_f19_128_2, LPF);
            rx_f19 = [rx_f19_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f19_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f19;
                          


        case 20 % 频点20
            rx_f20_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan5_2*t_240_pulse)); % 取差频
            rx_f20_128_temp1 = conv(rx_f20_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f20_128_temp2 = conv(rx_f20_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f20_128_1 = downsample(rx_f20_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128
            rx_f20_128_2 = downsample(rx_f20_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f20_temp_1 = conv(rx_f20_128_1, LPF);
            rx_f20_temp_2 = conv(rx_f20_128_2, LPF);
            rx_f20 = [rx_f20_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f20_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f20;



        case 21 % 频点21
            rx_f21_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan5_3*t_240_pulse)); % 取差频
            rx_f21_128_temp1 = conv(rx_f21_temp1(1:280*oversamp_IF), LPF_2);  % 低通滤波(前280chip)
            rx_f21_128_temp2 = conv(rx_f21_temp1(280*oversamp_IF+1:end), LPF_2);  % 低通滤波(后24chip)
            rx_f21_128_1 = downsample(rx_f21_128_temp1(S_lpf2+1:S_lpf2+280*oversamp_IF), 8, 7);  % 8倍抽取， 采样率降为128
            rx_f21_128_2 = downsample(rx_f21_128_temp2(S_lpf2+1:S_lpf2+24*oversamp_IF), 8, 7);
            rx_f21_temp_1 = conv(rx_f21_128_1, LPF);
            rx_f21_temp_2 = conv(rx_f21_128_2, LPF);
            rx_f21 = [rx_f21_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_f21_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f21;
  
    end
                
    
end