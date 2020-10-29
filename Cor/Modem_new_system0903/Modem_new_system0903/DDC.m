function rx_pulse_128_LPF = DDC(wav)

% 根据跳频图案和一帧中包含的总脉冲个数，将接收的中频信号波形下变频到基带并从64倍符号速率降为8倍符号速率
% 用于2Mbps A模式

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
load('LPF_1024.mat'); % 采样率1024MHz，通带15MHz,阻带60MHz
load('LPF_128.mat'); % 采样率128MHz，通带5MHz，阻带8MHz

S_lpf = 30;
S_lpf2 = 127;

% 21个频点下变频至零中频所需各点频率
freq_21=[240-40/3*10:40/3:240+40/3*10];%定义调试5个频点
freq_5=[240-80/3,240-40/3,240,240+40/3,240+80/3];%定义调试5个频点
freq_1=[240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240,240];%定义调试1个频点
f_trans=freq_21*1000000;
% f_trans=freq_1*1000000;

% rx_pulse_BB = zeros(21, length(wav));

for pulse_idx = 1:21
    
    dt=(0:length(wav)-1);
    rx_pulse_BB(pulse_idx,:)=wav.*exp(-j*(f_trans(pulse_idx)/1024/1000000)*dt*2*pi);
    rx_pulse_LPF(pulse_idx,:)=conv(rx_pulse_BB(pulse_idx,:),LPF_1024);
    rx_pulse_128(pulse_idx,:)=resample(rx_pulse_LPF(pulse_idx,:),128,1024);
    rx_pulse_128_LPF(pulse_idx,:)=conv(rx_pulse_128(pulse_idx,:),LPF_128);
    
%     figure;
%     plot(real(rx_pulse_128_LPF(pulse_idx,:)),'b');
%     hold on;
%     plot(imag(rx_pulse_128_LPF(pulse_idx,:)),'r');
%     close;
    
end