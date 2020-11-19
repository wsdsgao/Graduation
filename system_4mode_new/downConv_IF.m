function rx_chan_240 = downConv_IF(rx, BPF_CHAN_1, BPF_CHAN_2, f_chan_240, indi)

% 射频信号下变频至中频
% indi 用于指示当前波形是否属于通道5 （实际实现可能不需要）

% 参数定义
bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间
fs_IF = 1024e6;  % 射频、中频信号采样速率
oversamp_IF = T * fs_IF;  % 射频、中频信号过采样倍数
t_240_chip = -length(rx)/oversamp_IF*T/2:1/fs_IF:length(rx)/oversamp_IF*T/2-1/fs_IF;
t_240 = t_240_chip(1:end) + (T/oversamp_IF/2);  % 中心对称;

S_bpf = 253;

rx_chan_temp = conv(rx, BPF_CHAN_1);  %第一次带通滤波
rx_chan = rx_chan_temp(S_bpf+1:S_bpf+length(rx));
% 变换至中频240MHz（，再带通滤波）
% 通道1～4(indi==0): 取和频        
% 通道5(indi==1): 取差频
if (~indi) 
    rx_chan_240_temp = conv(rx_chan .* exp(1i*(2*pi*f_chan_240*t_240)), BPF_CHAN_2);
else
    rx_chan_240_temp = conv(rx_chan .* exp(1i*(-2*pi*f_chan_240*t_240)), BPF_CHAN_2);
end
    
rx_chan_240 = rx_chan_240_temp(S_bpf+1:S_bpf+length(rx_chan));
