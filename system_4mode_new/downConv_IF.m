function rx_chan_240 = downConv_IF(rx, BPF_CHAN_1, BPF_CHAN_2, f_chan_240, indi)

% ��Ƶ�ź��±�Ƶ����Ƶ
% indi ����ָʾ��ǰ�����Ƿ�����ͨ��5 ��ʵ��ʵ�ֿ��ܲ���Ҫ��

% ��������
bit_rate = 16e6;  % ��������
T = 1/bit_rate;  % ����ʱ��
fs_IF = 1024e6;  % ��Ƶ����Ƶ�źŲ�������
oversamp_IF = T * fs_IF;  % ��Ƶ����Ƶ�źŹ���������
t_240_chip = -length(rx)/oversamp_IF*T/2:1/fs_IF:length(rx)/oversamp_IF*T/2-1/fs_IF;
t_240 = t_240_chip(1:end) + (T/oversamp_IF/2);  % ���ĶԳ�;

S_bpf = 253;

rx_chan_temp = conv(rx, BPF_CHAN_1);  %��һ�δ�ͨ�˲�
rx_chan = rx_chan_temp(S_bpf+1:S_bpf+length(rx));
% �任����Ƶ240MHz�����ٴ�ͨ�˲���
% ͨ��1��4(indi==0): ȡ��Ƶ        
% ͨ��5(indi==1): ȡ��Ƶ
if (~indi) 
    rx_chan_240_temp = conv(rx_chan .* exp(1i*(2*pi*f_chan_240*t_240)), BPF_CHAN_2);
else
    rx_chan_240_temp = conv(rx_chan .* exp(1i*(-2*pi*f_chan_240*t_240)), BPF_CHAN_2);
end
    
rx_chan_240 = rx_chan_240_temp(S_bpf+1:S_bpf+length(rx_chan));
