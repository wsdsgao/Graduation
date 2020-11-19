function rx_pulse_mat = downConv_BB(wav, num_pulses)

% ��������
bit_rate = 16e6;  % ��������
T = 1/bit_rate;  % ����ʱ��
fs_BB = 128e6;  % �����źŲ�������
oversamp_BB = T * fs_BB;  % ��������Ϊ 8����������
num_bits_pulse = 304;

% ��ͨ�˲���
load('lib/filter/LPF.mat'); % ������128MHz��ͨ��5MHz
S_lpf = 30;

rx_pulse_mat = zeros(num_pulses, num_bits_pulse * oversamp_BB);
for pulse_idx = 1:num_pulses
    
    rx_temp_1 = conv(wav(pulse_idx, 1:280*oversamp_BB), LPF);
    rx_temp_2 = conv(wav(pulse_idx, 280*oversamp_BB+1:end), LPF);

    rx = [rx_temp_1(S_lpf+1:S_lpf+280*oversamp_BB), rx_temp_2(S_lpf+1:S_lpf+24*oversamp_BB)];
    rx_pulse_mat(pulse_idx,:) = rx;    
    
end
