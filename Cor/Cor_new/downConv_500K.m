function rx_pulse_mat = downConv_500K(wav, num_pulses, fh_pat)

% ������Ƶͼ����һ֡�а���������������������յ���Ƶ�źŲ����±�Ƶ����������64���������ʽ�Ϊ8����������
% ����500Kbpsģʽ

% ��������
bit_rate = 16e6;  % ��������
T = 1/bit_rate;  % ����ʱ��
fs_IF = 1024e6;  % ��Ƶ�źŲ�������
fs_BB = 128e6;  % �����źŲ�������
oversamp_BB = T * fs_BB;  % ��������Ϊ 8����������
oversamp_IF = T * fs_IF;
num_bits_pulse = 304;

% ��ͨ�˲���
load('lib/filter/LPF_2.mat'); % ������1024MHz��ͨ��8MHz
load('lib/filter/LPF.mat'); % ������128MHz��ͨ��5MHz
S_lpf = 30;
S_lpf2 = 127;

% 21��Ƶ���±�Ƶ������Ƶ�������Ƶ��
f_chan12_1 = 213.34e6;  % ͨ��1��2 ��1��Ƶ���Ӧ��Ƶ��
f_chan12_2 = 226.67e6;  % ͨ��1��2 ��2��Ƶ���Ӧ��Ƶ��
f_chan12_3 = 240e6;     % ͨ��1��2 ��3��Ƶ���Ӧ��Ƶ��
f_chan12_4 = 253.33e6;  % ͨ��1��2 ��4��Ƶ���Ӧ��Ƶ��
f_chan12_5 = 266.66e6;  % ͨ��1��2 ��5��Ƶ���Ӧ��Ƶ��
f_chan34_1 = 220.005e6; % ͨ��3��4 ��1��Ƶ���Ӧ��Ƶ��
f_chan34_2 = 233.335e6; % ͨ��3��4 ��2��Ƶ���Ӧ��Ƶ��
f_chan34_3 = 246.665e6; % ͨ��3��4 ��3��Ƶ���Ӧ��Ƶ��
f_chan34_4 = 259.995e6; % ͨ��3��4 ��4��Ƶ���Ӧ��Ƶ��
f_chan5_1 = f_chan12_2; % ͨ��5 ��1��Ƶ���Ӧ��Ƶ��
f_chan5_2 = f_chan12_3; % ͨ��5 ��2��Ƶ���Ӧ��Ƶ��
f_chan5_3 = f_chan12_4; % ͨ��5 ��3��Ƶ���Ӧ��Ƶ��

t_240_chip = -num_bits_pulse*T/2:1/fs_IF:num_bits_pulse*T/2-1/fs_IF;
t_240_pulse = t_240_chip(1:end) + (T/oversamp_IF/2);  % ���ĶԳ�;

rx_pulse_mat = zeros(num_pulses, num_bits_pulse * oversamp_BB);
for pulse_idx = 1:num_pulses
    
    f_idx = fh_pat(pulse_idx);  % ��������ź���Ƶͼ�����ҵ���ǰ�����Ӧ��Ƶ��
    
    switch f_idx % ����Ƶ���������Ӧ���˲�
                        
        % �任������Ƶ
        case 1  % Ƶ��1
            rx_f1_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_1*t_240_pulse)); % ȡ��Ƶ
            rx_f1_128_temp1 = conv(rx_f1_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f1_128_temp2 = conv(rx_f1_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f1_128_1 = downsample(rx_f1_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128 
            rx_f1_128_2 = downsample(rx_f1_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f1_temp_1 = conv(rx_f1_128_1, LPF);
            rx_f1_temp_2 = conv(rx_f1_128_2, LPF);
            rx_f1 = [rx_f1_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f1_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f1;



        case 2  % Ƶ��2
            rx_f2_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_2*t_240_pulse)); % ȡ��Ƶ
            rx_f2_128_temp1 = conv(rx_f2_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f2_128_temp2 = conv(rx_f2_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f2_128_1 = downsample(rx_f2_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128  
            rx_f2_128_2 = downsample(rx_f2_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f2_temp_1 = conv(rx_f2_128_1, LPF);
            rx_f2_temp_2 = conv(rx_f2_128_2, LPF);
            rx_f2 = [rx_f2_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f2_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f2;



        case 3  % Ƶ��3
            rx_f3_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_3*t_240_pulse)); % ȡ��Ƶ
            rx_f3_128_temp1 = conv(rx_f3_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f3_128_temp2 = conv(rx_f3_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f3_128_1 = downsample(rx_f3_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128  
            rx_f3_128_2 = downsample(rx_f3_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);  
            rx_f3_temp_1 = conv(rx_f3_128_1, LPF);
            rx_f3_temp_2 = conv(rx_f3_128_2, LPF);
            rx_f3 = [rx_f3_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f3_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f3;



        case 4  % Ƶ��4
            rx_f4_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_4*t_240_pulse)); % ȡ��Ƶ
            rx_f4_128_temp1 = conv(rx_f4_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f4_128_temp2 = conv(rx_f4_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f4_128_1 = downsample(rx_f4_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128  
            rx_f4_128_2 = downsample(rx_f4_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);  
            rx_f4_temp_1 = conv(rx_f4_128_1, LPF);
            rx_f4_temp_2 = conv(rx_f4_128_2, LPF);
            rx_f4 = [rx_f4_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f4_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f4;



        case 5  % Ƶ��5
            rx_f5_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_5*t_240_pulse)); % ȡ��Ƶ
            rx_f5_128_temp1 = conv(rx_f5_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f5_128_temp2 = conv(rx_f5_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f5_128_1 = downsample(rx_f5_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128
            rx_f5_128_2 = downsample(rx_f5_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);  
            rx_f5_temp_1 = conv(rx_f5_128_1, LPF);
            rx_f5_temp_2 = conv(rx_f5_128_2, LPF);
            rx_f5 = [rx_f5_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f5_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f5;



        case 6  % Ƶ��6                            
            rx_f6_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_1*t_240_pulse)); % ȡ��Ƶ
            rx_f6_128_temp1 = conv(rx_f6_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f6_128_temp2 = conv(rx_f6_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f6_128_1 = downsample(rx_f6_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128        
            rx_f6_128_2 = downsample(rx_f6_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f6_temp_1 = conv(rx_f6_128_1, LPF);
            rx_f6_temp_2 = conv(rx_f6_128_2, LPF);
            rx_f6 = [rx_f6_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f6_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f6;



        case 7  % Ƶ��7
            rx_f7_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_2*t_240_pulse)); % ȡ��Ƶ
            rx_f7_128_temp1 = conv(rx_f7_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f7_128_temp2 = conv(rx_f7_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f7_128_1 = downsample(rx_f7_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128        
            rx_f7_128_2 = downsample(rx_f7_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f7_temp_1 = conv(rx_f7_128_1, LPF);
            rx_f7_temp_2 = conv(rx_f7_128_2, LPF);
            rx_f7 = [rx_f7_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f7_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f7;



        case 8  % Ƶ��8    
            rx_f8_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_3*t_240_pulse)); % ȡ��Ƶ
            rx_f8_128_temp1 = conv(rx_f8_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f8_128_temp2 = conv(rx_f8_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f8_128_1 = downsample(rx_f8_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128        
            rx_f8_128_2 = downsample(rx_f8_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f8_temp_1 = conv(rx_f8_128_1, LPF);
            rx_f8_temp_2 = conv(rx_f8_128_2, LPF);
            rx_f8 = [rx_f8_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f8_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f8;
     


        case 9  % Ƶ��9  
            rx_f9_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_4*t_240_pulse)); % ȡ��Ƶ
            rx_f9_128_temp1 = conv(rx_f9_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f9_128_temp2 = conv(rx_f9_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f9_128_1 = downsample(rx_f9_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128  
            rx_f9_128_2 = downsample(rx_f9_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f9_temp_1 = conv(rx_f9_128_1, LPF);
            rx_f9_temp_2 = conv(rx_f9_128_2, LPF);
            rx_f9 = [rx_f9_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f9_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f9;



        case 10 % Ƶ��10   
            rx_f10_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_5*t_240_pulse)); % ȡ��Ƶ
            rx_f10_128_temp1 = conv(rx_f10_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f10_128_temp2 = conv(rx_f10_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f10_128_1 = downsample(rx_f10_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128   
            rx_f10_128_2 = downsample(rx_f10_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f10_temp_1 = conv(rx_f10_128_1, LPF);
            rx_f10_temp_2 = conv(rx_f10_128_2, LPF);
            rx_f10 = [rx_f10_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f10_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f10;



        case 11 % Ƶ��11                           
            rx_f11_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_1*t_240_pulse)); % ȡ��Ƶ
            rx_f11_128_temp1 = conv(rx_f11_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f11_128_temp2 = conv(rx_f11_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f11_128_1 = downsample(rx_f11_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128   
            rx_f11_128_2 = downsample(rx_f11_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7); 
            rx_f11_temp_1 = conv(rx_f11_128_1, LPF);
            rx_f11_temp_2 = conv(rx_f11_128_2, LPF);
            rx_f11 = [rx_f11_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f11_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f11;



        case 12 % Ƶ��12
            rx_f12_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_2*t_240_pulse)); % ȡ��Ƶ
            rx_f12_128_temp1 = conv(rx_f12_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f12_128_temp2 = conv(rx_f12_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f12_128_1 = downsample(rx_f12_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128
            rx_f12_128_2 = downsample(rx_f12_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f12_temp_1 = conv(rx_f12_128_1, LPF);
            rx_f12_temp_2 = conv(rx_f12_128_2, LPF);
            rx_f12 = [rx_f12_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f12_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f12;



        case 13 % Ƶ��13
            rx_f13_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_3*t_240_pulse)); % ȡ��Ƶ
            rx_f13_128_temp1 = conv(rx_f13_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f13_128_temp2 = conv(rx_f13_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f13_128_1 = downsample(rx_f13_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128 
            rx_f13_128_2 = downsample(rx_f13_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f13_temp_1 = conv(rx_f13_128_1, LPF);
            rx_f13_temp_2 = conv(rx_f13_128_2, LPF);
            rx_f13 = [rx_f13_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f13_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f13;



        case 14 % Ƶ��14
            rx_f14_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_4*t_240_pulse)); % ȡ��Ƶ
            rx_f14_128_temp1 = conv(rx_f14_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f14_128_temp2 = conv(rx_f14_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f14_128_1 = downsample(rx_f14_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128        
            rx_f14_128_2 = downsample(rx_f14_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7); 
            rx_f14_temp_1 = conv(rx_f14_128_1, LPF);
            rx_f14_temp_2 = conv(rx_f14_128_2, LPF);
            rx_f14 = [rx_f14_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f14_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f14;



        case 15 % Ƶ��15
            rx_f15_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_1*t_240_pulse)); % ȡ��Ƶ
            rx_f15_128_temp1 = conv(rx_f15_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f15_128_temp2 = conv(rx_f15_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f15_128_1 = downsample(rx_f15_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128        
            rx_f15_128_2 = downsample(rx_f15_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f15_temp_1 = conv(rx_f15_128_1, LPF);
            rx_f15_temp_2 = conv(rx_f15_128_2, LPF);
            rx_f15 = [rx_f15_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f15_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f15;



        case 16 % Ƶ��16
            rx_f16_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_2*t_240_pulse)); % ȡ��Ƶ
            rx_f16_128_temp1 = conv(rx_f16_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f16_128_temp2 = conv(rx_f16_temp1(259*oversamp_IF+1:end), LPF_2);    % ��ͨ�˲�(��45chip)
            rx_f16_128_1 = downsample(rx_f16_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128
            rx_f16_128_2 = downsample(rx_f16_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f16_temp_1 = conv(rx_f16_128_1, LPF);
            rx_f16_temp_2 = conv(rx_f16_128_2, LPF);
            rx_f16 = [rx_f16_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f16_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f16;



        case 17 % Ƶ��17
            rx_f17_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_3*t_240_pulse)); % ȡ��Ƶ
            rx_f17_128_temp1 = conv(rx_f17_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f17_128_temp2 = conv(rx_f17_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f17_128_1 = downsample(rx_f17_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128        
            rx_f17_128_2 = downsample(rx_f17_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f17_temp_1 = conv(rx_f17_128_1, LPF);
            rx_f17_temp_2 = conv(rx_f17_128_2, LPF);
            rx_f17 = [rx_f17_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f17_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f17;
                            


        case 18 % Ƶ��18
            rx_f18_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_4*t_240_pulse)); % ȡ��Ƶ
            rx_f18_128_temp1 = conv(rx_f18_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f18_128_temp2 = conv(rx_f18_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f18_128_1 = downsample(rx_f18_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128        
            rx_f18_128_2 = downsample(rx_f18_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f18_temp_1 = conv(rx_f18_128_1, LPF);
            rx_f18_temp_2 = conv(rx_f18_128_2, LPF);
            rx_f18 = [rx_f18_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f18_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f18;



        case 19 % Ƶ��19
            rx_f19_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan5_1*t_240_pulse)); % ȡ��Ƶ
            rx_f19_128_temp1 = conv(rx_f19_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f19_128_temp2 = conv(rx_f19_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f19_128_1 = downsample(rx_f19_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128      
            rx_f19_128_2 = downsample(rx_f19_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f19_temp_1 = conv(rx_f19_128_1, LPF);
            rx_f19_temp_2 = conv(rx_f19_128_2, LPF);
            rx_f19 = [rx_f19_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f19_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f19;
                          


        case 20 % Ƶ��20
            rx_f20_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan5_2*t_240_pulse)); % ȡ��Ƶ
            rx_f20_128_temp1 = conv(rx_f20_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f20_128_temp2 = conv(rx_f20_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f20_128_1 = downsample(rx_f20_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128
            rx_f20_128_2 = downsample(rx_f20_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f20_temp_1 = conv(rx_f20_128_1, LPF);
            rx_f20_temp_2 = conv(rx_f20_128_2, LPF);
            rx_f20 = [rx_f20_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f20_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f20;



        case 21 % Ƶ��21
            rx_f21_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan5_3*t_240_pulse)); % ȡ��Ƶ
            rx_f21_128_temp1 = conv(rx_f21_temp1(1:259*oversamp_IF), LPF_2);  % ��ͨ�˲�(ǰ259chip)
            rx_f21_128_temp2 = conv(rx_f21_temp1(259*oversamp_IF+1:end), LPF_2);  % ��ͨ�˲�(��45chip)
            rx_f21_128_1 = downsample(rx_f21_128_temp1(S_lpf2+1:S_lpf2+259*oversamp_IF), 8, 7);  % 8����ȡ�� �����ʽ�Ϊ128
            rx_f21_128_2 = downsample(rx_f21_128_temp2(S_lpf2+1:S_lpf2+45*oversamp_IF), 8, 7);
            rx_f21_temp_1 = conv(rx_f21_128_1, LPF);
            rx_f21_temp_2 = conv(rx_f21_128_2, LPF);
            rx_f21 = [rx_f21_temp_1(S_lpf+1:S_lpf+259*oversamp_BB), rx_f21_temp_2(S_lpf+1:S_lpf+45*oversamp_BB)];
            rx_pulse_mat(pulse_idx,:) = rx_f21;
  
    end
                
    
end