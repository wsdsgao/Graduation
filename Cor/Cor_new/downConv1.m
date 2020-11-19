function rx_pulse_mat = downConv(wav, num_pulses, fh_pat)

% ������Ƶͼ����һ֡�а���������������������յ���Ƶ�źŲ����±�Ƶ����������64���������ʽ�Ϊ4����������
% ����2Mbps Aģʽ

% ��������
bit_rate = 16e6;  % ��������
T = 1/bit_rate;  % ����ʱ��
fs_IF = 1024e6;  % ��Ƶ�źŲ�������
fs_BB = 64e6;  % �����źŲ�������
oversamp_BB = T * fs_BB;  % �����źŹ���������
oversamp_IF = T * fs_IF;  % ��Ƶ�źŹ���������
num_bits_pulse = 304;  % 2Mbps A\500Kbps\250Kbps һ��������ĳ���
                       % 2Mbps B ÿ������ȥ��ǰ�������ݲ��ֺ�ĳ���

% ��ͨ�˲���
load('filter/LPF_1.mat'); % ������1024MHz��ͨ��20MHz
load('filter/LPF_2.mat'); % ������256MHz��ͨ��10MHz
S_lpf1 = 32;
S_lpf2 = 32;

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
            rx_f1_256 = downsample(rx_f1_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 2  % Ƶ��2
            rx_f2_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_2*t_240_pulse)); % ȡ��Ƶ
            rx_f2_256 = downsample(rx_f2_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f2_64 = downsample(rx_f2_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f2_64;



        case 3  % Ƶ��3
            rx_f3_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_3*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f3_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 4  % Ƶ��4
            rx_f4_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_4*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f4_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 5  % Ƶ��5
            rx_f5_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_5*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f5_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 6  % Ƶ��6                            
            rx_f6_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_1*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f6_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 7  % Ƶ��7
            rx_f7_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_2*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f7_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 8  % Ƶ��8    
            rx_f8_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_3*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f8_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;
     


        case 9  % Ƶ��9  
            rx_f9_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_4*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f9_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 10 % Ƶ��10   
            rx_f10_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan12_5*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f10_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 11 % Ƶ��11                           
            rx_f11_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_1*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f11_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 12 % Ƶ��12
            rx_f12_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_2*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f12_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 13 % Ƶ��13
            rx_f13_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_3*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f13_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 14 % Ƶ��14
            rx_f14_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_4*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f14_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 15 % Ƶ��15
            rx_f15_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_1*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f15_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 16 % Ƶ��16
            rx_f16_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_2*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f16_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 17 % Ƶ��17
            rx_f17_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_3*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f17_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;
                            


        case 18 % Ƶ��18
            rx_f18_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan34_4*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f18_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 19 % Ƶ��19
            rx_f19_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan5_1*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f19_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;
                          


        case 20 % Ƶ��20
            rx_f20_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan5_2*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f20_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;



        case 21 % Ƶ��21
            rx_f21_temp1 = wav(pulse_idx,:) .* exp(1i*(-2*pi*f_chan5_3*t_240_pulse)); % ȡ��Ƶ
            rx_f1_256 = downsample(rx_f21_temp1(1:304*oversamp_IF), 4);  % 4����ȡ�� �����ʽ�Ϊ256��3��ƫ������
            rx_f1_64 = downsample(rx_f1_256(1:304*oversamp_IF/4), 4);  % 4����ȡ�� �����ʽ�Ϊ64��3��ƫ������
            rx_pulse_mat(pulse_idx,:) = rx_f1_64;
  
    end
                
    
end