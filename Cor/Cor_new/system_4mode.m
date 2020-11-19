clear all; clc;

% GMSK ������ܷ��� ��״̬��Viterbi������ܲ��� ��ʱ+Ƶƫ(0.001*Rb) 
%       ����ģʽ        4�֣�2M@A, 2M@B, 500K, 250K
%       ֡�ṹ          2M@A:      24bit��ͬ��ͷS1�� + 256bit�����ݣ� + 24bit��ͬ��ͷS2�� 
%                      2M@B:      0~Ln/2bit(����) + 24bit��ͬ��ͷS1�� + 256bit�����ݣ� + 24bit��ͬ��ͷS2�� + 0~Ln/2bit(����)
%                      500K\250K:  24bit��ͬ��ͷS1��+ 21bit��ͬ��ͷS3��+ 214bit�����ݣ�+ 21bit��ͬ��ͷS4�� + 24bit��ͬ��ͷS2��
%       ��������        2M@A\2M@B:   12           
%                      500K:  48
%                      250K:  96
%       һ֡����λ��     2M@A:   256*12 = 3072 chips
%                      2M@B:   256*12 + 512*6 = 6144 chips
%                      500K:    214*48 = 10272 chips������32bit 0��
%                      250K:    214*96 = 20544 chips������64bit 0��

% ����˵���� 
%       mode:  ==1: 2M@A;  
%              ==2: 2M@B;
%              ==3: 500K;
%              ==4: 250K
%

% ����������:  5֡���� 
%            2Mbps@A ģʽ: 5*3072 = 15360 bit
%            2Mbps@B ģʽ: 5*6144 = 30720 bit
%            500Kbps: 5*10272 = 51360 bit
%            250Kbps: 5*20544 = 102720 bit
%

mode = 1;  %����ʱ���Ը�

frame_num = 5;  % �������֡����
num_bits_pn = 24;  % ͬ��ͷS1\S2����
num_bits_pn_2 = 21;  % ͬ��ͷS3\S4����

bit_rate = 16e6;  % �������� 
T = 1/bit_rate;  % ����ʱ��
num_bits_pulse = 304;  % 2M_A\500K\250K���峤�ȣ�����ǰ��ͬ��ͷ��, 2M_B ����ȥ��ǰ�������ݲ��ֺ�ĳ���

fs_IF = 1024e6;  % ��Ƶ����Ƶ�źŲ�������
fs_BB = 128e6;  % �����źŲ�������
oversamp_BB = T * fs_BB;  % �����źŹ���������
oversamp_IF = T * fs_IF;  % ��Ƶ����Ƶ�źŹ���������
T_s_BB = 1/fs_BB;  % �����������
T_s_IF = 1/fs_IF;  % ��Ƶ����Ƶ�������
BER = zeros(1,15); %�������
Eb_N0 = 1: 1: 14;

% �����Ѵ�����
load('lib/g_1024.mat');  % GMSK���� g����(��˹��ͨ�˲���)
load('lib/f_trans.mat');  % 21��Ƶ�㣨��PN���еã�
% ���ն��˲���
load('lib/filter/BPF_CHAN1_UD.mat');  % ��ͨ ͨ��1
load('lib/filter/BPF_CHAN2_UD.mat');  % ��ͨ ͨ��2
load('lib/filter/BPF_CHAN3_UD.mat');  % ��ͨ ͨ��3
load('lib/filter/BPF_CHAN4_UD.mat');  % ��ͨ ͨ��4
load('lib/filter/BPF_CHAN5_UD.mat');  % ��ͨ ͨ��5
load('lib/filter/BPF_CHAN12_2_UD.mat');  % ��ͨ ͨ��1��2 �ڶ����˲�
load('lib/filter/BPF_CHAN34_2_UD.mat');  % ��ͨ ͨ��3��4 �ڶ����˲�
load('lib/filter/BPF_CHAN5_2_UD.mat');  % ��ͨ ͨ��5 �ڶ����˲�
load('lib/filter/c0_f_128.mat');  % ���ƥ���˲���1
load('lib/filter/c1_f_128.mat');  % ���ƥ���˲���2
load('lib/filter/LPF.mat');  % ��ͨ�˲��� 61�� ͨ��: 4.8MHz ����Ƶ��128MHz
load('lib/filter/LPF_2.mat');  % ��ͨ�˲��� 507�� ͨ��: 5MHz  ����Ƶ��1024MHz
S_lpf = 30;
S_lpf2 = 127;
S_bpf = 253;

% ÿ֡�ֳɵ���������� ��������ģʽѡ��
if mode == 1
    num_pulses = 12;
elseif mode == 2
    num_pulses = 12;
elseif mode == 3
    num_pulses = 48;
else 
    num_pulses = 96;
end

% �������������� (һ��frame_num֡)
if mode == 2
    mat_row = frame_num;  %����ģʽ2������֡��=5
else
    mat_row = num_pulses * frame_num;  
end

% ����PN�⡢��Ƶͼ������ʱͼ�����������ݾ��� 
bits = data_gen(mat_row, mode);   % ˫�����루��������ͬ��ͷ��ֻ��Bû�У�
[th_pat_lib_1, fh_pat_lib_1] = TF_gen;   % ��Ƶ����ʱ��ͼ��1
[th_pat_lib_2, fh_pat_lib_2] = TF_gen;   % ��Ƶ����ʱ��ͼ��2

[pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1] = pn_gen;  % 0\1 ��
[pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2] = pn_gen;

% ͬ��ͷ����
% mode 1
[wav_S1_1_mode1, wav_S2_1_mode1] = wav_gen(pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1, 1);
[wav_S1_2_mode1, wav_S2_2_mode1] = wav_gen(pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2, 1);
% mode 2
[wav_S1_1_mode2, wav_S2_1_mode2] = wav_gen(pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1, 2);
[wav_S1_2_mode2, wav_S2_2_mode2] = wav_gen(pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2, 2);
% mode 3
[wav_S1_S3_1_mode3, wav_S4_S2_1_mode3] = wav_gen(pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1, 3);
[wav_S1_S3_2_mode3, wav_S4_S2_2_mode3] = wav_gen(pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2, 3);
% mode 4
[wav_S1_S3_1_mode4, wav_S4_S2_1_mode4] = wav_gen(pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1, 4);
[wav_S1_S3_2_mode4, wav_S4_S2_2_mode4] = wav_gen(pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2, 4); 

% ͬ��ͷ��ͨ�˲�����
% mode 1
[wav_S1_1_LPF_mode1, wav_S2_1_LPF_mode1, wav_S1_2_LPF_mode1, wav_S2_2_LPF_mode1] = wavLPF(wav_S1_1_mode1, wav_S2_1_mode1, wav_S1_2_mode1, wav_S2_2_mode1, LPF, S_lpf);
% mode 2
[wav_S1_1_LPF_mode2, wav_S2_1_LPF_mode2, wav_S1_2_LPF_mode2, wav_S2_2_LPF_mode2] = wavLPF(wav_S1_1_mode2, wav_S2_1_mode2, wav_S1_2_mode2, wav_S2_2_mode2, LPF, S_lpf);
% mode 3
[wav_S1_S3_1_LPF_mode3, wav_S4_S2_1_LPF_mode3, wav_S1_S3_2_LPF_mode3, wav_S4_S2_2_LPF_mode3] = wavLPF(wav_S1_S3_1_mode3, wav_S4_S2_1_mode3, wav_S1_S3_2_mode3, wav_S4_S2_2_mode3, LPF, S_lpf);
% mode 4
[wav_S1_S3_1_LPF_mode4, wav_S4_S2_1_LPF_mode4, wav_S1_S3_2_LPF_mode4, wav_S4_S2_2_LPF_mode4] = wavLPF(wav_S1_S3_1_mode4, wav_S4_S2_1_mode4, wav_S1_S3_2_mode4, wav_S4_S2_2_mode4, LPF, S_lpf);

% ������ �����ڽ����ͳ�������ʣ�
if mode == 1
    for i = 1:mat_row/num_pulses  %֡��
        for j = 1:num_pulses
            bit_frame((j-1)*256+1:j*256) = bits((i-1)*num_pulses+j,25:280);  %�൱����������λ���������г�һ��
        end
        bit_check(i,:) = bit_frame;
    end
elseif mode == 2
    bit_check = bits;   
elseif mode == 3
    for i = 1:mat_row/num_pulses
        for j = 1:num_pulses
            bit_frame((j-1)*214+1:j*214) = bits((i-1)*num_pulses+j,46:259);
        end
        bit_check(i,:) = bit_frame;
    end
else
    for i = 1:mat_row/num_pulses
        for j = 1:num_pulses
            bit_frame((j-1)*214+1:j*214) = bits((i-1)*num_pulses+j,46:259);
        end
        bit_check(i,:) = bit_frame;
    end
end


% �����0:17 ѭ��
for SNR_idx = 10:10 %1:length(Eb_N0)

    N0 = 1/10^(Eb_N0(SNR_idx)/10);  % ���㵱ǰ����ȶ�Ӧ�ĸ�˹�������������ܶ�ֵ N0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                 ���䲿��
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % ����GMSK����
    signal_trans = transmitter(bits, fh_pat_lib_1, th_pat_lib_1, fh_pat_lib_2, th_pat_lib_2, pn_lib_S1_1, pn_lib_S1_2,...
    pn_lib_S2_1, pn_lib_S2_2, pn_lib_S3_1, pn_lib_S3_2, pn_lib_S4_1, pn_lib_S4_2, mode);

    % ǰ����
    signal_trans_2 = [zeros(1,0.5*oversamp_IF), signal_trans, zeros(1,100000*oversamp_IF)];  % ��0.5�����ŵĿհ�
    
    % �Ӹ�˹������
    White_N = sqrt(N0*oversamp_IF/2)*(randn(1, length(signal_trans_2))+1i*randn(1, length(signal_trans_2)));
    sig_N = signal_trans_2 + White_N; 
    
    % ��Ƶƫ (0.016MHz)(ǧ��֮һ����������)
    t = 0:1/oversamp_IF:length(signal_trans_2)/oversamp_IF-1/oversamp_IF;
    t = t(1:end) + (T/oversamp_IF/2);  % ���ĶԳ�
    pha = 0.016e6*2*pi*t*T; 
    pha = mod(pha, 2*pi);
    pha_cplx = complex(cos(pha), -sin(pha));
    rx = sig_N .* pha_cplx;  % rx: Ϊ������ģ��Ĵ���Ƶƫ����Ƶ�ź�

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               ���ղ���
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    result = receiver(rx, fh_pat_lib_1, th_pat_lib_1, fh_pat_lib_2, th_pat_lib_2, ...
                        wav_S1_1_LPF_mode1, wav_S1_2_LPF_mode1, wav_S2_1_LPF_mode1, wav_S2_2_LPF_mode1, ...
                        wav_S1_1_LPF_mode2, wav_S1_2_LPF_mode2, wav_S2_1_LPF_mode2, wav_S2_2_LPF_mode2, ...
                        wav_S1_S3_1_LPF_mode3, wav_S4_S2_1_LPF_mode3, wav_S1_S3_2_LPF_mode3, wav_S4_S2_2_LPF_mode3, ...
                        wav_S1_S3_1_LPF_mode4, wav_S4_S2_1_LPF_mode4, wav_S1_S3_2_LPF_mode4, wav_S4_S2_2_LPF_mode4, ...
                        frame_num);

      
    [a,b] = size(bit_check); 
    error = sum(sum(abs(bit_check-result)/2 == 1));
    miss = sum(sum(abs(bit_check-result)/2 == 0.5));
    total = a*b-miss;
    pe(SNR_idx)=error/(total);
    
end




