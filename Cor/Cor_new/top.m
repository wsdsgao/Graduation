clear all;
close all;

% 基本参数定义
bit_rate = 16e6; % 符号速率
T = 1/bit_rate;  % 符号时间
fs_IF = 1024e6;  % 中频信号采样速率
num_frame = 1;   % �??发�?�的帧数�??
num_data_frame = 1024; % �??帧有效数据长�??
oversamp_IF = T * fs_IF; % 射频过采样率

% 生成PN �?? 跳频跳时�??
% To Do: 每过�??段时间进行更�??
[th_pat_lib, fh_pat_lib] = TF_gen;   % 跳频、跳时�?�图�??
[pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4] = pn_gen;  % 0\1 �??

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                 发射部分
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 模式1
% 生成数据帧并进行LDPC编码
[data1, ~] = data_gen(num_frame, 1);
% 产生射频发生波形, 采样�??1024MHz
signal_trans_1 = transmitter(data1, num_frame, fh_pat_lib, th_pat_lib, pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 1);

% plot(abs(signal_trans_1));

% 模式2
[data2, ~] = data_gen(num_frame, 2);
signal_trans_2 = transmitter(data2, num_frame, fh_pat_lib, th_pat_lib, pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 2);

% 模式3
[data3, ~] = data_gen(num_frame, 3);
signal_trans_3 = transmitter(data3, num_frame, fh_pat_lib, th_pat_lib, pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 3);

% 模式4
[data4, ~] = data_gen(num_frame, 4);
signal_trans_4 = transmitter(data4, num_frame, fh_pat_lib, th_pat_lib, pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 4);

% 随机延时后相�??
% To Do: 给不同信号不同能�??
t_delay = round(rand(4, 1) * 1e-3 * fs_IF);

% signal_trans_D_1 = [zeros(1, t_delay(1)), signal_trans_1];
signal_trans_D_1 = [zeros(1, 100), signal_trans_1]; %为了仿真
signal_trans_D_2 = [zeros(1, t_delay(2)), signal_trans_2];
signal_trans_D_3 = [zeros(1, t_delay(3)), signal_trans_3];
signal_trans_D_4 = [zeros(1, t_delay(4)), signal_trans_4];

len_total = max([length(signal_trans_D_1), length(signal_trans_D_2), length(signal_trans_D_3), length(signal_trans_D_4)]);

signal_trans_D_1 = [signal_trans_D_1, zeros(1, len_total - length(signal_trans_D_1))];
signal_trans_D_2 = [signal_trans_D_2, zeros(1, len_total - length(signal_trans_D_2))];
signal_trans_D_3 = [signal_trans_D_3, zeros(1, len_total - length(signal_trans_D_3))];
signal_trans_D_4 = [signal_trans_D_4, zeros(1, len_total - length(signal_trans_D_4))];

signal_trans_mixed = zeros(1, len_total);
signal_trans_mixed = signal_trans_mixed + signal_trans_D_1 + signal_trans_D_2 + signal_trans_D_3 + signal_trans_D_4;

% % 时域
% figure;
% plot(real(signal_trans_mixed));
% % 频域
% figure;
% plot(20*log10(abs(fft(signal_trans_mixed))));

% 加噪�??
Es_N0 = 0;
SNRdB = Es_N0 - 10*log10(oversamp_IF);
signal_trans_mixed_noise = awgn(signal_trans_mixed, SNRdB, 'measured');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               接收部分
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

result = receiver(signal_trans_D_1, fh_pat_lib, th_pat_lib, pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 4);
