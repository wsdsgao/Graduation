clear all; clc;

% GMSK 解调性能仿真 减状态数Viterbi解调性能测试 跳时+频偏(0.001*Rb) 
%       速率模式        4种：2M@A, 2M@B, 500K, 250K
%       帧结构          2M@A:      24bit（同步头S1） + 256bit（数据） + 24bit（同步头S2） 
%                      2M@B:      0~Ln/2bit(数据) + 24bit（同步头S1） + 256bit（数据） + 24bit（同步头S2） + 0~Ln/2bit(数据)
%                      500K\250K:  24bit（同步头S1）+ 21bit（同步头S3）+ 214bit（数据）+ 21bit（同步头S4） + 24bit（同步头S2）
%       脉冲数量        2M@A\2M@B:   12           
%                      500K:  48
%                      250K:  96
%       一帧数据位数     2M@A:   256*12 = 3072 chips
%                      2M@B:   256*12 + 512*6 = 6144 chips
%                      500K:    214*48 = 10272 chips（包含32bit 0）
%                      250K:    214*96 = 20544 chips（包含64bit 0）

% 参数说明： 
%       mode:  ==1: 2M@A;  
%              ==2: 2M@B;
%              ==3: 500K;
%              ==4: 250K
%

% 仿真数据量:  2帧数据 
%            2Mbps@A 模式: 5*3072 = 15360 bit
%            2Mbps@B 模式: 5*6144 = 30720 bit
%            500Kbps: 5*10272 = 51360 bit
%            250Kbps: 5*20544 = 102720 bit
%

mode = 1;  

frame_num = 5;  % 仿真多少帧数据
num_bits_pn = 24;  % 同步头S1\S2长度
num_bits_pn_2 = 21;  % 同步头S3\S4长度

bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间
num_bits_pulse = 304;  % 2M_A\500K\250K脉冲长度（包含前后同步头）, 2M_B 脉冲去掉前后跳数据部分后的长度

fs_IF = 1024e6;  % 射频、中频信号采样速率
fs_BB = 128e6;  % 基带信号采样速率
oversamp_BB = T * fs_BB;  % 基带信号过采样速率
oversamp_IF = T * fs_IF;  % 射频、中频信号过采样速率
T_s_BB = 1/fs_BB;  % 基带采样间隔
T_s_IF = 1/fs_IF;  % 射频、中频采样间隔
BER = zeros(1,15); %误比特率
Eb_N0 = 1: 1: 14;

% 载入已存数据
load('lib/g_1024.mat');  % GMSK调制 g函数(高斯低通滤波器)
load('lib/f_trans.mat');  % 21个频点（由PN序列得）
% 接收端滤波器
load('lib/filter/BPF_CHAN1_UD.mat');  % 带通 通道1
load('lib/filter/BPF_CHAN2_UD.mat');  % 带通 通道2
load('lib/filter/BPF_CHAN3_UD.mat');  % 带通 通道3
load('lib/filter/BPF_CHAN4_UD.mat');  % 带通 通道4
load('lib/filter/BPF_CHAN5_UD.mat');  % 带通 通道5
load('lib/filter/BPF_CHAN12_2_UD.mat');  % 带通 通道1、2 第二级滤波
load('lib/filter/BPF_CHAN34_2_UD.mat');  % 带通 通道3、4 第二级滤波
load('lib/filter/BPF_CHAN5_2_UD.mat');  % 带通 通道5 第二级滤波
load('lib/filter/c0_f_128.mat');  % 解调匹配滤波器1
load('lib/filter/c1_f_128.mat');  % 解调匹配滤波器2
load('lib/filter/LPF.mat');  % 低通滤波器 61阶 通带: 4.8MHz 采样频率128MHz
load('lib/filter/LPF_2.mat');  % 低通滤波器 507阶 通带: 5MHz  采样频率1024MHz
S_lpf = 30;
S_lpf2 = 127;
S_bpf = 253;

% 每帧分成的脉冲包个数 根据速率模式选择
if mode == 1
    num_pulses = 12;
elseif mode == 2
    num_pulses = 12;
elseif mode == 3
    num_pulses = 48;
else 
    num_pulses = 96;
end

% 仿真的总脉冲个数 (一共frame_num帧)
if mode == 2
    mat_row = frame_num;
else
    mat_row = num_pulses * frame_num;  
end

% 生成PN库、跳频图案、跳时图案、生成数据矩阵 
bits = data_gen(mat_row, mode);   % 双极性码（还包含了同步头，只有B没有）
[th_pat_lib_1, fh_pat_lib_1] = TF_gen;   % 跳频、跳时总图案1
[th_pat_lib_2, fh_pat_lib_2] = TF_gen;   % 跳频、跳时总图案2

[pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1] = pn_gen;  % 0\1 码
[pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2] = pn_gen;

% 同步头波形
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

% 同步头低通滤波后波形
% mode 1
[wav_S1_1_LPF_mode1, wav_S2_1_LPF_mode1, wav_S1_2_LPF_mode1, wav_S2_2_LPF_mode1] = wavLPF(wav_S1_1_mode1, wav_S2_1_mode1, wav_S1_2_mode1, wav_S2_2_mode1, LPF, S_lpf);
% mode 2
[wav_S1_1_LPF_mode2, wav_S2_1_LPF_mode2, wav_S1_2_LPF_mode2, wav_S2_2_LPF_mode2] = wavLPF(wav_S1_1_mode2, wav_S2_1_mode2, wav_S1_2_mode2, wav_S2_2_mode2, LPF, S_lpf);
% mode 3
[wav_S1_S3_1_LPF_mode3, wav_S4_S2_1_LPF_mode3, wav_S1_S3_2_LPF_mode3, wav_S4_S2_2_LPF_mode3] = wavLPF(wav_S1_S3_1_mode3, wav_S4_S2_1_mode3, wav_S1_S3_2_mode3, wav_S4_S2_2_mode3, LPF, S_lpf);
% mode 4
[wav_S1_S3_1_LPF_mode4, wav_S4_S2_1_LPF_mode4, wav_S1_S3_2_LPF_mode4, wav_S4_S2_2_LPF_mode4] = wavLPF(wav_S1_S3_1_mode4, wav_S4_S2_1_mode4, wav_S1_S3_2_mode4, wav_S4_S2_2_mode4, LPF, S_lpf);

% 检查矩阵 （用于解调后统计误码率）
if mode == 1
    for i = 1:mat_row/num_pulses  %帧数
        for j = 1:num_pulses
            bit_frame((j-1)*256+1:j*256) = bits((i-1)*num_pulses+j,25:280);  %相当于挑出数据位，重新排列成一行
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


% 信噪比0:17 循环
for SNR_idx = 10:10 %1:length(Eb_N0)
   
    SNR_idx

    N0 = 1/10^(Eb_N0(SNR_idx)/10);  % 计算当前信噪比对应的高斯白噪声功率谱密度值 N0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                 发射部分
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % 生成GMSK波形
    signal_trans = transmitter(bits, fh_pat_lib_1, th_pat_lib_1, fh_pat_lib_2, th_pat_lib_2, pn_lib_S1_1, pn_lib_S1_2,...
    pn_lib_S2_1, pn_lib_S2_2, pn_lib_S3_1, pn_lib_S3_2, pn_lib_S4_1, pn_lib_S4_2, mode);

    % 前后补零
    signal_trans_2 = [zeros(1,0.5*oversamp_IF), signal_trans, zeros(1,100000*oversamp_IF)];  % 加0.5个符号的空白
    
    % 加高斯白噪声
    White_N = sqrt(N0*oversamp_IF/2)*(randn(1, length(signal_trans_2))+1i*randn(1, length(signal_trans_2)));
    sig_N = signal_trans_2 + White_N; 
    
    % 加频偏 (0.016MHz)(千分之一倍数据速率)
    t = 0:1/oversamp_IF:length(signal_trans_2)/oversamp_IF-1/oversamp_IF;
    t = t(1:end) + (T/oversamp_IF/2);  % 中心对称
    pha = 0.016e6*2*pi*t*T; 
    pha = mod(pha, 2*pi);
    pha_cplx = complex(cos(pha), -sin(pha));
    rx = sig_N .* pha_cplx;  % rx: 为送入解调模块的带有频偏的中频信号

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               接收部分
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




