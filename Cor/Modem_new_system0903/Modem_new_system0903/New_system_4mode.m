clear all; 
close all;

% GMSK 解调性能仿真 减状态数Viterbi解调性能测试 跳时+频偏(0.001*Rb) 
%       速率模式        4种：2M@A, 2M@B, 500K, 250K
%       帧结构          2M@A:      24bit（同步头S1） + 256bit（数据） + 24bit（同步头S2） 
%                      2M@B:      0~Ln/2bit(数据) + 24bit（同步头S1） + 256bit（数据） + 24bit（同步头S2） + 0~Ln/2bit(数据)
%                      500K\250K:  24bit（同步头S1）+ 21bit（同步头S3）+ 214bit（数据）+ 21bit（同步头S4） + 24bit（同步头S2）
%       脉冲数量        2M@A\2M@B:   12           
%                      500K:  48
%                      250K:  96
%       一帧数据位数    2M@A:   256*12 = 3072 chips
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

frame_num = 1;  % 仿真多少帧数据
num_bits_pn = 24;  % 同步头S1\S2长度
num_bits_pn_2 = 21;  % 同步头S3\S4长度

bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间
num_bits_pulse = 304;  % 2M_A\500K\250K脉冲长度（包含前后同步头）, 2M_B 脉冲去掉前后跳数据部分后的长度

fs_IF = 1024e6;  % 射频、中频频信号采样速率
fs_BB = 128e6;  % 基带信号采样速率
oversamp_BB = T * fs_BB;  % 基带信号过采样速率
oversamp_IF = T * fs_IF;  % 射频、中频信号过采样速率
T_s_BB = 1/fs_BB;  % 基带采样间隔
T_s_IF = 1/fs_IF;  % 射频、中频采样间隔
BER = zeros(1,15);
Eb_N0 = 0: 1: 14;

% 载入已存数据
load('lib/g_1024.mat');  % GMSK调制 g函数 
load('lib/f_trans.mat');  % 21个频点
load('lpf_coe.mat');  % 21个频点
f_trans=80000000*ones(1,21);

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

% LPF=[0.00139356394566424,0.00115641684212128,0.000757467920388439,-0.000469440499186030,-0.00232677464278785,-0.00414949860891374,-0.00494457523887240,-0.00377093447043267,-0.000247545649598074,0.00501881429030421,0.0103288738191063,0.0132636528084217,0.0115211209418782,0.00398526254669740,-0.00836588509437527,-0.0221342250426825,-0.0320312770507913,-0.0322494464576837,-0.0183331547312239,0.0110468950202278,0.0530627576908241,0.100906311211473,0.145239660105501,0.176580265525652,0.187891535968573,0.176580265525652,0.145239660105501,0.100906311211473,0.0530627576908241,0.0110468950202278,-0.0183331547312239,-0.0322494464576837,-0.0320312770507913,-0.0221342250426825,-0.00836588509437527,0.00398526254669740,0.0115211209418782,0.0132636528084217,0.0103288738191063,0.00501881429030421,-0.000247545649598074,-0.00377093447043267,-0.00494457523887240,-0.00414949860891374,-0.00232677464278785,-0.000469440499186030,0.000757467920388439,0.00115641684212128,0.00139356394566424];
% S_lpf=24;


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
bits = data_gen(mat_row, mode);   % 双极性码
% [th_pat_lib_1, fh_pat_lib_1] = TF_gen;   % 跳频、跳时总图案1
% [th_pat_lib_2, fh_pat_lib_2] = TF_gen;   % 跳频、跳时总图案2
% 
% [pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1] = pn_gen;  % 0\1 码
% [pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2] = pn_gen;

% 固定跳频和调试图案，调试用
% th_pat_lib_1=[294,354,248,138,200,258,218,158,264,374,312,254,84,460,162,124,162,318,428,52,350,388,350,194,110,420,90,244,92,140,402,92,422,268,420,372,256,240,34,124,120,164,256,272,478,388,392,348,294,354,248,138,200,258,218,158,264,374,312,254,84,460,162,124,162,318,428,52,350,388,350,194,110,420,90,244,92,140,402,92,422,268,420,372,256,240,34,124,120,164,256,272,478,388,392,348];
th_pat_lib_1 = [220,216,32,372,208,164,292,296,480,140,304,348,262,408,426,446,300,266,250,104,86,66,212,246,172,126,126,296,214,44,340,386,386,216,298,468,422,32,264,328,64,356,90,480,248,184,448,156,366,86,278,130,214,346,146,426,234,382,298,166,354,90,92,222,108,72,158,422,420,290,404,440,344,344,432,352,326,378,168,168,80,160,186,134,166,174,476,106,72,422,346,338,36,406,440,90];
th_pat_lib_1_12=th_pat_lib_1/2;

th_pat_lib_2=[294,354,248,138,200,258,218,158,264,374,312,254,84,460,162,124,162,318,428,52,350,388,350,194,110,420,90,244,92,140,402,92,422,268,420,372,256,240,34,124,120,164,256,272,478,388,392,348,294,354,248,138,200,258,218,158,264,374,312,254,84,460,162,124,162,318,428,52,350,388,350,194,110,420,90,244,92,140,402,92,422,268,420,372,256,240,34,124,120,164,256,272,478,388,392,348];
th_pat_lib_2_12=th_pat_lib_2/2;

% fh_pat_lib_1=[1,9,21,20,19,10,16,20,8,3,12,8,18,4,5,12,2,6,14,3,13,19,11,7,17,21,14,11,10,4,10,11,19,16,14,2,3,13,10,12,1,13,14,11,4,16,17,7,11,5,5,17,2,5,21,13,13,20,3,19,19,5,4,14,4,3,8,9,2,15,12,21,5,20,14,16,8,10,3,19,10,15,16,6,13,8,17,2,21,17,1,19,19,15,19,10];
% fh_pat_lib_1=[0,4,5,3,1,7,3,5,1,6,2,6,6,7,4,0,2,9,1,8,5,9,0,4,1,9,0,7,8,8,0,3,2,8,4,9,1,2,1,1,8,5,5,1,8,6,3,5,4,0,2,1,1,2,4,0,9,9,4,4,3,9,3,1,7,3,2,4,0,1,9,9,5,0,2,3,8,0,0,1,6,7,6,4,5,2,7,1,6,1,3,6,7,0,9,7];
% fh_pat_lib_1 = [9,1,11,4,10,2,3,6,14,12,8,5,13,15,7,9,1,11,4,10,2,3,6,14,12,8,5,13,15,7,9,1,11,4,10,2,3,6,14,12,8,5,13,15,7,9,1,11,4,10,2,3,6,14,12,8,5,13,15,7,9,1,11,4,10,2,3,6,14,12,8,5,13,15,7,9,1,11,4,10,2,3,6,14,12,8,5,13,15,7,9,1,11,4,10,2];
fh_pat_lib_1 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,1,2,3,4,5,6];
fh_pat_lib_1 = fh_pat_lib_1;
fh_pat_lib_1_fpga = fh_pat_lib_1-1;
% fh_pat_lib_1=randint([1,21],1,96);
% fh_pat_lib_2=[1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1];
fh_pat_lib_2=fh_pat_lib_1;

% fh_pat_lib_1=[10,1,19,12,18,21,13,7,11,16,8,15,14,2,6,5,20,3,17,9,4,19,4,9,7,10,15,3,18,21,2,11,16,1,5,12,6,20,17,8,14,13,21,10,6,17,14,5,10,1,19,12,18,21,13,7,11,16,8,15,14,2,6,5,20,3,17,9,4,19,4,9,7,10,15,3,18,21,2,11,16,1,5,12,6,20,17,8,14,13,21,10,6,17,14,5];
% fh_pat_lib_2=[10,1,19,12,18,21,13,7,11,16,8,15,14,2,6,5,20,3,17,9,4,19,4,9,7,10,15,3,18,21,2,11,16,1,5,12,6,20,17,8,14,13,21,10,6,17,14,5,10,1,19,12,18,21,13,7,11,16,8,15,14,2,6,5,20,3,17,9,4,19,4,9,7,10,15,3,18,21,2,11,16,1,5,12,6,20,17,8,14,13,21,10,6,17,14,5];

% [pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1] = pn_gen;  % 0\1 码
% [pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2] = pn_gen;
sync_gen_poly=[0,0,1,1,0,0,0,0,0,0,0,0,1]; %x13+x4+x3+x+1
sync_int_phase=[0,1,0,0,1,1,1,0,0,1,0,1,0];
sync_m_seq=m_sequence( sync_gen_poly,sync_int_phase);

% s1=[sync_m_seq(1:1152),sync_m_seq(1:1152)];     %生成S1序列 15*24 = 360
% % s2=[sync_m_seq(1:1152),sync_m_seq(1:1152)];
% s2=[sync_m_seq(1153:2304),sync_m_seq(1153:2304)];   
% s3=[sync_m_seq(2305:3312),sync_m_seq(2305:3312)];
% % s4=[sync_m_seq(2305:3312),sync_m_seq(2305:3312)];
% s4=[sync_m_seq(3313:4320),sync_m_seq(3313:4320)];

s1=[sync_m_seq(1:360)];     %生成S1序列 15*24 = 360
% s2=[sync_m_seq(1:1152),sync_m_seq(1:1152)];
s2=[sync_m_seq(361:720)];	 %生成S2序列 15*24 = 360
s3=[sync_m_seq(721:1035)];   %生成S3序列 15*21 = 315
% s4=[sync_m_seq(2305:3312),sync_m_seq(2305:3312)];
s4=[sync_m_seq(1036:1350)];  %生成S4序列 15*21 = 315

pn_lib_S1_1=reshape(s1,[24,15])';
pn_lib_S2_1=reshape(s1,[24,15])';
pn_lib_S3_1=reshape(s3,[21,15])';
pn_lib_S4_1=reshape(s3,[21,15])';

pn_lib_S1_2=reshape(s1,[24,15])';
pn_lib_S2_2=reshape(s1,[24,15])';
pn_lib_S3_2=reshape(s3,[21,15])';
pn_lib_S4_2=reshape(s3,[21,15])';

% % 同步头波形
% % mode 1
% [wav_S1_1_mode1, wav_S2_1_mode1] = wav_gen(pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1, 1);
% [wav_S1_2_mode1, wav_S2_2_mode1] = wav_gen(pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2, 1);
% % mode 2
% [wav_S1_1_mode2, wav_S2_1_mode2] = wav_gen(pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1, 2);
% [wav_S1_2_mode2, wav_S2_2_mode2] = wav_gen(pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2, 2);
% % mode 3
% [wav_S1_S3_1_mode3, wav_S4_S2_1_mode3] = wav_gen(pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1, 3);
% [wav_S1_S3_2_mode3, wav_S4_S2_2_mode3] = wav_gen(pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2, 3);

% % mode 4  生成同步头波形
[wav_S1_S3_1_mode4, wav_S2_S4_1_mode4] = wav_gen_new(pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1, 4);
[wav_S1_S3_2_mode4, wav_S4_S2_2_mode4] = wav_gen_new(pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2, 4); 

% figure;
% plot(real(wav_S1_S3_1_mode4(1,:)));
% hold on;
% plot(imag(wav_S1_S3_1_mode4(1,:)));

% figure;
% subplot(2,1,1)
% plot(real(wav_S1_S3_1_mode4(2,:)),'b');
% hold on;
% subplot(2,1,2)
% plot(real(wav_S2_S4_1_mode4(2,:)),'r');

% % 同步头低通滤波后波形
% % mode 1
% [wav_S1_1_LPF_mode1, wav_S2_1_LPF_mode1, wav_S1_2_LPF_mode1, wav_S2_2_LPF_mode1] = wavLPF(wav_S1_1_mode1, wav_S2_1_mode1, wav_S1_2_mode1, wav_S2_2_mode1, LPF, S_lpf);
% % mode 2
% [wav_S1_1_LPF_mode2, wav_S2_1_LPF_mode2, wav_S1_2_LPF_mode2, wav_S2_2_LPF_mode2] = wavLPF(wav_S1_1_mode2, wav_S2_1_mode2, wav_S1_2_mode2, wav_S2_2_mode2, LPF, S_lpf);
% % mode 3
% [wav_S1_S3_1_LPF_mode3, wav_S4_S2_1_LPF_mode3, wav_S1_S3_2_LPF_mode3, wav_S4_S2_2_LPF_mode3] = wavLPF(wav_S1_S3_1_mode3, wav_S4_S2_1_mode3, wav_S1_S3_2_mode3, wav_S4_S2_2_mode3, LPF, S_lpf);
% % mode 4
% [wav_S1_S3_1_LPF_mode4, wav_S4_S2_1_LPF_mode4, wav_S1_S3_2_LPF_mode4, wav_S4_S2_2_LPF_mode4] = wavLPF(wav_S1_S3_1_mode4, wav_S4_S2_1_mode4, wav_S1_S3_2_mode4, wav_S4_S2_2_mode4, LPF, S_lpf);
% wav_S1_S3_1_LPF_mode4 = wavLPF_new(wav_S1_S3_1_mode4, LPF, S_lpf);

% 检查矩阵 （用于解调后统计误码率）
% if mode == 1
%     for i = 1:mat_row/num_pulses
%         for j = 1:num_pulses
%             bit_frame((j-1)*256+1:j*256) = bits((i-1)*num_pulses+j,25:280);
%         end
%         bit_check(i,:) = bit_frame;
%     end
% elseif mode == 2
%     bit_check = bits;   
% elseif mode == 3
%     for i = 1:mat_row/num_pulses
%         for j = 1:num_pulses
%             bit_frame((j-1)*214+1:j*214) = bits((i-1)*num_pulses+j,46:259);
%         end
%         bit_check(i,:) = bit_frame;
%     end
% else
%     for i = 1:mat_row/num_pulses
%         for j = 1:num_pulses
%             bit_frame((j-1)*214+1:j*214) = bits((i-1)*num_pulses+j,46:259);
%         end
%         bit_check(i,:) = bit_frame;
%     end
% end

%%%%调试用，注释
% 信噪比0:17 循环
% for SNR_idx = 1:length(Eb_N0)
%    
%     SNR_idx
% 
%     N0 = 1/10^(Eb_N0(SNR_idx)/10);  % 计算当前信噪比对应的高斯白噪声功率谱密度值 N0
%%%%调试用，注释

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                                 发射部分
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% 生成GMSK波形
[signal_trans,signal_BB_out]  = transmitter(bits, fh_pat_lib_1, th_pat_lib_1, fh_pat_lib_2, th_pat_lib_2, pn_lib_S1_1, pn_lib_S1_2,...
pn_lib_S2_1, pn_lib_S2_2, pn_lib_S3_1, pn_lib_S3_2, pn_lib_S4_1, pn_lib_S4_2, mode);

signal_BB_out=round(signal_BB_out*256);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               生成仿真文件

    figure;
    plot(20*log10(abs(fft(signal_trans(1,1:1024*64)))));

    %%%时频图
    figure;
    X=real(signal_trans);
    T=0:1/1024:length(X)/1024;
    [S,F,T,P] = spectrogram(X,256,32,256,1024);
    mesh(T,F,10*log10(P)); axis tight;
    colorbar;  
    view(0,90);
    xlabel('Time (us)'); ylabel('MHz');

    
%生成本地同步字，量化
signal_tx = signal_trans;
    
wav_mode4_int_1 = round(wav_S1_S3_1_mode4*127);     %量化8bit

wav_mode4_int_real_1 = reshape(real(wav_mode4_int_1)',1,45*15);    
wav_mode4_int_imag_1 = reshape(imag(wav_mode4_int_1)',1,45*15);        

wav_mode4_int_6_1 = round(wav_S1_S3_1_mode4*7);     %量化4bit

wav_mode4_int_real_6_1 = reshape(real(wav_mode4_int_6_1)',1,45*15);    
wav_mode4_int_imag_6_1 = reshape(imag(wav_mode4_int_6_1)',1,45*15); 

wav_mode4_int_2 = round(wav_S1_S3_2_mode4*127);     %量化8bit

wav_mode4_int_real_2 = reshape(real(wav_mode4_int_2)',1,45*15);    
wav_mode4_int_imag_2 = reshape(imag(wav_mode4_int_2)',1,45*15);        

wav_mode4_int_6_2 = round(wav_S1_S3_2_mode4*7);     %量化4bit

wav_mode4_int_real_6_2 = reshape(real(wav_mode4_int_6_2)',1,45*15);    
wav_mode4_int_imag_6_2 = reshape(imag(wav_mode4_int_6_2)',1,45*15); 

figure;
plot(wav_mode4_int_real_1);
hold on;
plot(wav_mode4_int_real_6_1*20);

figure;
plot(wav_mode4_int_real_2);
hold on;
plot(wav_mode4_int_real_6_2*20);



%     IF_signal_tx = resample(signal_tx,128,1024);    %重采样到320MHz

%     
%     df=0;   %设置频率偏差
%     dp=0; %设置相位偏差
%     dt=[0:length(IF_signal_tx)-1];
%     IF_signal_tx_phase=IF_signal_tx.*exp(-j*(dp)*pi);       %加相位偏差
%     IF_signal_tx_freq=IF_signal_tx_phase.*exp(-j*(df/320000000)*2*pi*dt);   %加频偏
%     
%     IF_signal_tx_noise=awgn(real(IF_signal_tx_freq),100,'measured');        %加噪声
% 
%     IF_signal_tx_noise_1=round(IF_signal_tx/max(IF_signal_tx)*2048);    %量化
%     IF_signal_tx_noise_1=round(IF_signal_tx*256);    %量化
%     
%     IF_signal_tx_p2 = reshape(IF_signal_tx_noise_1',2,length(IF_signal_tx_noise_1)/2);
%     
%     IF_signal_tx_p2_AD = IF_signal_tx_p2';    
%     
%     m_file_name = 'NewFramedata_mode_A_320M_100dB_2048_F0_P0_dly10.txt';
%     
%     csvwrite(m_file_name,IF_signal_tx_p2_AD);

    IF_signal_tx = signal_tx(1:16:end);     %按4倍采样抽取
    IF_signal_tx_noise = awgn(IF_signal_tx,100,'measured'); %加噪声，100为信噪比
    IF_signal_tx_lpf = conv(IF_signal_tx_noise,lpf_coe);    %过低通滤波器
    IF_signal_tx_int = round((IF_signal_tx_lpf/max(IF_signal_tx_lpf))*256);  %量化
    
    figure;
    plot(real(IF_signal_tx_int));
    figure;
    plot(20*log10(abs(fft(IF_signal_tx_int))));
    
    IF_signal_tx_p2 = [real(IF_signal_tx_int)' imag(IF_signal_tx_int)'];    
    
    m_file_name = 'NewFramedata_mode_A_PNbond_64M_BB_100.txt';
%     csvwrite(m_file_name,IF_signal_tx_p2);
    dlmwrite(m_file_name,IF_signal_tx_p2,'delimiter','\t')

    
%     agc_base = IF_signal_tx;
%     agc_base_res=resample(agc_base,5,1);
%     ddt=[0:length(agc_base_res)-1];
%     agc_base_up=agc_base_res.*exp(-j*1/2*pi*ddt);
%     agc_base_ad=round(real(agc_base_up)*4000);
%     
%     noise_base=ones(1,length(agc_base_up))*2048;
%     noise_base_sig=awgn(noise_base,20,'measured');
%     noise=noise_base_sig-noise_base;
%     
%     agc_base_noise=round(agc_base_ad+noise);
%     agc_base_ad_p2=reshape(agc_base_noise,2,length(agc_base_noise)/2)';
%     figure;
%     plot(real(agc_base_noise));
% 
%     m_file_name = 'AD_agc_20dB_4000.txt';
%     csvwrite(m_file_name,agc_base_ad_p2);
%    
% 
%     dlmwrite(m_file_name,agc_base_ad_p2,'delimiter','\t')
  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                               接收部分
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rx=real(signal_trans);

rx=round(real(signal_trans)/max(real(signal_trans))*2048);    %量化
rx_noise=round(awgn(rx,-18,'measured'));        %加噪声

result = receiver_new(rx_noise, fh_pat_lib_1, th_pat_lib_1, ...
                    wav_S1_S3_1_mode4, wav_S2_S4_1_mode4,...
                    frame_num,signal_BB_out);

      
    [a,b] = size(bit_check); 
    error = sum(sum(abs(bit_check-result)/2 == 1));
    miss = sum(sum(abs(bit_check-result)/2 == 0.5));
    total = a*b-miss;
    pe(SNR_idx)=error/(total);
    
% end

a=load('resamp_out_100dB_2048_F5K_P0p5pi_dly10.txt');
resamp_out=complex(a(1:65000,1),a(1:65000,2));
resamp_pos=(1:length(resamp_out)-1);

resamp_out_phase=resamp_out*exp(-j*(-0.1)*pi);


% figure;
% plot(real(resamp_out(1:1520)));
% hold on;
% plot(imag(resamp_out(1:1520)));

abs_resamp=abs(resamp_out);
% 
% figure;
% plot(abs_resamp);


for n=1332
    st_pos=n;
    sync_out(1,:)=resamp_out(st_pos:8:st_pos+8*304-1);
%     resamp_dt(1,:)=resamp_pos(st_pos:8:st_pos+8*304-1);
    for i=2:12
        x(i,:)=(sum(th_pat_lib_1_12(1:i-1))+(304+103)*(i-1)+sum(th_pat_lib_1_12(2:i)))*8;
        sync_out(i,:)=resamp_out(st_pos+x(i,1):8:st_pos+8*304-1+x(i,1));
%         resamp_dt(1,:)=resamp_pos(st_pos+x(i,1):8:st_pos+8*304-1+x(i,1));
    end
 
    
    
%     figure;
%     plot(real(sync_out(1,1:24)));
%     hold on;
%     plot(imag(sync_out(1,1:24)));
% 
%     figure;
%     plot(real(sync_word_s1(1,:)));
%     hold on;
%     plot(imag(sync_word_s1(1,:)));
    
    rx_sync_word_s1=sync_out(:,3:24);
    rx_sync_word_s2=sync_out(:,283:304);

%     for i=1:12
%         rx_sw_xorr_s1(i,:)=rx_sync_word_s1(i,1:22).*conj(wav_S1_1_LPF_mode1(i,1:22));
%         rx_sw_xorr_s2(i,:)=rx_sync_word_s2(i,1:22).*conj(wav_S2_1_LPF_mode1(i,1:22));
%     end

%     rx_corr_S1_pat_1_mode1 = zeros(1,12);
%     rx_corr_S2_pat_1_mode1 = zeros(1,12);
%     for j = 1:12
%         rx_wav_S1_pat_1_mode1 = rx_sync_word_s1(j,1:22) .* conj(sync_word_s1(j,1:22));
%         rx_wav_S2_pat_1_mode1 = rx_sync_word_s2(j,1:22) .* conj(sync_word_s2(j,1:22));   
% 
%         rx_corr_S1_pat_1_mode1(j) = abs(sum(rx_wav_S1_pat_1_mode1));
%         rx_corr_S2_pat_1_mode1(j) = abs(sum(rx_wav_S2_pat_1_mode1));
%     end   
% 
%     rx_corr_S1_sum(n)=sum(rx_corr_S1_pat_1_mode1);
%     rx_corr_S2_sum(n)=sum(rx_corr_S2_pat_1_mode1);
%     
%     rx_corr_s1_S2_sum(n)=rx_corr_S1_sum(n)+rx_corr_S2_sum(n);
    
    for i=1:12
        rx_sw_xorr_s1(i,:)=rx_sync_word_s1(i,1:22).*conj(sync_word_s1(i,1:22));
        rx_sw_xorr_s2(i,:)=rx_sync_word_s2(i,1:22).*conj(sync_word_s2(i,1:22));
        rx_sw_xorr_delta(i,:)=rx_sw_xorr_s1(i,:).*conj(rx_sw_xorr_s2(i,:));
        real_sum(i)=real(sum(rx_sw_xorr_delta(i,:)));
        imag_sum(i)=-imag(sum(rx_sw_xorr_delta(i,:)));
        deltat_hat_S12(i) = (atan2(imag_sum(i),real_sum(i))); 
        deltat_hat_S12_1(i) = deltat_hat_S12(i)/pi*256
        delta_f_hat_S12(i) = deltat_hat_S12(i)/(2*pi*280);
%         rx_corr_S1_pat_1_mode1(i,n) = abs(sum(rx_sw_xorr_s1(i,:)));
%         rx_corr_S2_pat_1_mode1(i,n) = abs(sum(rx_sw_xorr_s2(i,:)));     
    end
    
    DF_hat_S12 = mean(delta_f_hat_S12); 

%     for i = 1:12
%         close all;
%         delta_theta_cplx_mat_S1 = rx_sync_word_s1(i,1:22) .* conj(sync_word_s1(i,:));
%         delta_theta_cplx_mat_S2 = rx_sync_word_s2(i,1:22) .* conj(sync_word_s2(i,:));  
%         delta_theta_cplx_S1 = sum(delta_theta_cplx_mat_S1(1:22)) / 22;
%         delta_theta_cplx_S2 = sum(delta_theta_cplx_mat_S2(1:22)) / 22;
%         delta_theta_S1(i) = atan2(-imag(delta_theta_cplx_S1), real(delta_theta_cplx_S1));
%         delta_theta_S2(i) = atan2(-imag(delta_theta_cplx_S2), real(delta_theta_cplx_S2));
% %         figure
% %         plot(real(delta_theta_cplx_mat_S1))
% %         hold on
% %         plot(imag(delta_theta_cplx_mat_S1))
% %         figure
% %         plot(real(delta_theta_cplx_mat_S2))
% %         hold on
% %         plot(imag(delta_theta_cplx_mat_S2))
%         phase_out(i,:) = sync_out(i,:).*exp(j*delta_theta_S1(i));
% %         figure
% %         plot(real(sync_out(i,1:24)))
% %         hold on
% %         plot(imag(sync_out(i,1:24)))
% %         figure
% %         plot(real(phase_out(i,1:24)))
% %         hold on
% %         plot(imag(phase_out(i,1:24)))
%     end    
    
    
    %
%     rx_sw_xorr_s1_sum(1,n)=sum(rx_corr_S1_pat_1_mode1(:,n));
%     rx_sw_xorr_s2_sum(1,n)=sum(rx_corr_S2_pat_1_mode1(:,n));
%     
%     rx_sw_xorr_sum(1,n)=rx_sw_xorr_s1_sum(1,n)+rx_sw_xorr_s2_sum(1,n);
end
