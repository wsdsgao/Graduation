function result = receiver(rx, fh_pat_lib_1, th_pat_lib_1, ...
    wav_S1_S3_1_mode4,wav_S2_S4_1_mode4, ...
    frame_num,signal_BB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          参数定义
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode_sel = 2;  % 用于记录最终判断当前接收波形属于哪种速率模式

num_bits_pn = 24;  % 同步头S1\S2长度
num_bits_pn_2 = 21;  % 同步头S3\S4长度

bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间
fs_IF = 1024e6;  % 射频、中频信号采样速率
fs_BB = 128e6;  % 基带信号采样速率
oversamp_BB = T * fs_BB;  % 基带信号过采样倍数
oversamp_IF = T * fs_IF;  % 中频信号过采样倍数

num_bits_pulse = 304;  % 2Mbps A\500Kbps\250Kbps 一个脉冲包的长度
                       % 2Mbps B 每个脉冲去掉前后跳数据部分后的长度
frame_counter = 1;  % 记录已检测到多少帧数据（仿真用）
flag_frame = 0;  % 解调完一帧，该标志位置1，等待一帧时间后再捕获下一帧（仿真用）
time_counter_mode1 = 6720*oversamp_IF - 0.5*oversamp_IF;  % 2Mbps A模式 等待时长
time_counter_mode2 = 7956*oversamp_IF - 0.5*oversamp_IF;  % 2Mbps B模式 等待时长
time_counter_mode3 = 26880*oversamp_IF + 10*oversamp_IF - 0.5*oversamp_IF;  % 500Kbps模式 等待时长
time_counter_mode4 = 53760*oversamp_IF - 0.5*oversamp_IF;  % 250Kbps模式 等待时长
counter = 0;

flag_Capture_C = 0;  % 指示是否捕获到一帧

% 各模式一帧总长度
time_frame_mode1 = (304*12+512*6+103*12);
time_frame_mode2 = (304*12+512*6+103*12);
time_frame_mode3 = (304*12+512*6)*4;
time_frame_mode4 = (304*12+512*6)*8;

% 各模式一帧有效数据长度
length_frame_mode1 = 3072;
length_frame_mode2 = 6144;
length_frame_mode3 = 10272;
length_frame_mode4 = 20544;

% 各模式一帧脉冲数量
num_pulses_mode1 = 12;
num_pulses_mode2 = 12;
num_pulses_mode3 = 48;
num_pulses_mode4 = 96;

% 各模式跳频图案
fh_pat_1_mode1 = fh_pat_lib_1(1:num_pulses_mode1);
% fh_pat_2_mode1 = fh_pat_lib_2(1:num_pulses_mode1);
fh_pat_1_mode2 = fh_pat_lib_1(1:num_pulses_mode2);
% fh_pat_2_mode2 = fh_pat_lib_2(1:num_pulses_mode2);
fh_pat_1_mode3 = fh_pat_lib_1(1:num_pulses_mode3);
% fh_pat_2_mode3 = fh_pat_lib_2(1:num_pulses_mode3);
fh_pat_1_mode4 = fh_pat_lib_1(1:num_pulses_mode4);
% fh_pat_2_mode4 = fh_pat_lib_2(1:num_pulses_mode4);

% 各模式跳时图案
th_pat_1_mode1 = th_pat_lib_1(1:num_pulses_mode1);
% th_pat_2_mode1 = th_pat_lib_2(1:num_pulses_mode1);
th_pat_1_mode2 = th_pat_lib_1(1:num_pulses_mode2);
% th_pat_2_mode2 = th_pat_lib_2(1:num_pulses_mode2);
th_pat_1_mode3 = th_pat_lib_1(1:num_pulses_mode3);
% th_pat_2_mode3 = th_pat_lib_2(1:num_pulses_mode3);
th_pat_1_mode4 = th_pat_lib_1(1:num_pulses_mode4);
% th_pat_2_mode4 = th_pat_lib_2(1:num_pulses_mode4);


% 载入已存数据
% load('lib/f_trans.mat');  % 21个频点
freq_21=[240-40/3*10:40/3:240+40/3*10];%定义调试5个频点
freq_5=[240-80/3,240-40/3,240,240+40/3,240+80/3];%定义调试5个频点
freq_1=[240,240,240,240,240];%定义调试1个频点
f_trans=freq_21*1000000;

% 下变频

rx_pulse_128_LPF = DDC(rx);


%跳频图案映射
for i=1:96
    rx_pulse(i,:)=rx_pulse_128_LPF(fh_pat_lib_1(i),:);
end 

rx_pulse_mat=repmat(rx_pulse,1,8);

%跳时图案映射

%     for i = 1:48
%         signal_BB_d8(i,1:4500*8)=[zeros(1,500*8),signal_BB(i,th_pat_lib_1(i)*4*8+1:8:th_pat_lib_1(i)*4*8+1+8*4000*8-1)];
%     end 

    for i = 1:12
        signal_BB_d8(i,:)=[zeros(1,100*8),signal_BB(i,1:8:end),zeros(1,2000*8)];
    end     

% for i=1:24
%     figure;
%     plot(real(signal_BB(i,1:64:45*64)));
%     hold on;
%     plot(real(signal_BB(i,259*64+1:64:304*64)));
% %     pause;
%     close;
% end 
    
% for m=2035%1:8000
for m=1:8000
    pos_stat=m;
    pos_last=m;
    
%     for i=1:48
%         rx_pulse_fh(i,:)=rx_pulse_mat(i,pos_stat:pos_stat+304*8);
%         pos_last=pos_last+(th_pat_lib_1(i)/2+th_pat_lib_1(i+1)/2+103+304)*8;
%         pos_stat=pos_last+1;
%     end 

%     for i = 1:12
%         signal_BB_d8(i,1:4500*8)=[zeros(1,500*8),signal_BB(i,th_pat_lib_1(i)*4*8+1:8:th_pat_lib_1(i)*4*8+1+8*4000*8-1)];
%     end 
    rx_pulse_fh=signal_BB_d8(:,m:end);

    %2MA 求相关峰
    for i=1:12
%     for i=1       
        xorr_rx_mode1_1(1,:) = rx_pulse_fh(i,1:8:8*8) .*(wav_S1_S3_1_mode4(i,1:8));
        xorr_rx_mode1_2(1,:) = rx_pulse_fh(i,8*8+1:8:8*16) .* conj(wav_S1_S3_1_mode4(i,9:16));
        xorr_rx_mode1_3(1,:) = rx_pulse_fh(i,8*16+1:8:8*24) .* conj(wav_S1_S3_1_mode4(i,17:24));
        xorr_rx_mode1_sum(1,i)=sum(xorr_rx_mode1_1(1,:))+sum(xorr_rx_mode1_2(1,:))+sum(xorr_rx_mode1_3(1,:));
    end 

    xorr_rx_abs_mode1(1,m)=sum(abs(xorr_rx_mode1_sum(1,:)));
    
    %2MB 求相关峰
    for i=1:12
%     for i=1        
        xorr_rx_mode2_1(1,:) = rx_pulse_fh(i,1:8:8*8) .*conj(wav_S1_S3_1_mode4(i,1:8));
        xorr_rx_mode2_2(1,:) = rx_pulse_fh(i,8*8+1:8:8*16) .*(wav_S1_S3_1_mode4(i,9:16));
        xorr_rx_mode2_3(1,:) = rx_pulse_fh(i,8*16+1:8:8*24) .*conj(wav_S1_S3_1_mode4(i,17:24));
        xorr_rx_mode2_sum(1,i)=sum(xorr_rx_mode2_1(1,:))+sum(xorr_rx_mode2_2(1,:))+sum(xorr_rx_mode2_3(1,:));
    end 

    xorr_rx_abs_mode2(1,m)=sum(abs(xorr_rx_mode2_sum(1,:)));    

%     %500K 求相关峰
%     for i=1:24
%         xorr_rx_mode3_1(1,:) = rx_pulse_fh(i,1:8:24*8) .*conj(wav_S1_S3_1_mode4(i,1:24));
%         xorr_rx_mode3_2(1,:) = rx_pulse_fh(i,24*8+1:8:24*8+21*8) .*(wav_S1_S3_1_mode4(i,25:45));
%         xorr_rx_mode3_sum(1,i)=sum(xorr_rx_mode3_1(1,:))+sum(xorr_rx_mode3_2(1,:));
%     end 
% 
%     xorr_rx_abs_mode3(1,m)=sum(abs(xorr_rx_mode3_sum(1,:)));    
% 
%     %250K 求相关峰
%     for i=1:24
%         xorr_rx_mode4_1(1,:) = rx_pulse_fh(i,1:8:24*8) .*conj(wav_S1_S3_1_mode4(i,1:24));
%         xorr_rx_mode4_2(1,:) = rx_pulse_fh(i,24*8+1:8:24*8+21*8) .*conj(wav_S1_S3_1_mode4(i,25:45));
%         xorr_rx_mode4_sum(1,i)=sum(xorr_rx_mode4_1(1,:))+sum(xorr_rx_mode4_2(1,:));
%     end 
% 
%     xorr_rx_abs_mode4(1,m)=sum(abs(xorr_rx_mode4_sum(1,:)));    
%     
end 

% for i = 1:12
%     figure;
%     plot(real(signal_BB(i,1:8:end))*1000);
%     hold on;
%     plot(real(rx_pulse_fh(i,1:end)));
% end 

figure;
subplot(4,1,1)
plot(xorr_rx_abs_mode1);
subplot(4,1,2)
plot(xorr_rx_abs_mode2);
% subplot(4,1,3)
% plot(xorr_rx_abs_mode3);
% subplot(4,1,4)
% plot(xorr_rx_abs_mode4);

% 
xorr_rx_abs_mode1_12=xorr_rx_abs_mode1(1,1:end-2240)+xorr_rx_abs_mode1(1,2241:end);
xorr_rx_abs_mode2_12=xorr_rx_abs_mode2(1,1:end-2240)+xorr_rx_abs_mode2(1,2241:end);
% xorr_rx_abs_mode3_12=xorr_rx_abs_mode3(1,1:end-2072)+xorr_rx_abs_mode3(1,2073:end);
% xorr_rx_abs_mode4_12=xorr_rx_abs_mode4(1,1:end-2072)+xorr_rx_abs_mode4(1,2073:end);
%     
figure;
subplot(4,1,1)
plot(xorr_rx_abs_mode1_12);
subplot(4,1,2)
plot(xorr_rx_abs_mode2_12);
% subplot(4,1,3)
% plot(xorr_rx_abs_mode3_12);
% subplot(4,1,4)
% plot(xorr_rx_abs_mode4_12);

% 
% figure;
% plot(xorr_rx_abs_mode1);
% hold on;
% plot(xorr_rx_abs_mode2);
% figure;
% plot(xorr_rx_abs_mode3_12);
% figure;
% plot(xorr_rx_abs_mode4_12);
% 
% xorr_rx_abs_mode34_12=xorr_rx_abs_mode3_12;

result=xorr_rx_abs_mode1_12;







