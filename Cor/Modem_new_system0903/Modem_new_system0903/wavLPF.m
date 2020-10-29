function [wav_S1_1_f, wav_S2_1_f, wav_S1_2_f, wav_S2_2_f] = wavLPF(wav_S1_1, wav_S2_1, wav_S1_2, wav_S2_2, LPF, S_lpf)

% 将本地生成的对应各同步序列的波形经过接收部分的基带信号低通滤波器，将最终经过低通滤波后且一个符号对应1个最佳采样点的波形，用于与同步头接收波形做互相关运算

num_pulses = size(wav_S1_1, 1);
oversamp = 8;

wav_S1_1_f = zeros(size(wav_S1_1,1), size(wav_S1_1, 2)/8);
wav_S2_1_f = zeros(size(wav_S2_1,1), size(wav_S2_1, 2)/8);
wav_S1_2_f = zeros(size(wav_S1_2,1), size(wav_S1_2, 2)/8);
wav_S2_2_f = zeros(size(wav_S2_2,1), size(wav_S2_2, 2)/8);

for i = 1:num_pulses
    S1_wav_1_temp = conv(wav_S1_1(i,:), LPF);
    S1_wav_1_f_temp = S1_wav_1_temp(S_lpf:S_lpf+length(wav_S1_1(i,:))-1);
    wav_S1_1_f(i,:) = S1_wav_1_f_temp(1:oversamp:end);

    S2_wav_1_temp = conv(wav_S2_1(i,:), LPF);
    S2_wav_1_f_temp = S2_wav_1_temp(S_lpf:S_lpf+length(wav_S2_1(i,:))-1);
    wav_S2_1_f(i,:) = S2_wav_1_f_temp(1:oversamp:end);
    
    S1_wav_2_temp = conv(wav_S1_2(i,:), LPF);
    S1_wav_2_f_temp = S1_wav_2_temp(S_lpf:S_lpf+length(wav_S1_2(i,:))-1);
    wav_S1_2_f(i,:) = S1_wav_2_f_temp(1:oversamp:end);

    S2_wav_2_temp = conv(wav_S2_2(i,:), LPF);
    S2_wav_2_f_temp = S2_wav_2_temp(S_lpf:S_lpf+length(wav_S2_2(i,:))-1);
    wav_S2_2_f(i,:) = S2_wav_2_f_temp(1:oversamp:end);
end
