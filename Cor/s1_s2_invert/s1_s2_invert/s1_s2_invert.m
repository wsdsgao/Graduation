clear all;
close all;
load('lib/g_1024.mat');  % g函数
sync_gen_poly=[0,0,1,1,0,0,0,0,0,0,0,0,1]; %x13+x4+x3+x+1
sync_int_phase=[0,1,0,0,1,1,1,0,0,1,0,1,0];
sync_m_seq=m_sequence( sync_gen_poly,sync_int_phase);

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

[wav_S1_S3_1_mode4, wav_S2_S4_1_mode4] = wav_gen_new(pn_lib_S1_1, pn_lib_S2_1, pn_lib_S3_1, pn_lib_S4_1, 4);

pn_lib_S1_2=1-reshape(s1,[24,15])';
pn_lib_S2_2=1-reshape(s1,[24,15])';
pn_lib_S3_2=1-reshape(s3,[21,15])';
pn_lib_S4_2=1-reshape(s3,[21,15])';

[wav_S1_S3_2_mode4, wav_S2_S4_2_mode4] = wav_gen_new(pn_lib_S1_2, pn_lib_S2_2, pn_lib_S3_2, pn_lib_S4_2, 4);

% figure;
% plot(pn_lib_S1_1(1,:));
% hold on;
% plot(pn_lib_S1_2(1,:));
% 
% 
% figure;
% plot(real(wav_S1_S3_1_mode4(1,:)));
% hold on;
% plot(real(wav_S1_S3_2_mode4(1,:)));
% 
% figure;
% plot(imag(wav_S1_S3_1_mode4(1,:)));
% hold on;
% plot(imag(wav_S1_S3_2_mode4(1,:)));

wav_S1_S3_1=wav_S1_S3_1_mode4(1,:);
wav_S1_S3_2=wav_S1_S3_2_mode4(1,:);

noise_wav_S1_S3_1=awgn(wav_S1_S3_1,100,'measured'); %加噪声，100为信噪比

x=abs(wav_S1_S3_1*wav_S1_S3_1');

y=abs(wav_S1_S3_1*[wav_S1_S3_2(2:end),0]');

pn_lib_S1 = pn_lib_S1_1(1,:);
pn_lib_S2 = 1-pn_lib_S2_1(1,:);

preamble_S1 = 2*(pn_lib_S1)-1;
preamble_S2 = 2*(pn_lib_S2)-1;

f_IF=0;
phi_last=0;
oversamp_IF=4;

for jj = 1:length(preamble_S1)   % 每一个脉冲: (304+th+3)bit长度的波形

    % 对不同位置，取5bit数据，准备送入调制模块

    % 前数据段
    if (jj == 1)
        bit_5 = [preamble_S1(jj), preamble_S1(jj), preamble_S1(jj:jj+2)];
    elseif (jj == 2)
        bit_5 = [preamble_S1(jj-1), preamble_S1(jj-1:jj+2)];
    elseif (jj == 23)
        bit_5 = [preamble_S1(jj-2:jj+1),0];
    elseif (jj == 24)
        bit_5 = [preamble_S1(jj-2:jj), 0, 0];
    else
        bit_5 = preamble_S1(jj-2:jj+2);
    end

    [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, f_IF, phi_last, g);
    s1_BB(jj) = complex(I_sig(1), Q_sig(1));

end

for xx=1:12

    bits0=randi([0,1],1,128);
    bits1=randi([0,1],1,128);
    bits2=randi([0,1],1,128);
    preamble_bits0=2*bits0-1;
    preamble_bits1=2*bits1-1;
    preamble_bits2=2*bits2-1;

    bits_trans = [preamble_bits0, preamble_S1,preamble_bits1, preamble_S2,preamble_bits2]; 
    
    phi_last=0;
    phi_int=0;
    
    fp=rand(1,1)*2-1
    
    for jj = 1:length(bits_trans)   % 每一个脉冲: (304+th+3)bit长度的波形
        % 对不同位置，取5bit数据，准备送入调制模块

        % 前数据段
        if (jj == 1)
            bit_5 = [bits_trans(jj), bits_trans(jj), bits_trans(jj:jj+2)];
        elseif (jj == 2)
            bit_5 = [bits_trans(jj-1), bits_trans(jj-1:jj+2)];
        elseif (jj == 127)
            bit_5 = [bits_trans(jj-2:jj+1),0];
        elseif (jj == 128)
            bit_5 = [bits_trans(jj-2:jj), 0, 0];

        % S1
        elseif (jj == 129)
            phi_last = 0;
            bit_5 = [bits_trans(jj), bits_trans(jj), bits_trans(jj:jj+2)];
        elseif (jj == 130)
            bit_5 = [bits_trans(jj-1), bits_trans(jj-1:jj+2)];
        elseif (jj == 128+24-1)
            bit_5 = [bits_trans(jj-2:jj+1),0];
        elseif (jj == 128+24)
            bit_5 = [bits_trans(jj-2:jj), 0, 0];

        % 中间数据段    
        elseif (jj == 128+24+1)  % 首同步头第1bit
            phi_last = 0;
            bit_5 = [bits_trans(jj), bits_trans(jj), bits_trans(jj:jj+2)];
        elseif (jj == 128+24+2)  % 首同步头第2bit
            bit_5 = [bits_trans(jj-1), bits_trans(jj-1:jj+2)];
        elseif (jj == 128+24+128-1)  % 数据部分倒数第2bit
            bit_5 = [bits_trans(jj-2:jj+1), 0];
        elseif (jj == 128+24+128)  % 数据部分最后1bit
            bit_5 = [bits_trans(jj-2:jj), 0, 0];                   

        % S2    
        elseif (jj == 128+24+128+1)  % 首同步头第1bit
            phi_last = 0;
            bit_5 = [bits_trans(jj), bits_trans(jj), bits_trans(jj:jj+2)];
        elseif (jj == 128+24+128+2)  % 首同步头第2bit
            bit_5 = [bits_trans(jj-1), bits_trans(jj-1:jj+2)];
        elseif (jj == 128+24+128+24-1)  % 数据部分倒数第2bit
            bit_5 = [bits_trans(jj-2:jj+1), 0];
        elseif (jj == 128+24+128+24)  % 数据部分最后1bit
            bit_5 = [bits_trans(jj-2:jj), 0, 0];  

        % 尾数据段     
        elseif (jj == 128+128+24+24+1)  % 尾同步头第1bit 
            phi_last = 0;
            bit_5 = [bits_trans(jj),bits_trans(jj),bits_trans(jj:jj+2)];
        elseif (jj == 128+128+24+24+2)  % 尾同步头第2bit 
            bit_5 = [bits_trans(jj-1), bits_trans(jj-1:jj+2)];
        elseif (jj == 128+128+24+24+128-1)  % 数据段倒数第2bit
            bit_5 = [bits_trans(jj-2:jj+1), 0];
        elseif (jj == 128+128+24+24+128)  % 数据段最后1bit
            bit_5 = [bits_trans(jj-2:jj), 0, 0];    

        else
            bit_5 = bits_trans(jj-2:jj+2);
        end

        [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, f_IF, phi_last, g);
        signal_trans_temp_BB(xx,(jj-1)*oversamp_IF+1:(jj)*oversamp_IF) = complex(I_sig(1:16:end), Q_sig(1:16:end))*exp(-j*fp);

    end
    
    for mm = 1:length(signal_trans_temp_BB)-24*4
       xorr(xx,mm)= signal_trans_temp_BB(xx,mm:4:mm+4*24-1)*conj(s1_BB)';
    end
    
%     figure;
%     plot(abs(xorr(xx,:)));    
    
end

% figure;
% plot(real(s1_BB(:)));
% hold on;
% plot(real(signal_trans_temp_BB(128*4+1:4:128*4+24*4)));
% 
% figure;
% plot(real(signal_trans_temp_BB));

xorr_sum=sum(abs(xorr).^2);


figure;
plot(abs(xorr_sum));




