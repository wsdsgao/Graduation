function result = receiver(rx, fh_pat_lib, th_pat_lib, pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, frame_num)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          ��������
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mode_sel = 0;  % ���ڼ�¼�����жϵ�ǰ���ղ���������������ģʽ

num_bits_pn = 24;  % ͬ��ͷS1\S2����
num_bits_pn_2 = 21;  % ͬ��ͷS3\S4����

bit_rate = 16e6;  % ��������
T = 1/bit_rate;  % ����ʱ��
fs_IF = 1024e6;  % ��Ƶ����Ƶ�źŲ�������
fs_BB = 64e6;  % �����źŲ�������
oversamp_BB = T * fs_BB;  % �����źŹ���������
oversamp_IF = T * fs_IF;  % ��Ƶ�źŹ���������

num_bits_pulse = 304;  % 2Mbps A\500Kbps\250Kbps һ��������ĳ���
                       % 2Mbps B ÿ������ȥ��ǰ�������ݲ��ֺ�ĳ���
frame_counter = 1;  % ��¼�Ѽ�⵽����֡���ݣ������ã�
% flag_frame = 0;  % �����һ֡���ñ�־λ��1���ȴ�һ֡ʱ����ٲ�����һ֡�������ã�
% time_counter_mode1 = 6720*oversamp_IF - 0.5*oversamp_IF;  % 2Mbps Aģʽ �ȴ�ʱ��
% time_counter_mode2 = 7956*oversamp_IF - 0.5*oversamp_IF;  % 2Mbps Bģʽ �ȴ�ʱ��
% time_counter_mode3 = 26880*oversamp_IF + 10*oversamp_IF - 0.5*oversamp_IF;  % 500Kbpsģʽ �ȴ�ʱ��
% time_counter_mode4 = 53760*oversamp_IF - 0.5*oversamp_IF;  % 250Kbpsģʽ �ȴ�ʱ��
counter = 0;

flag_Capture_C = 0;  % ָʾ�Ƿ񲶻�һ֡

% ��ģʽһ֡�ܳ���
time_frame_mode1 = (304*12+512*6+103*12);  %ͳһ��������ʱ�͹̶��ļ��
time_frame_mode2 = (304*12+512*6+103*12);
time_frame_mode3 = (304*12+512*6)*4+103*48;
time_frame_mode4 = (304*12+512*6)*8+103*96;

% ��ģʽһ֡��Ч���ݳ���
length_frame_mode1 = 3072;
length_frame_mode2 = 6144;
length_frame_mode3 = 10272;
length_frame_mode4 = 20544;

% ��ģʽһ֡��������
num_pulses_mode1 = 12;
num_pulses_mode2 = 12;
num_pulses_mode3 = 48;
num_pulses_mode4 = 96;

% ��ģʽ��Ƶͼ����ֻ��һ��ͼ������ͬ��
fh_pat_mode1 = fh_pat_lib(1:num_pulses_mode1);
fh_pat_mode2 = fh_pat_lib(1:num_pulses_mode2);
fh_pat_mode3 = fh_pat_lib(1:num_pulses_mode3);
fh_pat_mode4 = fh_pat_lib(1:num_pulses_mode4);

% ��ģʽ��ʱͼ��
th_pat_mode1 = th_pat_lib(1:num_pulses_mode1);
th_pat_mode2 = th_pat_lib(1:num_pulses_mode2);
th_pat_mode3 = th_pat_lib(1:num_pulses_mode3);
th_pat_mode4 = th_pat_lib(1:num_pulses_mode4);


% �����Ѵ�����
load('lib/f_trans.mat');  % 21��Ƶ��
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
load('filter/LPF_1.mat'); % ������1024MHz��ͨ��20MHz
load('filter/LPF_2.mat'); % ������256MHz��ͨ��10MHz
S_lpf1 = 32;
S_lpf2 = 32;
S_bpf = 253;

%���ɱ��ز��Σ���mode4����ͬ������;�ϸУ׼��
[wav_S1_mode1, wav_S2_mode1] = wav_gen(pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 1);

[wav_S1_mode2, wav_S2_mode2] = wav_gen(pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 2);
    
[wav_S1_S3_mode3, wav_S4_S2_mode3] = wav_gen(pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 3);

[wav_S1_S3_mode4, wav_S4_S2_mode4] = wav_gen(pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4, 4);

%����һ��ǰ���Ƿ��貹0���Լ�Ҫ��һ֡�����������˺ü�����ô�죿��

for i = 20:5:200   % �ȴ�һ֡��ʱ�䳤��
    
    %   % �����һ֡���ȴ�һ֡ʱ�䳤���������񣨷����ã�
    %   if (flag_frame)  
    %       if (counter < time_counter - 1)
    %           counter = counter + 1;
    %           continue;
    %       else
    %           counter = 0;
    %           flag_frame = 0;
    %       end
    %   end
    
    
     % δ����ɹ�
     if(~flag_Capture_C)
             
     % for mode1 
     
        % ��ȡ��Ӧһ֡���ȵĲ���
        temp_rx_mode1 = rx(i:i+time_frame_mode1*oversamp_IF-1);
        
        % ������Ƶͼ����Ӧ��Ƶ���ҵ���Ӧͨ���Ĳ���
        wav_temp_mode1 = zeros(num_pulses_mode1, num_bits_pulse * oversamp_IF);
        % ��1����ʱ��Ƶͼ��
        for pulse_idx = 1:num_pulses_mode1
            f_idx = fh_pat_mode1(pulse_idx); % ��ǰ�����ӦƵ��
            
            % ���㵱ǰ������һ֡�е���ʼλ��
            if pulse_idx == 1
                pos_pulse_mode1 = th_pat_mode1(1)/2;
            else
                pos_pulse_mode1 = sum(th_pat_mode1(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103) + th_pat_mode1(pulse_idx)/2;
            end

            % ��Ƶ -> ��Ƶ ��ֻȡ������λ��
            if (f_idx >= 1) && (f_idx <= 5)
                wav_temp1_mode1 = downConv_IF(temp_rx_mode1(pos_pulse_mode1*oversamp_IF+1:pos_pulse_mode1*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % ͨ��1
            elseif (f_idx >= 6) && (f_idx <= 10)
                wav_temp1_mode1 = downConv_IF(temp_rx_mode1(pos_pulse_mode1*oversamp_IF+1:pos_pulse_mode1*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % ͨ��2
            elseif (f_idx >= 11) && (f_idx <= 14)
                wav_temp1_mode1 = downConv_IF(temp_rx_mode1(pos_pulse_mode1*oversamp_IF+1:pos_pulse_mode1*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % ͨ��3
            elseif (f_idx >= 15) && (f_idx <= 18)
                wav_temp1_mode1 = downConv_IF(temp_rx_mode1(pos_pulse_mode1*oversamp_IF+1:pos_pulse_mode1*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % ͨ��4 
            elseif (f_idx >= 19) && (f_idx <= 21)
                wav_temp1_mode1 = downConv_IF(temp_rx_mode1(pos_pulse_mode1*oversamp_IF+1:pos_pulse_mode1*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % ͨ��5 
            end

            wav_temp_mode1(pulse_idx,:) = wav_temp1_mode1;

        end


        % ������Ƶͼ����Ӧ��Ƶ�㽫��Ӧͨ���Ĳ����˲�������Ƶ
        % wav_temp(������1024) -> rx_pulse_mat(������64)
        rx_pulse_mat_mode1 = downConv(wav_temp_mode1(:,:), num_pulses_mode1, fh_pat_mode1);  % ��Ӧǰ10ms
        % rx_pulse_mat_2_mode1 = downConv(wav_temp_mode1(:,:,2), num_pulses_mode1, fh_pat_2_mode1);  % ��Ӧ��10ms

        % plot(real(rx_pulse_mat_mode1(1, :)));
        % value = sum(abs(rx_pulse_mat_mode1(1, :)))


        % % Ԥȡǰ��24bitλ�õĲ��Σ���ͬ��ͷ�� �ȴ���������
        % ǰ10ms
        D_S1_mode1 = rx_pulse_mat_mode1(:,1:24*oversamp_BB);  %һ����һ������
        D_S1_one_mode1 = D_S1_mode1(:,4:oversamp_BB:end);  %������������ͬ��һ�����Ӧһ��bit
        D_S2_mode1 = rx_pulse_mat_mode1(:,1120+1:1120+24*oversamp_BB);  %2240=280*oversamp_BB
        D_S2_one_mode1 = D_S2_mode1(:,4:oversamp_BB:end);

        % % ��10ms
        % D_S1_2_mode1 = rx_pulse_mat_2_mode1(:,1:24*oversamp_BB);
        % D_S1_2_one_mode1 = D_S1_2_mode1(:,8:oversamp_BB:end);
        % D_S2_2_mode1 = rx_pulse_mat_2_mode1(:,2240+1:2240+24*oversamp_BB);
        % D_S2_2_one_mode1 = D_S2_2_mode1(:,8:oversamp_BB:end);

        % ͬ��ͷ�������                 
        % 1���������ʲ�������������
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% 
        %
        %                           ǰ10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%             

        %ԭ�㷨
        rx_corr_S1_pat_mode1 = zeros(1,num_pulses_mode1);
        rx_corr_S2_pat_mode1 = zeros(1,num_pulses_mode1);
        for j = 1:num_pulses_mode1

            rx_wav_S1_pat_mode1 = D_S1_one_mode1(j,1:22) .* conj(wav_S1_mode1(j,1:22));
            rx_wav_S2_pat_mode1 = D_S2_one_mode1(j,1:22) .* conj(wav_S2_mode1(j,1:22));   

            rx_corr_S1_pat_mode1(j) = abs(sum(rx_wav_S1_pat_mode1));
            rx_corr_S2_pat_mode1(j) = abs(sum(rx_wav_S2_pat_mode1));

            rx_corr_pat_mode1(j) =  rx_corr_S1_pat_mode1(j) + rx_corr_S2_pat_mode1(j);

        end

        corr_value_mode1 = sum(rx_corr_pat_mode1);

        corr_value_mode1_plot1((i-20)/5+1) = corr_value_mode1;

        % ��ǰ��ͬ��ͷ��ط�(�Ľ��㷨)
        corr_value_mode1_plot2((i-20)/5+1) = corr(rx_pulse_mat_mode1, wav_S1_S3_mode4, wav_S4_S2_mode4, 1);
        % rx_corr_S1_pat_mode1 = zeros(1,num_pulses_mode1);
        % rx_corr_S2_pat_mode1 = zeros(1,num_pulses_mode1);
        % rx_corr_pat_mode1 = zeros(1,num_pulses_mode1);
        % for j = 1:num_pulses_mode1

        %     rx_wav_S1_pat_mode1_1 = D_S1_one_mode1(j,1:8) .* (wav_S1_S3_mode4(j,1:8));
        %     rx_wav_S1_pat_mode1_2 = D_S1_one_mode1(j,9:16) .* conj(wav_S1_S3_mode4(j,9:16));
        %     rx_wav_S1_pat_mode1_3 = D_S1_one_mode1(j,17:24) .* conj(wav_S1_S3_mode4(j,17:24)); %����һ���Ƿ�Ϊ1-22�ȽϺ���
        %     rx_wav_S2_pat_mode1_1 = D_S2_one_mode1(j,1:8) .* (wav_S4_S2_mode4(j,22:29)); %����һ���Ƿ�Ϊ24-45�ȽϺ���
        %     rx_wav_S2_pat_mode1_2 = D_S2_one_mode1(j,9:16) .* conj(wav_S4_S2_mode4(j,30:37));
        %     rx_wav_S2_pat_mode1_3 = D_S2_one_mode1(j,17:24) .* conj(wav_S4_S2_mode4(j,38:45));   

        %     rx_corr_S1_pat_mode1(j) = abs(sum(rx_wav_S1_pat_mode1_1)+sum(rx_wav_S1_pat_mode1_2)+sum(rx_wav_S1_pat_mode1_3));
        %     rx_corr_S2_pat_mode1(j) = abs(sum(rx_wav_S2_pat_mode1_1)+sum(rx_wav_S2_pat_mode1_2)+sum(rx_wav_S2_pat_mode1_3));

        %     rx_corr_pat_mode1(j) =  rx_corr_S1_pat_mode1(j) + rx_corr_S2_pat_mode1(j);

        % end

        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                           ��10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % ��ǰ��ͬ��ͷ��ط�
        % rx_corr_S1_pat_2_mode1 = zeros(1,num_pulses_mode1);  
        % rx_corr_S2_pat_2_mode1 = zeros(1,num_pulses_mode1);
        % for j = 1:num_pulses_mode1

        %     rx_wav_S1_pat_2_mode1 = D_S1_2_one_mode1(j,1:22) .* conj(wav_S1_2_mode1(j,1:22));
        %     rx_wav_S2_pat_2_mode1 = D_S2_2_one_mode1(j,1:22) .* conj(wav_S2_2_mode1(j,1:22));   

        %     rx_corr_S1_pat_2_mode1(j) = abs(sum(rx_wav_S1_pat_2_mode1));  %��һ�еĺ�
        %     rx_corr_S2_pat_2_mode1(j) = abs(sum(rx_wav_S2_pat_2_mode1));

        % end
        
        % ͬ�������о�����
        % ǰ(���)ͬ��ͷ��ط�ֵ������ֵ
        if (corr_value_mode1 > 44*num_pulses_mode1) %�о�����Ҫ���ź�����ȥ�ȣ���û��
            mode_sel = 1;
            flag_Capture_C = 1;
            % if (sum(rx_corr_S1_pat_1_mode1) > sum(rx_corr_S2_pat_1_mode1))
            %     tag = 1;
            % else
            %     tag = 2;
            % end

        % elseif (sum(rx_corr_S1_pat_2_mode1)>20*num_pulses_mode1) || (sum(rx_corr_S2_pat_2_mode1)>20*num_pulses_mode1)
        %     mode_sel = 1;
        %     flag_Capture_C = 1;
        %     if (sum(rx_corr_S1_pat_2_mode1) > sum(rx_corr_S2_pat_2_mode1))
        %         tag = 3;
        %     else
        %         tag = 4;
        %     end
        end        

        
    % for mode2
    
        % ��ȡ��Ӧһ֡���ȵĲ���
        temp_rx_mode2 = rx(i:i+time_frame_mode2*oversamp_IF-1);
        % ������Ƶͼ����Ӧ��Ƶ���ҵ���Ӧͨ���Ĳ���
        % ��������10ms��ͼ��
        wav_temp_mode2 = zeros(num_pulses_mode2, num_bits_pulse * oversamp_IF);  %û������ʱ������
        % ��1����ʱ��Ƶͼ��
        for pulse_idx = 1:num_pulses_mode2
            f_idx = fh_pat_mode2(pulse_idx);  % ��ǰ�����Ӧ��Ƶ��

            % ���㵱ǰ������һ֡�е���ʼλ��
            if pulse_idx == 1
                pos_pulse_mode2 = 0;
            else
                pos_pulse_mode2 = sum(th_pat_mode2(1:pulse_idx-1)) + 103 * (pulse_idx-1) + (pulse_idx-1)*num_bits_pulse;
            end  

            % ���㵱ǰ����ĳ���
            num_bits_pulse_rx = th_pat_mode2(pulse_idx) + num_bits_pulse; %��������ʱ

            % ��Ƶ -> ��Ƶ
            if (f_idx >= 1) && (f_idx <= 5)    
                wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % ͨ��1
            elseif (f_idx >= 6) && (f_idx <= 10)
                wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % ͨ��2
            elseif (f_idx >= 11) && (f_idx <= 14)
                wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % ͨ��3
            elseif (f_idx >= 15) && (f_idx <= 18)
                wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % ͨ��4 
            elseif (f_idx >= 19) && (f_idx <= 21)
                wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % ͨ��5 
            end

            wav_temp_mode2(pulse_idx,:) = wav_temp1_mode2((th_pat_mode2(pulse_idx)/2)*oversamp_IF+1:(th_pat_mode2(pulse_idx)/2+num_bits_pulse)*oversamp_IF);  %û������ʱ������
        end


        % % ��2����ʱ��Ƶͼ��
        % for pulse_idx = 1:num_pulses_mode2
        %     f_idx = fh_pat_2_mode2(pulse_idx);  % ��ǰ�����Ӧ��Ƶ��

        %     % ���㵱ǰ������һ֡�е���ʼλ��
        %     if pulse_idx == 1
        %         pos_pulse_mode2 = 100;
        %     else
        %         pos_pulse_mode2 = sum(th_pat_2_mode2(1:pulse_idx-1)) + 3 * (pulse_idx-1) + (pulse_idx-1)*num_bits_pulse + 100 * pulse_idx;
        %     end
            
        %     % ���㵱ǰ���峤��
        %     num_bits_pulse_rx = th_pat_2_mode2(pulse_idx) + 3 + num_bits_pulse;

        %     % ��Ƶ -> ��Ƶ
        %     if (f_idx >= 1) && (f_idx <= 5)    
        %         wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % ͨ��1
        %     elseif (f_idx >= 6) && (f_idx <= 10)
        %         wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % ͨ��2
        %     elseif (f_idx >= 11) && (f_idx <= 14)
        %         wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % ͨ��3
        %     elseif (f_idx >= 15) && (f_idx <= 18)
        %         wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % ͨ��4 
        %     elseif (f_idx >= 19) && (f_idx <= 21)
        %         wav_temp1_mode2 = downConv_IF(temp_rx_mode2(pos_pulse_mode2*oversamp_IF+1:pos_pulse_mode2*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % ͨ��5 
        %     end

        %     wav_temp_mode2(pulse_idx,:,2) = wav_temp1_mode2((th_pat_2_mode2(pulse_idx)/2+3)*oversamp_IF+1:(th_pat_2_mode2(pulse_idx)/2+3+num_bits_pulse)*oversamp_IF);
        % end


        % ������Ƶͼ����Ӧ��Ƶ�㽫��Ӧͨ���Ĳ����˲�������Ƶ
        % wav_temp(������1024) -> rx_pulse_mat(������64)
        rx_pulse_mat_mode2 = downConv_2MLL(wav_temp_mode2(:,:), num_pulses_mode2, fh_pat_mode2, th_pat_mode2);  % ��Ӧǰ10ms
        % rx_pulse_mat_2_mode2 = downConv_2MLL(wav_temp_mode2(:,:,2), num_pulses_mode2, fh_pat_2_mode2, th_pat_2_mode2);  % ��Ӧ��10ms


        % % Ԥȡǰ��24bitλ�õĲ��� �ȴ���������
        % % ǰ10ms
        % D_S1_mode2 = rx_pulse_mat_mode2(:,1:num_bits_pn*oversamp_BB);
        % D_S1_one_mode2 = D_S1_mode2(:,4:oversamp_BB:end);
        % D_S2_mode2 = rx_pulse_mat_mode2(:,2240+1:2240+num_bits_pn*oversamp_BB);
        % D_S2_one_mode2 = D_S2_mode2(:,4:oversamp_BB:end);
        % % ��10ms
        % D_S1_2_mode2 = rx_pulse_mat_2_mode2(:,1:num_bits_pn*oversamp_BB);
        % D_S1_2_one_mode2 = D_S1_2_mode2(:,8:oversamp_BB:end);
        % D_S2_2_mode2 = rx_pulse_mat_2_mode2(:,2240+1:2240+num_bits_pn*oversamp_BB);
        % D_S2_2_one_mode2 = D_S2_2_mode2(:,8:oversamp_BB:end);


        % ͬ��ͷ�������                 
        % 1���������ʲ�������������
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                   ǰ10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % ��ǰ��ͬ��ͷ��ط�(��Ƶƫ����ƫ)
        corr_value_mode2 = corr(rx_pulse_mat_mode2, wav_S1_S3_mode4, wav_S4_S2_mode4, 2);
        % rx_corr_S1_pat_1_mode2 = zeros(1,num_pulses_mode2);
        % rx_corr_S2_pat_1_mode2 = zeros(1,num_pulses_mode2);
        % for j = 1:num_pulses_mode2

        %     rx_wav_S1_pat_1 = D_S1_1_one_mode2(j,1:22) .* conj(wav_S1_1_mode2(j,1:22));
        %     rx_wav_S2_pat_1 = D_S2_1_one_mode2(j,1:22) .* conj(wav_S2_1_mode2(j,1:22));   

        %     rx_corr_S1_pat_1_mode2(j) = abs(sum(rx_wav_S1_pat_1));
        %     rx_corr_S2_pat_1_mode2(j) = abs(sum(rx_wav_S2_pat_1));

        % end


        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                       ��10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % % ��ǰ��ͬ��ͷ��ط�(��Ƶƫ����ƫ)
        % rx_corr_S1_pat_2_mode2 = zeros(1,num_pulses_mode2);  
        % rx_corr_S2_pat_2_mode2 = zeros(1,num_pulses_mode2);
        % for j = 1:num_pulses_mode2

        %     rx_wav_S1_pat_2 = D_S1_2_one_mode2(j,1:22) .* conj(wav_S1_2_mode2(j,1:22));
        %     rx_wav_S2_pat_2 = D_S2_2_one_mode2(j,1:22) .* conj(wav_S2_2_mode2(j,1:22));   

        %     rx_corr_S1_pat_2_mode2(j) = abs(sum(rx_wav_S1_pat_2));
        %     rx_corr_S2_pat_2_mode2(j) = abs(sum(rx_wav_S2_pat_2));

        % end


         % ͬ�������о�����
         % ǰ(���)ͬ��ͷ��ط�ֵ������ֵ       
        if (corr_value_mode2 > 44*num_pulses_mode2) 
            mode_sel = 2;
            flag_Capture_C = 1; 
            % if (sum(rx_corr_S1_pat_1_mode2) > sum(rx_corr_S2_pat_1_mode2))
            %     tag = 1;
            % else
            %     tag = 2;
            % end
        %  elseif (sum(rx_corr_S1_pat_2_mode2)>20*num_pulses_mode2) || (sum(rx_corr_S2_pat_2_mode2)>20*num_pulses_mode2)
        %     mode_sel = 2;
        %     flag_Capture_C = 1;
        %     if (sum(rx_corr_S1_pat_2_mode2) > sum(rx_corr_S2_pat_2_mode2))
        %         tag = 3;
        %     else
        %         tag = 4;
        %     end
        end 

    % for mode3

        if (i+time_frame_mode3*oversamp_IF-1) > length(rx)
            continue;
        end
    
        % ��ȡ��Ӧ48�����峤�ȵĲ���
        temp_rx_mode3 = rx(i:i+time_frame_mode3*oversamp_IF-1);
        
        % ������Ƶͼ����Ӧ��Ƶ���ҵ���Ӧͨ���Ĳ���
        % ��������10ms��ͼ��
        wav_temp_mode3 = zeros(num_pulses_mode3, num_bits_pulse * oversamp_IF);
        % ��1����ʱ��Ƶͼ��
        for pulse_idx = 1:num_pulses_mode3
            f_idx = fh_pat_mode3(pulse_idx);  % ��ǰ�����Ӧ��Ƶ��
            
            % ���㵱ǰ������һ֡�е���ʼλ��
            if pulse_idx == 1
                pos_pulse_mode3 = th_pat_mode3(1)/2;
            else
                pos_pulse_mode3 = sum(th_pat_mode3(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103) + th_pat_mode3(pulse_idx)/2;
            end

            % ��Ƶ -> ��Ƶ
            if (f_idx >= 1) && (f_idx <= 5)    
                wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % ͨ��1
            elseif (f_idx >= 6) && (f_idx <= 10)
                wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % ͨ��2
            elseif (f_idx >= 11) && (f_idx <= 14)
                wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % ͨ��3
            elseif (f_idx >= 15) && (f_idx <= 18)
                wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % ͨ��4 
            elseif (f_idx >= 19) && (f_idx <= 21)
                wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % ͨ��5 
            end

            wav_temp_mode3(pulse_idx,:) = wav_temp1_mode3;

        end


        % ��2����ʱ��Ƶͼ��
        % for pulse_idx = 1:num_pulses_mode3
        %     f_idx = fh_pat_2_mode3(pulse_idx);  % ��ǰ�����Ӧ��Ƶ��
            
        %     % ���㵱ǰ������һ֡�е���ʼλ��
        %     if pulse_idx == 1
        %         pos_pulse_mode3 = th_pat_2_mode3(1)/2;
        %     else
        %         pos_pulse_mode3 = sum(th_pat_2_mode3(1:pulse_idx-1)) + (pulse_idx-1)*num_bits_pulse + th_pat_2_mode3(pulse_idx)/2;
        %     end

        %     % ��Ƶ -> ��Ƶ
        %     if (f_idx >= 1) && (f_idx <= 5)    
        %         wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % ͨ��1
        %     elseif (f_idx >= 6) && (f_idx <= 10)
        %         wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % ͨ��2
        %     elseif (f_idx >= 11) && (f_idx <= 14)
        %         wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % ͨ��3
        %     elseif (f_idx >= 15) && (f_idx <= 18)
        %         wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % ͨ��4 
        %     elseif (f_idx >= 19) && (f_idx <= 21)
        %         wav_temp1_mode3 = downConv_IF(temp_rx_mode3(pos_pulse_mode3*oversamp_IF+1:pos_pulse_mode3*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % ͨ��5 
        %     end

        %     wav_temp_mode3(pulse_idx,:,2) = wav_temp1_mode3;

        % end



        % ������Ƶͼ����Ӧ��Ƶ�㽫��Ӧͨ���Ĳ����˲�������Ƶ
        % wav_temp(������1024) -> rx_pulse_mat(������64)
        rx_pulse_mat_mode3 = downConv(wav_temp_mode3(:,:), num_pulses_mode3, fh_pat_mode3);  % ��Ӧǰ10ms
        % rx_pulse_mat_2_mode3 = downConv_500K(wav_temp_mode3(:,:,2), num_pulses_mode3, fh_pat_2_mode3);  % ��Ӧ��10ms


        % % Ԥȡǰ��24bitλ�õĲ��� �ȴ���������
        % % ǰ10ms
        % D_S1_mode3 = rx_pulse_mat_mode3(:,1:24*oversamp_BB);
        % D_S1_one_mode3 = D_S1_mode3(:,4:oversamp_BB:end);
        % D_S2_mode3 = rx_pulse_mat_mode3(:,280*oversamp_BB+1:304*oversamp_BB);
        % D_S2_one_mode3 = D_S2_mode3(:,4:oversamp_BB:end);
        % D_S3_mode3 = rx_pulse_mat_mode3(:,24*oversamp_BB+1:45*oversamp_BB);
        % D_S3_one_mode3 = D_S3_mode3(:,4:oversamp_BB:end);
        % D_S4_mode3 = rx_pulse_mat_mode3(:,259*oversamp_BB+1:280*oversamp_BB);
        % D_S4_one_mode3 = D_S4_mode3(:,4:oversamp_BB:end);

        % % ��10ms
        % D_S1_2_mode3 = rx_pulse_mat_2_mode3(:,1:24*oversamp_BB);
        % D_S1_2_one_mode3 = D_S1_2_mode3(:,8:oversamp_BB:end);
        % D_S2_2_mode3 = rx_pulse_mat_2_mode3(:,280*oversamp_BB+1:304*oversamp_BB);
        % D_S2_2_one_mode3 = D_S2_2_mode3(:,8:oversamp_BB:end);
        % D_S3_2_mode3 = rx_pulse_mat_2_mode3(:,24*oversamp_BB+1:45*oversamp_BB);
        % D_S3_2_one_mode3 = D_S3_2_mode3(:,8:oversamp_BB:end);
        % D_S4_2_mode3 = rx_pulse_mat_2_mode3(:,259*oversamp_BB+1:280*oversamp_BB);
        % D_S4_2_one_mode3 = D_S4_2_mode3(:,8:oversamp_BB:end);            


        % ͬ��ͷ�������    (S1_S3____S4_S2)             
        % 1���������ʲ�������������
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                   ǰ10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % ��ǰ��ͬ��ͷ��ط�
        corr_value_mode3 = corr(rx_pulse_mat_mode3, wav_S1_S3_mode4, wav_S4_S2_mode4, 3);
        % rx_corr_S1_pat_1_mode3 = zeros(1,num_pulses_mode3);
        % rx_corr_S2_pat_1_mode3 = zeros(1,num_pulses_mode3);
        % rx_corr_S3_pat_1_mode3 = zeros(1,num_pulses_mode3);
        % rx_corr_S4_pat_1_mode3 = zeros(1,num_pulses_mode3);
        % for j = 1:num_pulses_mode3

        %     rx_wav_S1_pat_1 = D_S1_1_one_mode3(j,1:24) .* conj(wav_S1_S3_1_mode3(j,1:24));     % S1 24bit
        %     rx_wav_S2_pat_1 = D_S2_1_one_mode3(j,1:22) .* conj(wav_S4_S2_1_mode3(j,21+1:43));  % S2 22bit  
        %     rx_wav_S3_pat_1 = D_S3_1_one_mode3(j,1:19) .* conj(wav_S1_S3_1_mode3(j,24+1:43));  % S3 19bit
        %     rx_wav_S4_pat_1 = D_S4_1_one_mode3(j,1:21) .* conj(wav_S4_S2_1_mode3(j,1:21));     % S4 21bit

        %     rx_corr_S1_pat_1_mode3(j) = abs(sum(rx_wav_S1_pat_1));
        %     rx_corr_S2_pat_1_mode3(j) = abs(sum(rx_wav_S2_pat_1));
        %     rx_corr_S3_pat_1_mode3(j) = abs(sum(rx_wav_S3_pat_1));
        %     rx_corr_S4_pat_1_mode3(j) = abs(sum(rx_wav_S4_pat_1));

        % end

        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                       ��10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % rx_corr_S1_pat_2_mode3 = zeros(1,num_pulses_mode3);  
        % rx_corr_S2_pat_2_mode3 = zeros(1,num_pulses_mode3);
        % rx_corr_S3_pat_2_mode3 = zeros(1,num_pulses_mode3);
        % rx_corr_S4_pat_2_mode3 = zeros(1,num_pulses_mode3);
        % for j = 1:num_pulses_mode3

        %     rx_wav_S1_pat_2 = D_S1_2_one_mode3(j,1:24) .* conj(wav_S1_S3_2_mode3(j,1:24));     % S1 24bit
        %     rx_wav_S2_pat_2 = D_S2_2_one_mode3(j,1:22) .* conj(wav_S4_S2_2_mode3(j,21+1:43));  % S2 22bit  
        %     rx_wav_S3_pat_2 = D_S3_2_one_mode3(j,1:19) .* conj(wav_S1_S3_2_mode3(j,24+1:43));  % S3 19bit
        %     rx_wav_S4_pat_2 = D_S4_2_one_mode3(j,1:21) .* conj(wav_S4_S2_2_mode3(j,1:21));     % S4 21bit
            
        %     rx_corr_S1_pat_2_mode3(j) = abs(sum(rx_wav_S1_pat_2));
        %     rx_corr_S2_pat_2_mode3(j) = abs(sum(rx_wav_S2_pat_2));
        %     rx_corr_S3_pat_2_mode3(j) = abs(sum(rx_wav_S3_pat_2));
        %     rx_corr_S4_pat_2_mode3(j) = abs(sum(rx_wav_S4_pat_2));

        % end


         % ͬ�������о�����
         % ǰ(���)ͬ��ͷ�����������12������ƽ��ֵ������4        
         if (corr_value_mode3 > 82*num_pulses_mode3)
             flag_Capture_C = 1;
             mode_sel = 3;
             
        %      if (sum(rx_corr_S1_pat_1_mode3) + sum(rx_corr_S3_pat_1_mode3) > sum(rx_corr_S2_pat_1_mode3) + sum(rx_corr_S4_pat_1_mode3))
        %          tag = 1
        %      else
        %          tag = 2
        %      end

        %  elseif ((sum(rx_corr_S1_pat_2_mode3)>22*num_pulses_mode3) && (sum(rx_corr_S3_pat_2_mode3)>17*num_pulses_mode3)) ...
        %          || ((sum(rx_corr_S2_pat_2_mode3)>20*num_pulses_mode3) && (sum(rx_corr_S4_pat_2_mode3)>19*num_pulses_mode3))
        %      flag_Capture_C = 1;
        %      mode_sel = 3;
        %      if (sum(rx_corr_S1_pat_2_mode3) + sum(rx_corr_S3_pat_2_mode3) > sum(rx_corr_S2_pat_2_mode3) + sum(rx_corr_S4_pat_2_mode3))
        %         tag = 3         
        %      else
        %         tag = 4
        %      end

         end         
    

    % for mode4
    
        if (i+time_frame_mode4*oversamp_IF-1) > length(rx)
            continue;
        end

        % ��ȡ��Ӧ96�����峤�ȵĲ���
        temp_rx_mode4 = rx(i:i+time_frame_mode4*oversamp_IF-1);
        
        % ������Ƶͼ����Ӧ��Ƶ���ҵ���Ӧͨ���Ĳ���
        % ��������10ms��ͼ��
        wav_temp_mode4 = zeros(num_pulses_mode4, num_bits_pulse * oversamp_IF);
        % ��1����ʱ��Ƶͼ��
        for pulse_idx = 1:num_pulses_mode4
            f_idx = fh_pat_mode4(pulse_idx);  % ��ǰ�����Ӧ��Ƶ��
            
            % ���㵱ǰ������һ֡�е���ʼλ��
            if pulse_idx == 1
                pos_pulse_mode4 = th_pat_mode4(1)/2;
            else
                pos_pulse_mode4 = sum(th_pat_mode4(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103) + th_pat_mode4(pulse_idx)/2;
            end
            
            % ��Ƶ -> ��Ƶ
            if (f_idx >= 1) && (f_idx <= 5)    
                wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % ͨ��1
            elseif (f_idx >= 6) && (f_idx <= 10)
                wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % ͨ��2
            elseif (f_idx >= 11) && (f_idx <= 14)
                wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % ͨ��3
            elseif (f_idx >= 15) && (f_idx <= 18)
                wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % ͨ��4 
            elseif (f_idx >= 19) && (f_idx <= 21)
                wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % ͨ��5 
            end

            wav_temp_mode4(pulse_idx,:) = wav_temp1_mode4;

        end


        % % ��2����ʱ��Ƶͼ��
        % for pulse_idx = 1:num_pulses_mode4
        %     f_idx = fh_pat_2_mode4(pulse_idx);  % ��ǰ����Ķ�ӦƵ��
            
        %     % ���㵱ǰ������һ֡�е���ʼλ��
        %     if pulse_idx == 1
        %         pos_pulse_mode4 = th_pat_2_mode4(1)/2;
        %     else
        %         pos_pulse_mode4 = sum(th_pat_2_mode4(1:pulse_idx-1)) + (pulse_idx-1)*num_bits_pulse + th_pat_2_mode4(pulse_idx)/2;
        %     end
            
        %     % ��Ƶ -> ��Ƶ
        %     if (f_idx >= 1) && (f_idx <= 5)    
        %         wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % ͨ��1
        %     elseif (f_idx >= 6) && (f_idx <= 10)
        %         wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % ͨ��2
        %     elseif (f_idx >= 11) && (f_idx <= 14)
        %         wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % ͨ��3
        %     elseif (f_idx >= 15) && (f_idx <= 18)
        %         wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % ͨ��4 
        %     elseif (f_idx >= 19) && (f_idx <= 21)
        %         wav_temp1_mode4 = downConv_IF(temp_rx_mode4(pos_pulse_mode4*oversamp_IF+1:pos_pulse_mode4*oversamp_IF+num_bits_pulse*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % ͨ��5 
        %     end

        %     wav_temp_mode4(pulse_idx,:,2) = wav_temp1_mode4;

        % end



        % ������Ƶͼ����Ӧ��Ƶ�㽫��Ӧͨ���Ĳ����˲�������Ƶ
        % wav_temp(������1024) -> rx_pulse_mat(������64)
        rx_pulse_mat_mode4 = downConv(wav_temp_mode4(:,:), num_pulses_mode4, fh_pat_mode4);  % ��Ӧǰ10ms
        % rx_pulse_mat_2_mode4 = downConv_500K(wav_temp_mode4(:,:,2), num_pulses_mode4, fh_pat_2_mode4);  % ��Ӧ��10ms


        % % Ԥȡǰ��24bitλ�õĲ��� �ȴ���������
        % % ǰ10ms
        % D_S1_mode4 = rx_pulse_mat_mode4(:,1:24*oversamp_BB);
        % D_S1_one_mode4 = D_S1_mode4(:,4:oversamp_BB:end);
        % D_S2_mode4 = rx_pulse_mat_mode4(:,280*oversamp_BB+1:304*oversamp_BB);
        % D_S2_one_mode4 = D_S2_mode4(:,4:oversamp_BB:end);
        % D_S3_mode4 = rx_pulse_mat_mode4(:,24*oversamp_BB+1:45*oversamp_BB);
        % D_S3_one_mode4 = D_S3_mode4(:,4:oversamp_BB:end);
        % D_S4_mode4 = rx_pulse_mat_mode4(:,259*oversamp_BB+1:280*oversamp_BB);
        % D_S4_one_mode4 = D_S4_mode4(:,4:oversamp_BB:end);

        % % ��10ms
        % D_S1_2_mode4 = rx_pulse_mat_2_mode4(:,1:24*oversamp_BB);
        % D_S1_2_one_mode4 = D_S1_2_mode4(:,8:oversamp_BB:end);
        % D_S2_2_mode4 = rx_pulse_mat_2_mode4(:,280*oversamp_BB+1:304*oversamp_BB);
        % D_S2_2_one_mode4 = D_S2_2_mode4(:,8:oversamp_BB:end);
        % D_S3_2_mode4 = rx_pulse_mat_2_mode4(:,24*oversamp_BB+1:45*oversamp_BB);
        % D_S3_2_one_mode4 = D_S3_2_mode4(:,8:oversamp_BB:end);
        % D_S4_2_mode4 = rx_pulse_mat_2_mode4(:,259*oversamp_BB+1:280*oversamp_BB);
        % D_S4_2_one_mode4 = D_S4_2_mode4(:,8:oversamp_BB:end);            


        % ͬ��ͷ�������    (S1_S3____S4_S2)             
        % 1���������ʲ�������������
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                   ǰ10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % ��ǰ��ͬ��ͷ��ط�
        corr_value_mode4 = corr(rx_pulse_mat_mode4, wav_S1_S3_mode4, wav_S4_S2_mode4, 4);
        % rx_corr_S1_pat_1_mode4 = zeros(1,num_pulses_mode4);
        % rx_corr_S2_pat_1_mode4 = zeros(1,num_pulses_mode4);
        % rx_corr_S3_pat_1_mode4 = zeros(1,num_pulses_mode4);
        % rx_corr_S4_pat_1_mode4 = zeros(1,num_pulses_mode4);
        % for j = 1:num_pulses_mode4

        %     rx_wav_S1_pat_1 = D_S1_1_one_mode4(j,1:24) .* conj(wav_S1_S3_1_mode4(j,1:24));     % S1 24bit
        %     rx_wav_S2_pat_1 = D_S2_1_one_mode4(j,1:22) .* conj(wav_S4_S2_1_mode4(j,21+1:43));  % S2 22bit  
        %     rx_wav_S3_pat_1 = D_S3_1_one_mode4(j,1:19) .* conj(wav_S1_S3_1_mode4(j,24+1:43));  % S3 19bit
        %     rx_wav_S4_pat_1 = D_S4_1_one_mode4(j,1:21) .* conj(wav_S4_S2_1_mode4(j,1:21));     % S4 21bit

        %     rx_corr_S1_pat_1_mode4(j) = abs(sum(rx_wav_S1_pat_1));
        %     rx_corr_S2_pat_1_mode4(j) = abs(sum(rx_wav_S2_pat_1));
        %     rx_corr_S3_pat_1_mode4(j) = abs(sum(rx_wav_S3_pat_1));
        %     rx_corr_S4_pat_1_mode4(j) = abs(sum(rx_wav_S4_pat_1));

        % end

        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%
        %
        %                                       ��10ms
        %
        %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%

        % % ��ǰ��ͬ��ͷ��ط�
        % rx_corr_S1_pat_2_mode4 = zeros(1,num_pulses_mode4);  
        % rx_corr_S2_pat_2_mode4 = zeros(1,num_pulses_mode4);
        % rx_corr_S3_pat_2_mode4 = zeros(1,num_pulses_mode4);
        % rx_corr_S4_pat_2_mode4 = zeros(1,num_pulses_mode4);
        % for j = 1:num_pulses_mode4

        %     rx_wav_S1_pat_2 = D_S1_2_one_mode4(j,1:24) .* conj(wav_S1_S3_2_mode4(j,1:24));     % S1 24bit
        %     rx_wav_S2_pat_2 = D_S2_2_one_mode4(j,1:22) .* conj(wav_S4_S2_2_mode4(j,21+1:43));  % S2 22bit  
        %     rx_wav_S3_pat_2 = D_S3_2_one_mode4(j,1:19) .* conj(wav_S1_S3_2_mode4(j,24+1:43));  % S3 19bit
        %     rx_wav_S4_pat_2 = D_S4_2_one_mode4(j,1:21) .* conj(wav_S4_S2_2_mode4(j,1:21));     % S4 21bit
            
        %     rx_corr_S1_pat_2_mode4(j) = abs(sum(rx_wav_S1_pat_2));
        %     rx_corr_S2_pat_2_mode4(j) = abs(sum(rx_wav_S2_pat_2));
        %     rx_corr_S3_pat_2_mode4(j) = abs(sum(rx_wav_S3_pat_2));
        %     rx_corr_S4_pat_2_mode4(j) = abs(sum(rx_wav_S4_pat_2));

        % end


         % ͬ�������о�����
        
         if (corr_value_mode4 > 82*num_pulses_mode4) 
             flag_Capture_C = 1;
             mode_sel = 4;
             
        %      if (sum(rx_corr_S1_pat_1_mode4) + sum(rx_corr_S3_pat_1_mode4) > sum(rx_corr_S2_pat_1_mode4) + sum(rx_corr_S4_pat_1_mode4))
        %          tag = 1
        %      else
        %          tag = 2
        %      end
          

        %  elseif ((sum(rx_corr_S1_pat_2_mode4)>22*num_pulses_mode4) && (sum(rx_corr_S3_pat_2_mode4)>17*num_pulses_mode4)) ...
        %          || ((sum(rx_corr_S2_pat_2_mode4)>20*num_pulses_mode4) && (sum(rx_corr_S4_pat_2_mode4)>19*num_pulses_mode4))
        %      flag_Capture_C = 1;
        %      mode_sel = 4;
             
        %      if (sum(rx_corr_S1_pat_2_mode4) + sum(rx_corr_S3_pat_2_mode4) > sum(rx_corr_S2_pat_2_mode4) + sum(rx_corr_S4_pat_2_mode4))
        %         tag = 3         
        %      else
        %         tag = 4
        %      end

        end
   
        disp(i);
        % figure(1);
        % plot(real(D_S1_one_mode1(1,:)));
        % figure(2);
        % plot(real(wav_S1_mode1(1,:)));

    end % (~flag_Capture_C)

    % ����ɹ�
    % �����ж�������ģʽ�����Ӧ�Ľ��ģ��
    if (flag_Capture_C == 1)
        if mode_sel == 1
            num_pulses_rec = num_pulses_mode1;  % һ֡����������
            length_frame = length_frame_mode1;  % һ֡����Ч���ݳ���
            % time_counter = time_counter_mode1;  % �����һ֡�ĵȴ�ʱ��
            % if tag == 1 || tag == 2             % ȷ����Ƶ\��ʱ��ͬ��ͷͼ��
                % fh_pat_rec = fh_pat_mode1;
                % th_pat_rec = th_pat_mode1;
                wav_S1_rec = wav_S1_mode1;    
                wav_S2_rec = wav_S2_mode1;    
            % elseif tag == 3 || tag == 4
            %     fh_pat_rec = fh_pat_2_mode1;
            %     th_pat_rec = th_pat_2_mode1;
            %     wav_S1_rec = wav_S1_2_mode1;    
            %     wav_S2_rec = wav_S2_2_mode1; 
            % end
            time_frame = time_frame_mode1;  % һ֡���ܳ���(��th��103)
            rx_Cap = rx(i:i+time_frame*oversamp_IF+oversamp_IF*2);  %��ȡ��2��bit 
            result_frame = receiver_Cap_mode1(rx_Cap, fh_pat_mode1, th_pat_mode1, wav_S1_S3_mode4, wav_S4_S2_mode4, wav_S1_rec, wav_S2_rec, time_frame, num_pulses_rec);  % ����2Mbps Aģʽ��Ӧ�Ľ��ģ����
            result(ceil((i+(time_frame-0.5)*oversamp_IF)/(time_frame*oversamp_IF)), 1:length_frame) = result_frame;  % ��¼������
            
        elseif mode_sel == 2
            num_pulses_rec = num_pulses_mode2;  % һ֡����������
            length_frame = length_frame_mode2;  % һ֡����Ч���ݳ���
            % time_counter = time_counter_mode2;  % �����һ֡�ĵȴ�ʱ��
            % if tag == 1 || tag == 2             % ȷ����Ƶ\��ʱ��ͬ��ͷͼ��
            %     fh_pat_rec = fh_pat_1_mode2;
            %     th_pat_rec = th_pat_1_mode2;
                wav_S1_rec = wav_S1_mode2;    
                wav_S2_rec = wav_S2_mode2;    
            % elseif tag == 3 || tag == 4
            %     fh_pat_rec = fh_pat_2_mode2;
            %     th_pat_rec = th_pat_2_mode2;
            %     wav_S1_rec = wav_S1_2_mode2;    
            %     wav_S2_rec = wav_S2_2_mode2; 
            % end 
            time_frame = time_frame_mode2;      % һ֡���ܳ���(����ʱ)
            rx_Cap = rx(i:i+time_frame*oversamp_IF+oversamp_IF*2);
            result_frame = receiver_Cap_mode2(rx_Cap, fh_pat_mode2, th_pat_mode2, wav_S1_S3_mode4, wav_S4_S2_mode4, wav_S1_rec, wav_S2_rec, time_frame, num_pulses_rec);  % ����2Mbps Bģʽ��Ӧ�Ľ��ģ����
            result(ceil((i+(time_frame-0.5)*oversamp_IF)/(time_frame*oversamp_IF)), 1:length_frame) = result_frame;  % ��¼��������ȡ����
            
        elseif mode_sel == 3
            num_pulses_rec = num_pulses_mode3;  % һ֡����������
            length_frame = length_frame_mode3;  % һ֡����Ч���ݳ���
            % time_counter = time_counter_mode3;  % �����һ֡�ĵȴ�ʱ��
            % if tag == 1 || tag == 2             % ȷ����Ƶ\��ʱ��ͬ��ͷͼ��
            %     fh_pat_rec = fh_pat_1_mode3;
            %     th_pat_rec = th_pat_1_mode3;
                wav_S1_rec = wav_S1_S3_mode3;    
                wav_S2_rec = wav_S4_S2_mode3;    
            % elseif tag == 3 || tag == 4
            %     fh_pat_rec = fh_pat_2_mode3;
            %     th_pat_rec = th_pat_2_mode3;
            %     wav_S1_rec = wav_S1_S3_2_mode3;    
            %     wav_S2_rec = wav_S4_S2_2_mode3; 
            % end 
            time_frame = time_frame_mode3;      % һ֡���ܳ���(����ʱ)
            rx_Cap = rx(i:i+time_frame*oversamp_IF+oversamp_IF*2);
            result_frame = receiver_Cap_mode3(rx_Cap, fh_pat_mode3, th_pat_mode3, wav_S1_S3_mode4, wav_S4_S2_mode4, wav_S1_rec, wav_S2_rec, time_frame, num_pulses_rec);   % ����500Kbpsģʽ��Ӧ�Ľ��ģ����
            result(ceil((i+(time_frame-0.5)*oversamp_IF)/((time_frame+10)*oversamp_IF)), 1:length_frame) = result_frame;  % ��¼������
            
        else
            num_pulses_rec = num_pulses_mode4;   % һ֡����������
            length_frame = length_frame_mode4;   % һ֡����Ч���ݳ���
            % time_counter = time_counter_mode4;   % �����һ֡�ĵȴ�ʱ��
            % if tag == 1 || tag == 2              % ȷ����Ƶ\��ʱ��ͬ��ͷͼ��
            %     fh_pat_rec = fh_pat_1_mode4;
            %     th_pat_rec = th_pat_1_mode4;
                % wav_S1_rec = wav_S1_S3_mode4;    
                % wav_S2_rec = wav_S4_S2_mode4;    
            % elseif tag == 3 || tag == 4
            %     fh_pat_rec = fh_pat_2_mode4;
            %     th_pat_rec = th_pat_2_mode4;
            %     wav_S1_rec = wav_S1_S3_2_mode4;    
            %     wav_S2_rec = wav_S4_S2_2_mode4; 
            % end 
            time_frame = time_frame_mode4;      % һ֡���ܳ���(����ʱ)
            rx_Cap = rx(i:i+time_frame*oversamp_IF+oversamp_IF*2);
            result_frame = receiver_Cap_mode4(rx_Cap, fh_pat_mode4, th_pat_mode4, wav_S1_S3_mode4, wav_S4_S2_mode4, time_frame, num_pulses_rec);   % ����250Kbpsģʽ��Ӧ�Ľ��ģ����
            result(ceil((i+(time_frame-0.5)*oversamp_IF)/(time_frame*oversamp_IF)), 1:length_frame) = result_frame;  % ��¼������
            
        end
        
        flag_Capture_C = 0;
        frame_counter = frame_counter + 1;  % ���񵽵�֡������1
        % flag_frame = 1; 
        % tag = 0;
       
        if frame_counter == frame_num+1  %4֡�����Ѿ�ȫ��������
            return 
        end
   
    end

end  % (for i = ...)

figure(3);
k=1:37;
plot(k, corr_value_mode1_plot1);
figure(4);
k=1:37;
plot(k, corr_value_mode1_plot2);