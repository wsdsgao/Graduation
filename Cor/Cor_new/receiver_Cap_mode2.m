function result = receiver_Cap_mode2(rx, fh_pat, th_pat, wav_S1_S3_mode4, wav_S4_S2_mode4, wav_S1, wav_S2, time_frame, num_pulses)

% 2Mbps Bģʽ �Ľ��ģ��

% ��������
oversamp_IF = 64;   % ��Ƶ����Ƶ�źŹ���������
oversamp_BB = 4;    % �����źŹ���������
num_bits_pulse = 304;  % 2Mbps A\500Kbps\250Kbps һ��������ĳ���
                       % 2Mbps B ÿ������ȥ��ǰ�������ݲ��ֺ�ĳ���
bit_rate = 16e6;  % ��������
T = 1/bit_rate;   % ����ʱ��

num_bits_pn = 24; % ͬ��ͷS1\S2����

t2 = 0:1/oversamp_IF:num_bits_pulse-1/oversamp_IF;
t2 = t2(1:end) + (T/oversamp_IF/2);  % ���ĶԳ�

nn = oversamp_IF / 8;

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
% load('lib/filter/LPF.mat');  % ��ͨ�˲��� 61�� ͨ��: 4.8MHz ����Ƶ��128MHz
% load('lib/filter/LPF_2.mat');  % ��ͨ�˲��� 507�� ͨ��: 5MHz  ����Ƶ��1024MHz
% S_lpf = 30;
% S_lpf2 = 127;
% S_bpf = 253;
            
Wav_str_Cap_mat = zeros(nn, time_frame*oversamp_IF);
% ���8����������ȡһ·���Σ���8·
for offset = 1:nn
    Wav_str_Cap_mat(offset,:) = rx(1+(offset-1)*8:time_frame*oversamp_IF+(offset-1)*8);
end

Wav_str_Cap_F_temp = zeros(num_pulses, num_bits_pulse*oversamp_IF, nn);
Corr_S1 = zeros(1, nn);
Corr_S2 = zeros(1, nn); 
for offset = 1:nn

    temp_rx = Wav_str_Cap_mat(offset,:);

    for pulse_idx = 1:num_pulses

        f_idx = fh_pat(pulse_idx); % ��ǰ�����Ӧ��Ƶ��

        % ���㵱ǰ������һ֡�е���ʼλ��
        if pulse_idx == 1
            pos_pulse = 0;
        else
            pos_pulse = sum(th_pat(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103);
        end

        % ���㵱ǰ����ĳ���
        num_bits_pulse_rx = th_pat(pulse_idx)  + num_bits_pulse;

        % ��Ƶ -> ��Ƶ
        if (f_idx >= 1) && (f_idx <= 5)    
            wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % ͨ��1
        elseif (f_idx >= 6) && (f_idx <= 10)
            wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % ͨ��2
        elseif (f_idx >= 11) && (f_idx <= 14)
            wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % ͨ��3
        elseif (f_idx >= 15) && (f_idx <= 18)
            wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % ͨ��4 
        elseif (f_idx >= 19) && (f_idx <= 21)
            wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % ͨ��5 
        end

        Wav_str_Cap_F_temp(pulse_idx,:,offset) = wav_temp1((th_pat(pulse_idx)/2)*oversamp_IF+1:(th_pat(pulse_idx)/2+num_bits_pulse)*oversamp_IF);
    end

    % ������Ƶͼ����Ӧ��Ƶ�㽫��Ӧͨ���Ĳ����˲�������Ƶ
    % wav_temp(������1024) -> rx_pulse_mat(������64)
    rx_pulse_mat = downConv(Wav_str_Cap_F_temp(:,:,offset), num_pulses, fh_pat);

    Corr_mode2(offset) = corr(rx_pulse_mat, wav_S1_S3_mode4, wav_S4_S2_mode4, 2);
    % % Ԥȡǰ��24bitλ�õĲ��� �ȴ���������
    % D_S1 = rx_pulse_mat(:,1:24*oversamp_BB);
    % D_S1_one = D_S1(:,4:oversamp_BB:end);
    % D_S2 = rx_pulse_mat(:,1120+1:1120+24*oversamp_BB);
    % D_S2_one = D_S2(:,4:oversamp_BB:end);


    % % ͬ��ͷ����ϸУ׼                
    % % ���ñ��ص�PN����������� 
    % rx_corr_S1 = zeros(1,num_pulses);
    % rx_corr_S2 = zeros(1,num_pulses);
    % for j = 1:num_pulses

    %     rx_wav_S1 = D_S1_one(j,1:22) .* conj(wav_S1(j,1:22));
    %     rx_wav_S2 = D_S2_one(j,1:22) .* conj(wav_S2(j,1:22));   

    %     rx_corr_S1(j) = abs(sum(rx_wav_S1));
    %     rx_corr_S2(j) = abs(sum(rx_wav_S2));

    % end

    % Corr_S1(offset) = sum(rx_corr_S1);
    % Corr_S2(offset) = sum(rx_corr_S2);

end

pos_Corr = find(Corr_mode2 == max(Corr_mode2));

% pos_Corr_S1 = find(Corr_S1 == max(Corr_S1));
% pos_Corr_S2 = find(Corr_S2 == max(Corr_S2));

% if (tag == 1) || (tag == 3)
%     pos_Corr = pos_Corr_S1;
% elseif (tag == 2) || (tag == 4)
%     pos_Corr = pos_Corr_S2;
% end


% ȡ9·���� ���1��������
for offset = 1:9
    Wav_str_Cap_mat(offset,:) = rx((pos_Corr-1)*8+(offset-1)+1-4:time_frame*oversamp_IF+(pos_Corr-1)*8+(offset-1)-4);
end            

Wav_str_Cap_F_temp2 = zeros(num_pulses, num_bits_pulse*oversamp_IF, 9);
Corr_S1_Cap_F = zeros(1, 9);
Corr_S2_Cap_F = zeros(1, 9);
for offset = 1:9

    temp_rx = Wav_str_Cap_mat(offset,:);

    for pulse_idx = 1:num_pulses

        f_idx = fh_pat(pulse_idx);  % ��ǰ�����Ӧ��Ƶ��

        % ���㵱ǰ������һ֡�е���ʼλ��
        if pulse_idx == 1
            pos_pulse = 0;
        else
            pos_pulse = sum(th_pat(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103);
        end
        
        % ���㵱ǰ���峤��
        num_bits_pulse_rx = th_pat(pulse_idx) + num_bits_pulse;

        % ��Ƶ -> ��Ƶ
        if (f_idx >= 1) && (f_idx <= 5)    
            wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % ͨ��1
        elseif (f_idx >= 6) && (f_idx <= 10)
            wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % ͨ��2
        elseif (f_idx >= 11) && (f_idx <= 14)
            wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % ͨ��3
        elseif (f_idx >= 15) && (f_idx <= 18)
            wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % ͨ��4 
        elseif (f_idx >= 19) && (f_idx <= 21)
            wav_temp1 = downConv_IF(temp_rx(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % ͨ��5 
        end

        Wav_str_Cap_F_temp2(pulse_idx,:,offset) = wav_temp1((th_pat(pulse_idx)/2)*oversamp_IF+1:(th_pat(pulse_idx)/2+num_bits_pulse)*oversamp_IF);
    end

    % ������Ƶͼ����Ӧ��Ƶ�㽫��Ӧͨ���Ĳ����˲�������Ƶ
    % wav_temp(������1024) -> rx_pulse_mat(������64)
    rx_pulse_mat = downConv(Wav_str_Cap_F_temp2(:,:,offset), num_pulses, fh_pat); 

    Corr_Cap_F(offset) = corr(rx_pulse_mat, wav_S1_S3_mode4, wav_S4_S2_mode4, 2);

    % % Ԥȡǰ��24bitλ�õĲ��� �ȴ���������
    % D_S1 = rx_pulse_mat(:,1:24*oversamp_BB);
    % D_S1_one = D_S1(:,8:oversamp_BB:end);
    % D_S2 = rx_pulse_mat(:,2240+1:2240+24*oversamp_BB);
    % D_S2_one = D_S2(:,8:oversamp_BB:end);


    % % ͬ��ͷ����ϸУ׼                
    % % ���ñ��ص�PN����Ƶƫ����         
    % % ��8·Ƶ�ʶ���֮��Ĳ����뱾�ز��ε����ֵ��ÿ·���8�������㣩
    % for j = 1:num_pulses

    %     rx_wav_S1 = D_S1_one(j,1:22) .* conj(wav_S1(j,1:22));
    %     rx_wav_S2 = D_S2_one(j,1:22) .* conj(wav_S2(j,1:22));   

    %     rx_corr_S1(j) = abs(sum(rx_wav_S1));
    %     rx_corr_S2(j) = abs(sum(rx_wav_S2));

    % end

    % Corr_S1_Cap_F(offset) = sum(rx_corr_S1);
    % Corr_S2_Cap_F(offset) = sum(rx_corr_S2);

end

pos_Corr_F = find(Corr_Cap_F == max(Corr_Cap_F));

offset_i = pos_Corr;
offset_j = pos_Corr_F;

% pos_Corr_S1_F = find(Corr_S1_Cap_F == max(Corr_S1_Cap_F));
% pos_Corr_S2_F = find(Corr_S2_Cap_F == max(Corr_S2_Cap_F));

% if (tag == 1) || (tag == 3)
%     offset_i = pos_Corr_S1
%     offset_j = pos_Corr_S1_F
%     choose_flag = 1;

% elseif (tag == 2) || (tag == 4)
%     offset_i = pos_Corr_S2
%     offset_j = pos_Corr_S2_F
%     choose_flag = 2;

% end

temp_rx_FNL = rx((offset_i-1)*8+(offset_j-1)+1-4:time_frame*oversamp_IF+(offset_i-1)*8+(offset_j-1)-4);
Wav_temp_FNL = zeros(num_pulses, num_bits_pulse * oversamp_IF);
Wav_store_FNL = zeros(num_pulses, 512 * oversamp_IF);
for pulse_idx = 1:num_pulses

    f_idx = fh_pat(pulse_idx); % ��ǰ�����Ӧ��Ƶ��

    % ���㵱ǰ������һ֡�е���ʼλ��
    if pulse_idx == 1
        pos_pulse = 0;
    else
        pos_pulse = sum(th_pat(1:pulse_idx-1)) + (pulse_idx-1)*(num_bits_pulse+103);
    end

    % ���㵱ǰ����ĳ���
    num_bits_pulse_rx = th_pat(pulse_idx) + num_bits_pulse;

    % ��Ƶ -> ��Ƶ
    if (f_idx >= 1) && (f_idx <= 5)    
        wav_temp1 = downConv_IF(temp_rx_FNL(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN1, BPF_CHAN12_2, 240e6 - f_trans(3), 0);  % ͨ��1
    elseif (f_idx >= 6) && (f_idx <= 10)
        wav_temp1 = downConv_IF(temp_rx_FNL(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN2, BPF_CHAN12_2, 240e6 - f_trans(8), 0);  % ͨ��2
    elseif (f_idx >= 11) && (f_idx <= 14)
        wav_temp1 = downConv_IF(temp_rx_FNL(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN3, BPF_CHAN34_2, 240e6 - (f_trans(12)+13.33e6/2), 0);   % ͨ��3
    elseif (f_idx >= 15) && (f_idx <= 18)
        wav_temp1 = downConv_IF(temp_rx_FNL(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN4, BPF_CHAN34_2, 240e6 - (f_trans(16)+13.33e6/2), 0);   % ͨ��4 
    elseif (f_idx >= 19) && (f_idx <= 21)
        wav_temp1 = downConv_IF(temp_rx_FNL(pos_pulse*oversamp_IF+1:pos_pulse*oversamp_IF+num_bits_pulse_rx*oversamp_IF), BPF_CHAN5, BPF_CHAN5_2, f_trans(20) - 240e6, 1);  % ͨ��5 
    end

    Wav_temp_FNL(pulse_idx,:) = wav_temp1((th_pat(pulse_idx)/2)*oversamp_IF+1:(th_pat(pulse_idx)/2+num_bits_pulse)*oversamp_IF); %�м����ݲ���
    Wav_store_FNL(pulse_idx, 1:(th_pat(pulse_idx))*oversamp_IF) = [wav_temp1(1:(th_pat(pulse_idx)/2)*oversamp_IF), wav_temp1((th_pat(pulse_idx)/2+num_bits_pulse)*oversamp_IF+1:end)]; %�����ݲ���

end


% ������Ƶͼ����Ӧ��Ƶ�㽫��Ӧͨ���Ĳ����˲�������Ƶ
% wav_temp(������1024) -> rx_pulse_mat(������64)
rx_pulse_mat_FNL = downConv_2MLL(Wav_temp_FNL, num_pulses, fh_pat, th_pat); %����˲�����������  

% Ԥȡǰ��24bitλ�õĲ��� �ȴ���������
D_S1 = rx_pulse_mat_FNL(:,1:num_bits_pn*oversamp_BB);  %���ǰ�����������ģ�飬�����ҲҪ��
D_S1_one = D_S1(:,4:oversamp_BB:end);
D_S2 = rx_pulse_mat_FNL(:,1120+1:1120+num_bits_pn*oversamp_BB);
D_S2_one = D_S2(:,4:oversamp_BB:end);

zk_S1 = zeros(num_pulses, (num_bits_pn-2));
zk_S2 = zeros(num_pulses, (num_bits_pn-2));
for j = 1:num_pulses                  
    zk_S1(j,:) = D_S1_one(j,1:num_bits_pn-2) .* conj(wav_S1(j,1:num_bits_pn-2));
    zk_S2(j,:) = D_S2_one(j,1:num_bits_pn-2) .* conj(wav_S2(j,1:num_bits_pn-2));
end

% ��������
for j = 1:num_pulses 
    
    zk_corr_S12(j) = sum(conj(zk_S1(j,:)) .* zk_S2(j,:));
    % zk_corr_S12 = conj(zk_S1(j,:)) .* zk_S2(j,:);
    % deltat_hat_S12(j) = (atan2(-imag(sum(zk_corr_S12)), real(sum(zk_corr_S12)))); 
    % delta_f_hat_S12(j) = deltat_hat_S12(j)/(2*pi*280*T);

end
DF_hat_S12 = atan2(-imag(mean(zk_corr_S12)), real(mean(zk_corr_S12)))/(2*pi*280*T);

% ͬ��ͷ����Ƶƫ
counter_f_S1 = repmat(DF_hat_S12,[num_pulses,1]) * 2 * pi * t2(16:16:24*oversamp_IF) * T;
counter_f_S2 = repmat(DF_hat_S12,[num_pulses,1]) * 2 * pi * t2(280*oversamp_IF+16:16:280*oversamp_IF+24*oversamp_IF) * T;
D_S1_half = D_S1 .* complex(cos(counter_f_S1), sin(counter_f_S1));  % ��ǰ��ͬ��ͷλ�õĲ�������Ƶƫ
D_S2_half = D_S2 .* complex(cos(counter_f_S2), sin(counter_f_S2));

% �ٴ����ñ��ص�PN���μ�����ƫ
delta_theta = zeros(num_pulses, 1);
% delta_theta_S2 = zeros(num_pulses, 1);
for j = 1:num_pulses
    delta_theta_cplx_mat_S1 = D_S1_half(j,4:oversamp_BB:end) .* conj(wav_S1(j,:));
    delta_theta_cplx_mat_S2 = D_S2_half(j,4:oversamp_BB:end) .* conj(wav_S2(j,:));
    delta_theta_cplx = (sum(delta_theta_cplx_mat_S1(1:22)) + sum(delta_theta_cplx_mat_S2(1:22)))/(22*2); %��22����
    delta_theta(j) = atan2(-imag(delta_theta_cplx), real(delta_theta_cplx));
    % delta_theta_cplx_S1 = sum(delta_theta_cplx_mat_S1(1:22)) / 22;
    % delta_theta_cplx_S2 = sum(delta_theta_cplx_mat_S2(1:22)) / 22;
    % delta_theta_S1(j) = atan2(-imag(delta_theta_cplx_S1), real(delta_theta_cplx_S1));
    % delta_theta_S2(j) = atan2(-imag(delta_theta_cplx_S2), real(delta_theta_cplx_S2));
end

delta_f = DF_hat_S12;
% if (choose_flag == 1)
%     delta_theta = delta_theta_S1;
%     choose_flag = 0;

% elseif (choose_flag == 2)
%     delta_theta = delta_theta_S2;
%     choose_flag = 0;
% end

rx_pulse_mat_mid = rx_pulse_mat_FNL(:, 1:280*oversamp_BB); %ǰͬ��ͷ+�м����ݲ���
rx_pulse_mat_pre = downConv_th(Wav_store_FNL(:,:), num_pulses, fh_pat, th_pat, 1);  % function downConv_th ֻ����2M@Bģʽ��ǰ��ʱ���ݲ��֣�
rx_pulse_mat_aft = [rx_pulse_mat_FNL(:,280*oversamp_BB+1:end), downConv_th(Wav_store_FNL(:,:), num_pulses, fh_pat, th_pat, 2)]; %��βͬ��ͷ+����ʱ���ݲ��֣�



% �Ե�ǰ֡�Ĳ�������֮ǰ�������õ�Ƶƫ����ƫ����ֵ������У����׼��������ģ��
counter_f_mid = repmat(delta_f, [num_pulses, 1]) * 2 * pi * t2(16:oversamp_IF/oversamp_BB:280*oversamp_IF) * T; %��16��ԭ����4
rx_pulse_mat_mid = rx_pulse_mat_mid .* complex(cos(counter_f_mid), sin(counter_f_mid)) .* repmat(complex(cos(delta_theta), sin(delta_theta)), [1, 280 * oversamp_BB]);            
for j = 1:num_pulses
    th = th_pat(j)/2;
    t_S1 = -th:1/oversamp_IF:0-1/oversamp_IF;
    t_S1 = t_S1(1:end) + (T/oversamp_IF/2);  
    t_S2 = 280:1/oversamp_IF:304+th-1/oversamp_IF;
    t_S2 = t_S2(1:end) + (T/oversamp_IF/2);  
    counter_f_pre = delta_f * 2 * pi * t_S1(16:oversamp_IF/oversamp_BB:end) * T;
    counter_f_aft = delta_f * 2 * pi * t_S2(16:oversamp_IF/oversamp_BB:end) * T;
    rx_pulse_mat_pre(j,1:th*oversamp_BB) = rx_pulse_mat_pre(j,1:th*oversamp_BB) .* complex(cos(counter_f_pre), sin(counter_f_pre)) .* repmat(complex(cos(delta_theta(j)), sin(delta_theta(j))), [1, th * oversamp_BB]); 
    rx_pulse_mat_aft(j,1:(th+24)*oversamp_BB) = rx_pulse_mat_aft(j,1:(th+24)*oversamp_BB) .* complex(cos(counter_f_aft), sin(counter_f_aft)) .* repmat(complex(cos(delta_theta(j)), sin(delta_theta(j))), [1, (th+24) * oversamp_BB]); 
end

%ֻȡ�����ݣ����ͬ��ͷ��4�����ʣ������ڽ��
rx_pulse_mat_FNL2 = cell(num_pulses, 1); %������cell��
for j = 1:num_pulses
    th = th_pat(j)/2;
    rx_pulse_mat_FNL2{j, 1} = zeros(1, (th*2+256)*oversamp_BB);
    rx_pulse_mat_FNL2{j, 1} = [rx_pulse_mat_pre(j, 1:th*oversamp_BB), rx_pulse_mat_mid(j, 24*oversamp_BB+1:280*oversamp_BB), rx_pulse_mat_aft(j, 24*oversamp_BB+1:(th+24)*oversamp_BB)];
end

 % GMSK ���
 % �����ݲ���
 out_temp_pre = zeros(num_pulses, 255);
 out_temp_aft = zeros(num_pulses, 255);
 for pulse_idx = 1:num_pulses
    th = th_pat(pulse_idx)/2;
    
    % ǰ�����ݲ���
    if mod(th+3,2)==1
        iter = th+3+1;
        flag_D = 1;
    else
        iter = th+3;
        flag_D = 0;
    end

    out_temp_pre(pulse_idx,1:(th+2)) = GMSK_demod(rx_pulse_mat_pre(pulse_idx, 1:(th+3)*oversamp_BB), c0_f, c1_f, oversamp_BB, flag_D, iter);                

    % �������ݲ���
    if mod(24+th,2)==1
        iter = 24+th+1;
        flag_D = 1;
    else
        iter = 24+th;
        flag_D = 0;
    end

    out_temp_aft(pulse_idx,1:(th+24-1)) = GMSK_demod(rx_pulse_mat_aft(pulse_idx, 1:(th+24)*oversamp_BB), c0_f, c1_f, oversamp_BB, flag_D, iter);                

 end           

 % �м䲿��
 out_temp_mid = zeros(num_pulses, 279);
 for pulse_idx = 1:num_pulses
    out_temp_mid(pulse_idx,:) = GMSK_demod(rx_pulse_mat_mid(pulse_idx, 1:(24+256)*oversamp_BB), c0_f, c1_f, oversamp_BB, 0, 280);
 end

 % ���������������һ֡
 result_last_bit = 0;
 for pulse_idx = 1:num_pulses
    th = th_pat(pulse_idx);
    result(result_last_bit+1:result_last_bit+th+256) = [out_temp_pre(pulse_idx, 3:th/2+2), out_temp_mid(pulse_idx, 24:end), out_temp_aft(pulse_idx, 24:th/2+24-1)];
    result_last_bit = result_last_bit + 256 + th;
 end

