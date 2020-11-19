clear all; 
close all;

% rng(0);
Ne = 10;

bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间

f_IF = 240e6; %射频频率
fs_IF = 64e6;  % 射频、中频频信号采样速率
fs_BB = 128e6;  % 基带信号采样速率
num_bits_pulse = 5000; 
oversamp_BB = T * fs_BB;  % 基带信号过采样速率
oversamp_IF = T * fs_IF;  % 射频、中频信号过采样速率
T_s_BB = 1/fs_BB;  % 基带采样间隔

load('lib/g_1024.mat');  % GMSK调制 g函数 
g = g(1:16:end);

error_index = 1;
error_cnt = zeros(1,5);
Eb_N0_cnt = zeros(1,5);

for Eb_N0 = -4:10

error_rate = 0;

for ll = 1 : Ne

% generate code

I_single = randi(2,1,num_bits_pulse);
I_single = I_single - 1;

% LDPC

% I_ldpc = encoder(I_single', num_bits_pulse, num_bits_pulse*3);

% load('test_I.mat');
% I = 2*(I_ldpc') - 1;

% I_single = [0,0,1,0,0,0,1,1,1,1];
% I = [1,-1,1,-1,-1,-1,1,1,1,1];

% 预编码

% I_pre = precode(I_single);
% I = 2*I_pre - 1;
I = 2*I_single - 1;

% coding
bit_5 = zeros(1,5);


phi_last = 0;



for i = 1:num_bits_pulse
    if i == 1
        bit_5 = [-1,-1,I(i:i+2)];
    elseif i == 2
        bit_5 = [-1,I(i-1:i+2)];
    elseif i == num_bits_pulse-1
        bit_5 = [I(i-2:i+1),-1];
    elseif i == num_bits_pulse
        bit_5 = [I(i-2:i),-1,-1];
    else
        bit_5 = I(i-2:i+2);
    end

    [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, phi_last, g);
    signal_trans_BB((i-1)*oversamp_IF+1:(i)*oversamp_IF) = complex(I_sig, Q_sig);
    phi_all((i-1)*oversamp_IF+1:(i)*oversamp_IF) = phi_int;
end

% figure;
% plot(mod(phi_all,2*pi))

% 射频信号

t = linspace(0, num_bits_pulse*T*3, oversamp_IF*num_bits_pulse);
signal_trans_IF = signal_trans_BB .* exp(1i*2*pi*0*t);

fre = t./(num_bits_pulse*T*3)*fs_IF;
% figure;
% plot(fre, 20*log10(abs(fft(signal_trans_IF))))

% 加频偏
% signal_trans_IF = signal_trans_IF .* exp(1i*2*pi*0.001*bit_rate*t);


% 加噪声
% Eb_N0 = 5;
SNRdB = Eb_N0 - 10*log10(oversamp_IF);
signal_recv_IF_noise = awgn(signal_trans_IF, SNRdB, 'measured');
% signal_recv_IF_noise = signal_trans_IF;
signal_recv_noise_IF_FFT = abs(fft(signal_recv_IF_noise));

% figure;
% plot(fre, 20*log10(signal_recv_noise_IF_FFT))

% 带通滤波

% 低通滤波

LPF = [0.000419533866398058,-0.000888354496183992,-0.000853437623768592,-0.000631298799865539,-2.40800381611274e-05,0.000757726612508301,0.00124279157153801,0.00100718942187309,6.11349594777944e-06,-0.00127966699253509,-0.00204136110437968,-0.00162139817984884,-7.51640527054630e-06,0.00200063649880008,0.00314685992485242,0.00246771090895958,8.74772351090950e-06,-0.00298552899026191,-0.00465094660214941,-0.00361541320946284,-1.02589693362414e-05,0.00431736672636526,0.00668350449440606,0.00516716692440808,1.16677091975264e-05,-0.00612641623913192,-0.00945499326702595,-0.00729373821942378,-1.26285747392181e-05,0.00864688543125808,0.0133562430870217,0.0103253357103565,1.39560789926429e-05,-0.0123588569838066,-0.0192224558203244,-0.0149957418809320,-1.48710949738265e-05,0.0184487731613767,0.0292381503162025,0.0233590153736128,1.56686296873214e-05,-0.0308372382510145,-0.0514618173035326,-0.0440872983709100,-1.56557954715510e-05,0.0744478850711668,0.158618847348271,0.224900901783822,0.250015876706213,0.224900901783822,0.158618847348271,0.0744478850711668,-1.56557954715510e-05,-0.0440872983709100,-0.0514618173035326,-0.0308372382510145,1.56686296873214e-05,0.0233590153736128,0.0292381503162025,0.0184487731613767,-1.48710949738265e-05,-0.0149957418809320,-0.0192224558203244,-0.0123588569838066,1.39560789926429e-05,0.0103253357103565,0.0133562430870217,0.00864688543125808,-1.26285747392181e-05,-0.00729373821942378,-0.00945499326702595,-0.00612641623913192,1.16677091975264e-05,0.00516716692440808,0.00668350449440606,0.00431736672636526,-1.02589693362414e-05,-0.00361541320946284,-0.00465094660214941,-0.00298552899026191,8.74772351090950e-06,0.00246771090895958,0.00314685992485242,0.00200063649880008,-7.51640527054630e-06,-0.00162139817984884,-0.00204136110437968,-0.00127966699253509,6.11349594777944e-06,0.00100718942187309,0.00124279157153801,0.000757726612508301,-2.40800381611274e-05,-0.000631298799865539,-0.000853437623768592,-0.000888354496183992,0.000419533866398058];
% LPF = [-0.000802208597432932,0.000274012677990304,0.000642510473715215,0.000946999415062327,0.000891916050444873,0.000343793840903912,-0.000536240328039012,-0.00131477084764724,-0.00148902995381680,-0.000790675452817240,0.000581761681437109,0.00196476830994624,0.00252374305662866,0.00172239576761462,-0.000280356461407567,-0.00256991435306819,-0.00385674604167267,-0.00316803805115370,-0.000490591828793483,0.00302142690249472,0.00549129901414302,0.00527836650663518,0.00198141321136188,-0.00308987510505797,-0.00735036957404938,-0.00819657411948721,-0.00450206260885427,0.00246334247774722,0.00932162509793802,0.0121291137400687,0.00851827609598724,-0.000667699637068033,-0.0112575451563172,-0.0174798363625150,-0.0149003745646248,-0.00316325538761452,0.0129951262378319,0.0253189360751059,0.0258185755732741,0.0111788546136904,-0.0143738689839372,-0.0394190876857844,-0.0492672140398497,-0.0320968020127244,0.0152607586838193,0.0845460387353867,0.157862956932812,0.213627410283510,0.234433292058803,0.213627410283510,0.157862956932812,0.0845460387353867,0.0152607586838193,-0.0320968020127244,-0.0492672140398497,-0.0394190876857844,-0.0143738689839372,0.0111788546136904,0.0258185755732741,0.0253189360751059,0.0129951262378319,-0.00316325538761452,-0.0149003745646248,-0.0174798363625150,-0.0112575451563172,-0.000667699637068033,0.00851827609598724,0.0121291137400687,0.00932162509793802,0.00246334247774722,-0.00450206260885427,-0.00819657411948721,-0.00735036957404938,-0.00308987510505797,0.00198141321136188,0.00527836650663518,0.00549129901414302,0.00302142690249472,-0.000490591828793483,-0.00316803805115370,-0.00385674604167267,-0.00256991435306819,-0.000280356461407567,0.00172239576761462,0.00252374305662866,0.00196476830994624,0.000581761681437109,-0.000790675452817240,-0.00148902995381680,-0.00131477084764724,-0.000536240328039012,0.000343793840903912,0.000891916050444873,0.000946999415062327,0.000642510473715215,0.000274012677990304,-0.000802208597432932];
signal_recv_BB = filtfilt(LPF, 1, signal_recv_IF_noise);
% signal_recv_BB = signal_recv_IF_noise;
% figure;
% plot(fre, 20*log10(abs(fft(signal_recv_BB))))


% viterbi译码

% signal_recv_BB = signal_trans_IF;
signal_recv_BB = [signal_recv_BB, zeros(1,2*oversamp_IF)];
de_out_index = 1;

de_out = zeros(size(I));
soft_info = zeros(size(I));
viterbi_deep = 40;
% signal_trans_IF(1 : oversamp_IF*viterbi_deep)
path_record = cell(128, 2);
path_weight = cell(128, 2);
signal_recv_dif = signal_recv_BB(1:4) .* conj(signal_recv_BB(9:12));
for q = 1:32
    init_path = dec2bin(q-1, 5) - '0';
    init_path = [0,0,init_path];
    path_record{q,1} = init_path;
    bit_5_fst = init_path(end-6:end-2);
    bit_5_sec = init_path(end-5:end-1);
    bit_5_trd = init_path(end-4:end);
    [phi_last, I_sig_fst, Q_sig_fst, ~] = GMSK(2*bit_5_fst-1, 0, g);
    [phi_last, ~, ~, ~] = GMSK(2*bit_5_sec-1, phi_last, g);
    [~, I_sig_trd, Q_sig_trd, ~] = GMSK(2*bit_5_trd-1, phi_last, g);
    g_path_dif = complex(I_sig_fst, Q_sig_fst) .* conj(complex(I_sig_trd, Q_sig_trd));
    path_weight_m = path_weight{q, 1};
    path_weight_m = [path_weight_m, real(sum(signal_recv_dif .* conj(g_path_dif)))]; % 取最大值
    path_weight{q, 1} = path_weight_m;
end  

for n = 2:num_bits_pulse
    
    signal_recv_dif = signal_recv_BB(1+(n-1)*4:n*4) .* conj(signal_recv_BB(1+(n+1)*4:(n+2)*4));

    for j = 1 : 128
        if ~isempty(path_record{j,1})
            bit_record = path_record{j,1};
            bit_record_cur0 = [bit_record, 0];
            bit_record_cur1 = [bit_record, 1];
            
            bit_5_fst0 = bit_record_cur0(end-6:end-2);
            bit_5_sec0 = bit_record_cur0(end-5:end-1);
            bit_5_trd0 = bit_record_cur0(end-4:end);
            bit_5_fst1 = bit_record_cur1(end-6:end-2);
            bit_5_sec1 = bit_record_cur1(end-5:end-1);
            bit_5_trd1 = bit_record_cur1(end-4:end);
            
            [phi_last0, I_sig0_fst, Q_sig0_fst, ~] = GMSK(2*bit_5_fst0-1, 0, g);
            [phi_last0, ~, ~, ~] = GMSK(2*bit_5_sec0-1, phi_last0, g);
            [~, I_sig0_trd, Q_sig0_trd, ~] = GMSK(2*bit_5_trd0-1, phi_last0, g);

            [phi_last1, I_sig1_fst, Q_sig1_fst, ~] = GMSK(2*bit_5_fst1-1, 0, g);
            [phi_last1, ~, ~, ~] = GMSK(2*bit_5_sec1-1, phi_last1, g);
            [~, I_sig1_trd, Q_sig1_trd, ~] = GMSK(2*bit_5_trd1-1, phi_last1, g);
            
            state_index0 = bin2dec(num2str(bit_record_cur0(end-6:end)))+1;
            state_index1 = bin2dec(num2str(bit_record_cur1(end-6:end)))+1;
            
            g_path_dif0 = complex(I_sig0_fst, Q_sig0_fst) .* conj(complex(I_sig0_trd, Q_sig0_trd));
            g_path_dif1 = complex(I_sig1_fst, Q_sig1_fst) .* conj(complex(I_sig1_trd, Q_sig1_trd));
            
            path_weight0 = real(sum(signal_recv_dif .* conj(g_path_dif0)));
            path_weight1 = real(sum(signal_recv_dif .* conj(g_path_dif1)));

            path_weight0_new = [path_weight{j,1}, path_weight0];
            path_weight1_new = [path_weight{j,1}, path_weight1];

            if isempty(path_weight{state_index0,2}) || sum(path_weight0_new) > sum(path_weight{state_index0,2})
                path_weight{state_index0, 2} = path_weight0_new;
                path_record{state_index0, 2} = bit_record_cur0;                   
            end
            
            if isempty(path_weight{state_index1,2}) || sum(path_weight1_new) > sum(path_weight{state_index1,2})
                path_weight{state_index1, 2} = path_weight1_new;
                path_record{state_index1, 2} = bit_record_cur1;                                   
            end
        end
    end

    for t = 1:128
        path_record{t,1} = path_record{t,2};
        path_record{t,2} = [];
        path_weight{t,1} = path_weight{t,2};
        path_weight{t,2} = [];
    end

    % path_weight(:,1) = path_weight(:,2);
    % path_weight(:,2) = zeros(128,1);

    if n == num_bits_pulse
        path_weight_mat = cell2mat(path_weight);
        path_weight_sum = sum(path_weight_mat, 2);
        [~,index_max] = max(path_weight_sum);
        out = path_record{index_max,1};
        de_out(de_out_index : end) = out(3:end-4);

        for w = 3:size(out,2)-4
            z_m = 0;
            o_m = 0;
            for r = 1:128
                p_tmp = path_record{r,1};
                if p_tmp(w) == 0
                    if path_weight_sum(r,1) > z_m
                        z_m = path_weight_sum(r,1);
                    end
                else
                    if path_weight_sum(r,1) > o_m
                        o_m = path_weight_sum(r,1);
                    end                  
                end
            end
            soft_info(de_out_index) = z_m - o_m;
            de_out_index = de_out_index + 1;
        end
        
    else
        if size(path_record{1,1},2) == viterbi_deep
            path_weight_mat = cell2mat(path_weight);
            path_weight_sum = sum(path_weight_mat, 2);
            [~,index_max] = max(path_weight_sum);
            out = path_record{index_max,1};
            z_m = 0;
            o_m = 0;
            for r = 1:128
                p_tmp = path_record{r,1};
                if p_tmp(3) == 0
                    if path_weight_sum(r,1) > z_m
                        z_m = path_weight_sum(r,1);
                    end
                else
                    if path_weight_sum(r,1) > o_m
                        o_m = path_weight_sum(r,1);
                    end                  
                end
            end
            soft_info(de_out_index) = z_m - o_m;
            de_out(de_out_index) = out(3);
            de_out_index = de_out_index + 1;
            for t = 1:128
                tmp = path_record{t,1};
                w_tmp = path_weight{t,1};
                path_record{t,1} = tmp(2:end);
                path_weight{t,1} = w_tmp(2:end);
            end
        end
    end
end

% [Lde, ~] = Ldecoder2(soft_info', num_bits_pulse, num_bits_pulse*3);

% 解预编码

% de_out_pre = decode_pre(de_out);

error = I_single - de_out;

error(error~=0) = 1;
error(end-3:end) = 0;

error_rate = error_rate + sum(error);
end

error_rate/num_bits_pulse/Ne;

error_cnt(error_index) = error_rate/num_bits_pulse/Ne;
Eb_N0_cnt(error_index) = Eb_N0;
error_index = error_index + 1;

end

% BER plot

f = figure;
f.PaperUnits = 'centimeters';
f.PaperSize = [16, 12];
f.Units = 'centimeters';
f.Position = [0, 0, 16, 12];
semilogy(Eb_N0_cnt, error_cnt, '-ks', 'LineWidth', 2);
xlim([min(Eb_N0_cnt)-1, max(Eb_N0_cnt)+1]);
grid on;


