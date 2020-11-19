clear all; 
close all;

% rng(0);
Ne = 1;
KK = 6; %编码倍数

bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间

f_IF = 240e6; %射频频率
fs_IF = 64e6;  % 射频、中频频信号采样速率
fs_BB = 128e6;  % 基带信号采样速率
num_bits_pulse = 1024; 
oversamp_BB = T * fs_BB;  % 基带信号过采样速率
oversamp_IF = T * fs_IF;  % 射频、中频信号过采样速率
T_s_BB = 1/fs_BB;  % 基带采样间隔

load('lib/g_1024.mat');  % GMSK调制 g函数 
g = g(1:16:end);

error_index = 1;
error_cnt = zeros(1,5);
Eb_N0_cnt = zeros(1,5);

for Eb_N0 = -1.4

error_rate = 0;

for ll = 1 : Ne

% generate code

I_single = randi(2,1,num_bits_pulse);
I_single = I_single - 1;

% load('test_I.mat');


% I_single = [0,0,1,0,0,0,1,1,1,1];
% I = [1,-1,1,-1,-1,-1,1,1,1,1];

% 预编码

I_pre = precode(I_single);
% I = 2*I_pre - 1;


% LDPC
I_ldpc = encoder(I_pre', num_bits_pulse, num_bits_pulse*KK/2);
I_ldpc6 = [I_ldpc; I_ldpc];
I = 2*(I_ldpc6') - 1;



% I = 2*I_single - 1;

% coding
bit_5 = zeros(1,5);


phi_last = 0;



for i = 1:num_bits_pulse*KK
    if i == 1
        bit_5 = [-1,-1,I(i:i+2)];
    elseif i == 2
        bit_5 = [-1,I(i-1:i+2)];
    elseif i == num_bits_pulse*KK-1
        bit_5 = [I(i-2:i+1),-1];
    elseif i == num_bits_pulse*KK
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

t = linspace(0, num_bits_pulse*T*KK, oversamp_IF*num_bits_pulse*KK);
signal_trans_IF = signal_trans_BB .* exp(1i*2*pi*0*t);

fre = t./(num_bits_pulse*KK*T)*fs_IF;
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

%6 ~ 9
% LPF = [-8.31034897815694e-05,5.49834554170762e-06,8.96876467242553e-05,0.000182775376084856,0.000204473627532599,8.74203445655418e-05,-0.000161331758553307,-0.000428209226566028,-0.000528376381192817,-0.000311046901275352,0.000215550062278594,0.000823214635593835,0.00113560596377343,0.000830741285692785,-0.000120282576005487,-0.00133057976730626,-0.00211224231674170,-0.00182795136500055,-0.000319397993091591,0.00184485354626735,0.00351571398626448,0.00352113505210602,0.00139291110843283,-0.00215170782116458,-0.00533881767105549,-0.00615205941116020,-0.00349940105073380,0.00189627645440703,0.00749063477936225,0.0100053932360894,0.00720823826116176,-0.000518486075301420,-0.00979344523556798,-0.0155243543466729,-0.0134820054404147,-0.00295639205309983,0.0120019610308354,0.0237581424851113,0.0245442258244239,0.0107828833143032,-0.0138425440097149,-0.0384047054131546,-0.0483823249881099,-0.0317475518716441,0.0150648890280452,0.0841128410987495,0.157546084932118,0.213572158128231,0.234506593320022,0.213572158128231,0.157546084932118,0.0841128410987495,0.0150648890280452,-0.0317475518716441,-0.0483823249881099,-0.0384047054131546,-0.0138425440097149,0.0107828833143032,0.0245442258244239,0.0237581424851113,0.0120019610308354,-0.00295639205309983,-0.0134820054404147,-0.0155243543466729,-0.00979344523556798,-0.000518486075301420,0.00720823826116176,0.0100053932360894,0.00749063477936225,0.00189627645440703,-0.00349940105073380,-0.00615205941116020,-0.00533881767105549,-0.00215170782116458,0.00139291110843283,0.00352113505210602,0.00351571398626448,0.00184485354626735,-0.000319397993091591,-0.00182795136500055,-0.00211224231674170,-0.00133057976730626,-0.000120282576005487,0.000830741285692785,0.00113560596377343,0.000823214635593835,0.000215550062278594,-0.000311046901275352,-0.000528376381192817,-0.000428209226566028,-0.000161331758553307,8.74203445655418e-05,0.000204473627532599,0.000182775376084856,8.96876467242553e-05,5.49834554170762e-06,-8.31034897815694e-05];

%6.7 ~8.7
% LPF = [-0.00105986524153879,-0.000516080434296824,-0.000144578495591309,0.000469023161791935,0.000992009860433083,0.00104535498099612,0.000448916377280878,-0.000599609325005874,-0.00154902636196073,-0.00176315084387256,-0.000913448937108964,0.000729100247319563,0.00232347381230053,0.00286014482715793,0.00176327599546244,-0.000650621952960946,-0.00317608897861725,-0.00428318593331421,-0.00299961327427968,0.000373766980623597,0.00419642545207672,0.00623270936864497,0.00490434828518391,0.000374957424087243,-0.00522180554163347,-0.00871066821855678,-0.00760114889544637,-0.00171135329239916,0.00630577776550046,0.0120338840728796,0.0116110633999896,0.00413276079122796,-0.00727121706809576,-0.0165307758001576,-0.0176649491481550,-0.00831351424643959,0.00815857462057874,0.0234442483711252,0.0280452285076553,0.0164517415020282,-0.00879883000843695,-0.0365325948593853,-0.0506811931626110,-0.0371332594708111,0.00924665036632640,0.0807673987528880,0.158368715861701,0.218184196315662,0.240635037339314,0.218184196315662,0.158368715861701,0.0807673987528880,0.00924665036632640,-0.0371332594708111,-0.0506811931626110,-0.0365325948593853,-0.00879883000843695,0.0164517415020282,0.0280452285076553,0.0234442483711252,0.00815857462057874,-0.00831351424643959,-0.0176649491481550,-0.0165307758001576,-0.00727121706809576,0.00413276079122796,0.0116110633999896,0.0120338840728796,0.00630577776550046,-0.00171135329239916,-0.00760114889544637,-0.00871066821855678,-0.00522180554163347,0.000374957424087243,0.00490434828518391,0.00623270936864497,0.00419642545207672,0.000373766980623597,-0.00299961327427968,-0.00428318593331421,-0.00317608897861725,-0.000650621952960946,0.00176327599546244,0.00286014482715793,0.00232347381230053,0.000729100247319563,-0.000913448937108964,-0.00176315084387256,-0.00154902636196073,-0.000599609325005874,0.000448916377280878,0.00104535498099612,0.000992009860433083,0.000469023161791935,-0.000144578495591309,-0.000516080434296824,-0.00105986524153879];
% 7 ~ 9
LPF = [0.000419533866398058,-0.000888354496183992,-0.000853437623768592,-0.000631298799865539,-2.40800381611274e-05,0.000757726612508301,0.00124279157153801,0.00100718942187309,6.11349594777944e-06,-0.00127966699253509,-0.00204136110437968,-0.00162139817984884,-7.51640527054630e-06,0.00200063649880008,0.00314685992485242,0.00246771090895958,8.74772351090950e-06,-0.00298552899026191,-0.00465094660214941,-0.00361541320946284,-1.02589693362414e-05,0.00431736672636526,0.00668350449440606,0.00516716692440808,1.16677091975264e-05,-0.00612641623913192,-0.00945499326702595,-0.00729373821942378,-1.26285747392181e-05,0.00864688543125808,0.0133562430870217,0.0103253357103565,1.39560789926429e-05,-0.0123588569838066,-0.0192224558203244,-0.0149957418809320,-1.48710949738265e-05,0.0184487731613767,0.0292381503162025,0.0233590153736128,1.56686296873214e-05,-0.0308372382510145,-0.0514618173035326,-0.0440872983709100,-1.56557954715510e-05,0.0744478850711668,0.158618847348271,0.224900901783822,0.250015876706213,0.224900901783822,0.158618847348271,0.0744478850711668,-1.56557954715510e-05,-0.0440872983709100,-0.0514618173035326,-0.0308372382510145,1.56686296873214e-05,0.0233590153736128,0.0292381503162025,0.0184487731613767,-1.48710949738265e-05,-0.0149957418809320,-0.0192224558203244,-0.0123588569838066,1.39560789926429e-05,0.0103253357103565,0.0133562430870217,0.00864688543125808,-1.26285747392181e-05,-0.00729373821942378,-0.00945499326702595,-0.00612641623913192,1.16677091975264e-05,0.00516716692440808,0.00668350449440606,0.00431736672636526,-1.02589693362414e-05,-0.00361541320946284,-0.00465094660214941,-0.00298552899026191,8.74772351090950e-06,0.00246771090895958,0.00314685992485242,0.00200063649880008,-7.51640527054630e-06,-0.00162139817984884,-0.00204136110437968,-0.00127966699253509,6.11349594777944e-06,0.00100718942187309,0.00124279157153801,0.000757726612508301,-2.40800381611274e-05,-0.000631298799865539,-0.000853437623768592,-0.000888354496183992,0.000419533866398058];
% LPF = [-0.000802208597432932,0.000274012677990304,0.000642510473715215,0.000946999415062327,0.000891916050444873,0.000343793840903912,-0.000536240328039012,-0.00131477084764724,-0.00148902995381680,-0.000790675452817240,0.000581761681437109,0.00196476830994624,0.00252374305662866,0.00172239576761462,-0.000280356461407567,-0.00256991435306819,-0.00385674604167267,-0.00316803805115370,-0.000490591828793483,0.00302142690249472,0.00549129901414302,0.00527836650663518,0.00198141321136188,-0.00308987510505797,-0.00735036957404938,-0.00819657411948721,-0.00450206260885427,0.00246334247774722,0.00932162509793802,0.0121291137400687,0.00851827609598724,-0.000667699637068033,-0.0112575451563172,-0.0174798363625150,-0.0149003745646248,-0.00316325538761452,0.0129951262378319,0.0253189360751059,0.0258185755732741,0.0111788546136904,-0.0143738689839372,-0.0394190876857844,-0.0492672140398497,-0.0320968020127244,0.0152607586838193,0.0845460387353867,0.157862956932812,0.213627410283510,0.234433292058803,0.213627410283510,0.157862956932812,0.0845460387353867,0.0152607586838193,-0.0320968020127244,-0.0492672140398497,-0.0394190876857844,-0.0143738689839372,0.0111788546136904,0.0258185755732741,0.0253189360751059,0.0129951262378319,-0.00316325538761452,-0.0149003745646248,-0.0174798363625150,-0.0112575451563172,-0.000667699637068033,0.00851827609598724,0.0121291137400687,0.00932162509793802,0.00246334247774722,-0.00450206260885427,-0.00819657411948721,-0.00735036957404938,-0.00308987510505797,0.00198141321136188,0.00527836650663518,0.00549129901414302,0.00302142690249472,-0.000490591828793483,-0.00316803805115370,-0.00385674604167267,-0.00256991435306819,-0.000280356461407567,0.00172239576761462,0.00252374305662866,0.00196476830994624,0.000581761681437109,-0.000790675452817240,-0.00148902995381680,-0.00131477084764724,-0.000536240328039012,0.000343793840903912,0.000891916050444873,0.000946999415062327,0.000642510473715215,0.000274012677990304,-0.000802208597432932];
signal_recv_BB = filtfilt(LPF, 1, signal_recv_IF_noise);
% signal_recv_BB = signal_recv_IF_noise;
% figure;
% plot(fre, 20*log10(abs(fft(signal_recv_BB))))


% viterbi译码

% signal_recv_BB = signal_trans_IF;
% signal_recv_BB = [signal_recv_BB, zeros(1,2*oversamp_IF)];
de_out_index = 1;

de_out = zeros(size(I));
soft_info = zeros(size(I));
viterbi_deep = 40;
% signal_trans_IF(1 : oversamp_IF*viterbi_deep)
path_record = cell(32, 2);
path_weight = cell(32, 2);
path_phi = zeros(32, 2);

signal_recv = signal_recv_BB(1:4);
for q = 1:8
    init_path = dec2bin(q-1, 3) - '0';
    init_path = [0,0,init_path];
    path_record{q,1} = init_path;
    [phi_last, I_sig, Q_sig, ~] = GMSK(2*init_path-1, 0, g);
    path_phi(q, 1) = phi_last;
    path_weight_m = path_weight{q, 1};
    path_weight_m = [path_weight_m, real(sum(signal_recv .* complex(I_sig, -Q_sig)))]; % 取最大值
    path_weight{q, 1} = path_weight_m;
end  

for n = 2:num_bits_pulse*KK
    
    signal_recv_dif = signal_recv_BB(1+(n-1)*4:n*4);

    for j = 1 : 32
        if ~isempty(path_record{j,1})
            bit_record = path_record{j,1};
            bit_record_cur0 = [bit_record, 0];
            bit_record_cur1 = [bit_record, 1];
            
            bit_5_0 = bit_record_cur0(end-4:end);
            bit_5_1 = bit_record_cur1(end-4:end);
  
            [phi_last0, I_sig0, Q_sig0, ~] = GMSK(2*bit_5_0-1, path_phi(j,1), g);
            [phi_last1, I_sig1, Q_sig1, ~] = GMSK(2*bit_5_1-1, path_phi(j,1), g);
 
            state_index0 = bin2dec(num2str(bit_record_cur0(end-4:end)))+1;
            state_index1 = bin2dec(num2str(bit_record_cur1(end-4:end)))+1;
            
            path_weight0 = real(sum(signal_recv_dif .* complex(I_sig0, -Q_sig0)));
            path_weight1 = real(sum(signal_recv_dif .* complex(I_sig1, -Q_sig1)));

            path_weight0_new = [path_weight{j,1}, path_weight0];
            path_weight1_new = [path_weight{j,1}, path_weight1];

            if n == num_bits_pulse-1 || n == num_bits_pulse
                if isempty(path_weight{state_index0,2}) || sum(path_weight0_new) > sum(path_weight{state_index0,2})
                    path_weight{state_index0, 2} = path_weight0_new;
                    path_record{state_index0, 2} = bit_record_cur0;    
                    path_phi(state_index0, 2) = phi_last0;         
                end
            else
                if isempty(path_weight{state_index0,2}) || sum(path_weight0_new) > sum(path_weight{state_index0,2})
                    path_weight{state_index0, 2} = path_weight0_new;
                    path_record{state_index0, 2} = bit_record_cur0;    
                    path_phi(state_index0, 2) = phi_last0;         
                end
                
                if isempty(path_weight{state_index1,2}) || sum(path_weight1_new) > sum(path_weight{state_index1,2})
                    path_weight{state_index1, 2} = path_weight1_new;
                    path_record{state_index1, 2} = bit_record_cur1; 
                    path_phi(state_index1, 2) = phi_last1;                                   
                end
            end
        end
    end


    len = size(path_weight{1,2},2);

    for t = 1:32
        path_record{t,1} = path_record{t,2};
        path_record{t,2} = [];
        path_weight{t,1} = path_weight{t,2};
        if isempty(path_weight{t,1})
            path_weight{t,1} = zeros(1,len);
        end
        path_weight{t,2} = [];
    end
    

    path_phi(:,1) = path_phi(:,2);
    path_phi(:,2) = path_phi(32,1);

    if n == num_bits_pulse*KK
        path_weight_mat = cell2mat(path_weight);
        path_weight_sum = sum(path_weight_mat, 2);
        [~,index_max] = max(path_weight_sum);
        out = path_record{index_max,1};
        de_out(de_out_index : end) = out(3:end-2);

        for w = 3:size(out,2)-2
            z_m = 0;
            o_m = 0;
            for r = 1:32
                if ~isempty(path_record{r,1})
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

            end
            soft_info(de_out_index) = z_m - o_m;
            de_out_index = de_out_index + 1;
            tmp_mat = path_weight_mat(:,2:end);
            path_weight_mat = [];
            path_weight_mat = tmp_mat;
            path_weight_sum = sum(path_weight_mat, 2);
        end
    else
        if size(path_record{1,1},2) == viterbi_deep
            path_weight_mat = cell2mat(path_weight);
            path_weight_sum = sum(path_weight_mat, 2);
            [~,index_max] = max(path_weight_sum);
            out = path_record{index_max,1};
            z_m = 0;
            o_m = 0;
            for r = 1:32
                if ~isempty(path_record{r,1})
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
            end
            soft_info(de_out_index) = z_m - o_m;
            de_out(de_out_index) = out(3);
            de_out_index = de_out_index + 1;
            for t = 1:32
                tmp = path_record{t,1};
                w_tmp = path_weight{t,1};
                path_record{t,1} = tmp(2:end);
                path_weight{t,1} = w_tmp(2:end);
            end
        end
    end
end




% 译LDPC
[Lde, ~] = Ldecoder2(soft_info', num_bits_pulse, num_bits_pulse*KK);

% 解预编码

de_out_pre = decode_pre(Lde');

error = I_single - de_out_pre;

error(error~=0) = 1;
% error(end-3:end) = 0;

error_rate = error_rate + sum(error);
end

error_rate/num_bits_pulse/Ne

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


