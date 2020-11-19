clear all; 
close all;

% rng(0);
Ne = 100;
KK = 3; %编码倍数

bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间

f_IF = 240e6; %射频频率
fs_IF = 64e6;  % 射频、中频频信号采样速率
fs_BB = 128e6;  % 基带信号采样速率
num_data_frame = 1024;
num_data_pulse = 256; 
num_bit_pulse = 256 + 24*2; 
num_pn = 24;
num_pulse = 12; % 脉冲个数
oversamp_BB = T * fs_BB;  % 基带信号过采样速率
oversamp_IF = T * fs_IF;  % 射频、中频信号过采样速率
T_s_BB = 1/fs_BB;  % 基带采样间隔

load('lib/g_1024.mat');  % GMSK调制 g函数 
g = g(1:16:end);

error_index = 1;
error_cnt = zeros(1,5);
Eb_N0_cnt = zeros(1,5);

for Eb_N0 = -1.6 : 0.2 : -1

error_rate = 0;
bulk_rate = 0;

for ll = 1 : Ne

% 随机产生单极性数据
I_single = randi(2,1,num_data_frame);
I_single = I_single - 1;

% LDPC编码
I_ldpc = encoder(I_single', num_data_frame, num_data_frame*KK);
I_ldpc = I_ldpc';

% 随机生成同步头
S1 = randi(2,num_pulse,num_pn);
S1 = S1 - 1;

S2 = randi(2,num_pulse,num_pn);
S2 = S2 - 1;

% GMSK调制
signal_trans_BB_total = zeros(num_pulse, num_bit_pulse*oversamp_IF);
for z = 1:num_pulse
    % 预编码
    I_pre = precode(I_ldpc(1+(z-1)*num_data_pulse:z*num_data_pulse));
    % 相位不连续
    [GMSK_S1, ~] = GMSK_mode(S1(z,:)*2-1, num_pn, oversamp_IF, 0, g);
    [GMSK_data, ~] = GMSK_mode(I_pre*2-1, num_data_pulse, oversamp_IF, 0, g);
    [GMSK_S2, ~] = GMSK_mode(S2(z,:)*2-1, num_pn, oversamp_IF, 0, g);
    signal_trans_BB_total(z,:) = [GMSK_S1, GMSK_data, GMSK_S2];
end

% figure;
% plot(mod(phi_all,2*pi))

% 射频信号
t = linspace(0, num_bit_pulse*T, size(signal_trans_BB_total,2));
signal_trans_IF = signal_trans_BB_total .* repmat(exp(1i*2*pi*0*t), [num_pulse, 1]);

fre = t./(num_pulse*T)*fs_IF;
% figure;
% plot(fre, 20*log10(abs(fft(signal_trans_IF))))

% 加频偏
signal_trans_IF = signal_trans_IF .* repmat(exp(1i*2*pi*0.001*bit_rate*t), [num_pulse, 1]);

% 加相偏
signal_trans_IF = signal_trans_IF .* repmat(exp(1i*2*pi*0.01), [num_pulse, 1]);

% 加噪声
SNRdB = Eb_N0 - 10*log10(oversamp_IF);
signal_recv_IF_noise = awgn(signal_trans_IF, SNRdB, 'measured');
% signal_recv_IF_noise = signal_trans_IF;

% figure;
% signal_recv_noise_IF_FFT = abs(fft(signal_recv_IF_noise));
% plot(fre, 20*log10(signal_recv_noise_IF_FFT))

% 低通滤波

%6 ~ 9
% LPF = [-8.31034897815694e-05,5.49834554170762e-06,8.96876467242553e-05,0.000182775376084856,0.000204473627532599,8.74203445655418e-05,-0.000161331758553307,-0.000428209226566028,-0.000528376381192817,-0.000311046901275352,0.000215550062278594,0.000823214635593835,0.00113560596377343,0.000830741285692785,-0.000120282576005487,-0.00133057976730626,-0.00211224231674170,-0.00182795136500055,-0.000319397993091591,0.00184485354626735,0.00351571398626448,0.00352113505210602,0.00139291110843283,-0.00215170782116458,-0.00533881767105549,-0.00615205941116020,-0.00349940105073380,0.00189627645440703,0.00749063477936225,0.0100053932360894,0.00720823826116176,-0.000518486075301420,-0.00979344523556798,-0.0155243543466729,-0.0134820054404147,-0.00295639205309983,0.0120019610308354,0.0237581424851113,0.0245442258244239,0.0107828833143032,-0.0138425440097149,-0.0384047054131546,-0.0483823249881099,-0.0317475518716441,0.0150648890280452,0.0841128410987495,0.157546084932118,0.213572158128231,0.234506593320022,0.213572158128231,0.157546084932118,0.0841128410987495,0.0150648890280452,-0.0317475518716441,-0.0483823249881099,-0.0384047054131546,-0.0138425440097149,0.0107828833143032,0.0245442258244239,0.0237581424851113,0.0120019610308354,-0.00295639205309983,-0.0134820054404147,-0.0155243543466729,-0.00979344523556798,-0.000518486075301420,0.00720823826116176,0.0100053932360894,0.00749063477936225,0.00189627645440703,-0.00349940105073380,-0.00615205941116020,-0.00533881767105549,-0.00215170782116458,0.00139291110843283,0.00352113505210602,0.00351571398626448,0.00184485354626735,-0.000319397993091591,-0.00182795136500055,-0.00211224231674170,-0.00133057976730626,-0.000120282576005487,0.000830741285692785,0.00113560596377343,0.000823214635593835,0.000215550062278594,-0.000311046901275352,-0.000528376381192817,-0.000428209226566028,-0.000161331758553307,8.74203445655418e-05,0.000204473627532599,0.000182775376084856,8.96876467242553e-05,5.49834554170762e-06,-8.31034897815694e-05];

%6.7 ~8.7
% LPF = [-0.00105986524153879,-0.000516080434296824,-0.000144578495591309,0.000469023161791935,0.000992009860433083,0.00104535498099612,0.000448916377280878,-0.000599609325005874,-0.00154902636196073,-0.00176315084387256,-0.000913448937108964,0.000729100247319563,0.00232347381230053,0.00286014482715793,0.00176327599546244,-0.000650621952960946,-0.00317608897861725,-0.00428318593331421,-0.00299961327427968,0.000373766980623597,0.00419642545207672,0.00623270936864497,0.00490434828518391,0.000374957424087243,-0.00522180554163347,-0.00871066821855678,-0.00760114889544637,-0.00171135329239916,0.00630577776550046,0.0120338840728796,0.0116110633999896,0.00413276079122796,-0.00727121706809576,-0.0165307758001576,-0.0176649491481550,-0.00831351424643959,0.00815857462057874,0.0234442483711252,0.0280452285076553,0.0164517415020282,-0.00879883000843695,-0.0365325948593853,-0.0506811931626110,-0.0371332594708111,0.00924665036632640,0.0807673987528880,0.158368715861701,0.218184196315662,0.240635037339314,0.218184196315662,0.158368715861701,0.0807673987528880,0.00924665036632640,-0.0371332594708111,-0.0506811931626110,-0.0365325948593853,-0.00879883000843695,0.0164517415020282,0.0280452285076553,0.0234442483711252,0.00815857462057874,-0.00831351424643959,-0.0176649491481550,-0.0165307758001576,-0.00727121706809576,0.00413276079122796,0.0116110633999896,0.0120338840728796,0.00630577776550046,-0.00171135329239916,-0.00760114889544637,-0.00871066821855678,-0.00522180554163347,0.000374957424087243,0.00490434828518391,0.00623270936864497,0.00419642545207672,0.000373766980623597,-0.00299961327427968,-0.00428318593331421,-0.00317608897861725,-0.000650621952960946,0.00176327599546244,0.00286014482715793,0.00232347381230053,0.000729100247319563,-0.000913448937108964,-0.00176315084387256,-0.00154902636196073,-0.000599609325005874,0.000448916377280878,0.00104535498099612,0.000992009860433083,0.000469023161791935,-0.000144578495591309,-0.000516080434296824,-0.00105986524153879];
% 7 ~ 9
LPF = [0.000419533866398058,-0.000888354496183992,-0.000853437623768592,-0.000631298799865539,-2.40800381611274e-05,0.000757726612508301,0.00124279157153801,0.00100718942187309,6.11349594777944e-06,-0.00127966699253509,-0.00204136110437968,-0.00162139817984884,-7.51640527054630e-06,0.00200063649880008,0.00314685992485242,0.00246771090895958,8.74772351090950e-06,-0.00298552899026191,-0.00465094660214941,-0.00361541320946284,-1.02589693362414e-05,0.00431736672636526,0.00668350449440606,0.00516716692440808,1.16677091975264e-05,-0.00612641623913192,-0.00945499326702595,-0.00729373821942378,-1.26285747392181e-05,0.00864688543125808,0.0133562430870217,0.0103253357103565,1.39560789926429e-05,-0.0123588569838066,-0.0192224558203244,-0.0149957418809320,-1.48710949738265e-05,0.0184487731613767,0.0292381503162025,0.0233590153736128,1.56686296873214e-05,-0.0308372382510145,-0.0514618173035326,-0.0440872983709100,-1.56557954715510e-05,0.0744478850711668,0.158618847348271,0.224900901783822,0.250015876706213,0.224900901783822,0.158618847348271,0.0744478850711668,-1.56557954715510e-05,-0.0440872983709100,-0.0514618173035326,-0.0308372382510145,1.56686296873214e-05,0.0233590153736128,0.0292381503162025,0.0184487731613767,-1.48710949738265e-05,-0.0149957418809320,-0.0192224558203244,-0.0123588569838066,1.39560789926429e-05,0.0103253357103565,0.0133562430870217,0.00864688543125808,-1.26285747392181e-05,-0.00729373821942378,-0.00945499326702595,-0.00612641623913192,1.16677091975264e-05,0.00516716692440808,0.00668350449440606,0.00431736672636526,-1.02589693362414e-05,-0.00361541320946284,-0.00465094660214941,-0.00298552899026191,8.74772351090950e-06,0.00246771090895958,0.00314685992485242,0.00200063649880008,-7.51640527054630e-06,-0.00162139817984884,-0.00204136110437968,-0.00127966699253509,6.11349594777944e-06,0.00100718942187309,0.00124279157153801,0.000757726612508301,-2.40800381611274e-05,-0.000631298799865539,-0.000853437623768592,-0.000888354496183992,0.000419533866398058];
% LPF = [-0.000802208597432932,0.000274012677990304,0.000642510473715215,0.000946999415062327,0.000891916050444873,0.000343793840903912,-0.000536240328039012,-0.00131477084764724,-0.00148902995381680,-0.000790675452817240,0.000581761681437109,0.00196476830994624,0.00252374305662866,0.00172239576761462,-0.000280356461407567,-0.00256991435306819,-0.00385674604167267,-0.00316803805115370,-0.000490591828793483,0.00302142690249472,0.00549129901414302,0.00527836650663518,0.00198141321136188,-0.00308987510505797,-0.00735036957404938,-0.00819657411948721,-0.00450206260885427,0.00246334247774722,0.00932162509793802,0.0121291137400687,0.00851827609598724,-0.000667699637068033,-0.0112575451563172,-0.0174798363625150,-0.0149003745646248,-0.00316325538761452,0.0129951262378319,0.0253189360751059,0.0258185755732741,0.0111788546136904,-0.0143738689839372,-0.0394190876857844,-0.0492672140398497,-0.0320968020127244,0.0152607586838193,0.0845460387353867,0.157862956932812,0.213627410283510,0.234433292058803,0.213627410283510,0.157862956932812,0.0845460387353867,0.0152607586838193,-0.0320968020127244,-0.0492672140398497,-0.0394190876857844,-0.0143738689839372,0.0111788546136904,0.0258185755732741,0.0253189360751059,0.0129951262378319,-0.00316325538761452,-0.0149003745646248,-0.0174798363625150,-0.0112575451563172,-0.000667699637068033,0.00851827609598724,0.0121291137400687,0.00932162509793802,0.00246334247774722,-0.00450206260885427,-0.00819657411948721,-0.00735036957404938,-0.00308987510505797,0.00198141321136188,0.00527836650663518,0.00549129901414302,0.00302142690249472,-0.000490591828793483,-0.00316803805115370,-0.00385674604167267,-0.00256991435306819,-0.000280356461407567,0.00172239576761462,0.00252374305662866,0.00196476830994624,0.000581761681437109,-0.000790675452817240,-0.00148902995381680,-0.00131477084764724,-0.000536240328039012,0.000343793840903912,0.000891916050444873,0.000946999415062327,0.000642510473715215,0.000274012677990304,-0.000802208597432932];
% signal_recv_BB = filtfilt(LPF, 1, signal_recv_IF_noise);
for j = 1:num_pulse
    signal_recv_tmp = conv(signal_recv_IF_noise(j,:), LPF);
    signal_recv_BB(j,:) = signal_recv_tmp(48+1:48+length(signal_recv_IF_noise));
end

% signal_recv_BB = signal_recv_IF_noise;
% figure;
% plot(fre, 20*log10(abs(fft(signal_recv_BB))))

% 提取接收信号的同步头与数据
S1_store = zeros(num_pulse, num_pn*oversamp_IF);
S2_store = zeros(num_pulse, num_pn*oversamp_IF);
data_store = zeros(num_pulse, num_data_pulse*oversamp_IF);
% t_S1 = zeros(size(S1_store));
% t_S2 = zeros(size(S2_store));
% t_data = zeros(size(data_store));

S1_store = signal_recv_BB(:, 1:num_pn*oversamp_IF);
data_store = signal_recv_BB(:, 1+num_pn*oversamp_IF:(num_pn+num_data_pulse)*oversamp_IF);
S2_store = signal_recv_BB(:, (num_pn+num_data_pulse)*oversamp_IF+1:(num_pn+num_data_pulse+num_pn)*oversamp_IF);

% 产生本地同步头
for e = 1:num_pulse
    [S1_wav(e,:), ~] = GMSK_mode(2*S1(e,:)-1, num_pn, oversamp_IF, 0, g);
    [S2_wav(e,:), ~] = GMSK_mode(2*S2(e,:)-1, num_pn, oversamp_IF, 0, g);    
    % 过低通滤波
    % S1_wav_tmp = conv(S1_wav(e,:), LPF);
    % S1_wav_lpf(e,:) = S1_wav_tmp(48+1:48+length(S1_wav(e,:)));
    % S2_wav_tmp = conv(S2_wav(e,:), LPF);
    % S2_wav_lpf(e,:) = S2_wav_tmp(48+1:48+length(S2_wav(e,:)));
end

% 相关估计频偏
zk_S1 = zeros(num_pulse, num_pn*oversamp_IF);
zk_S2 = zeros(num_pulse, num_pn*oversamp_IF);
for j = 1:num_pulse                  
    zk_S1(j,:) = S1_store(j,:) .* conj(S1_wav(j,:));
    zk_S2(j,:) = S2_store(j,:) .* conj(S2_wav(j,:)); 
    zk_corr_S12 = conj(zk_S1(j,:)) .* zk_S2(j,:);
    deltat_hat_S12(j) = (atan2(-imag(sum(zk_corr_S12)), real(sum(zk_corr_S12)))); 
    delta_f_hat_S12(j) = deltat_hat_S12(j)/(2*pi*280*T);
end

dif_w_f = mean(delta_f_hat_S12);

% dif_S1_w = zeros(num_pulse, num_pn*oversamp_IF);
% dif_S2_w = zeros(num_pulse, num_pn*oversamp_IF);
% for j = 1:num_pulse
%     for f = 1:num_pn*oversamp_IF/2
%         dif_S1_w(j,f) = conj(zk_S1(j, f)) .* zk_S1(j, f+num_pn*oversamp_IF/2);
%         dif_S2_w(j,f) = conj(zk_S2(j, f)) .* zk_S2(j, f+num_pn*oversamp_IF/2);
%     end
% end

% dif_w = sum(sum(dif_S1_w,2)) + sum(sum(dif_S2_w,2));
% dif_w_ang = atan2(-imag(dif_w), real(dif_w));
% dif_w_f = dif_w_ang/(2*pi*11*T);

% % S2的初始相位
% DF_hat_S12_phi = deltat_hat_S12 - dif_w_ang/11*282;

% 同步头消除频偏
counter_f_S1 = repmat(dif_w_f,[num_pulse, num_pn*oversamp_IF]) * 2 * pi .* repmat(t(1 : num_pn*oversamp_IF), [num_pulse, 1]);
counter_f_S2 = repmat(dif_w_f,[num_pulse, num_pn*oversamp_IF]) * 2 * pi .* repmat(t((num_bit_pulse-num_pn)*oversamp_IF+1:(num_bit_pulse)*oversamp_IF), [num_pulse, 1]);
D_S1_half = S1_store .* complex(cos(counter_f_S1), sin(counter_f_S1));
D_S2_half = S2_store .* complex(cos(counter_f_S2), sin(counter_f_S2));

% 利用本地波形计算相偏
delta_theta_S1 = zeros(num_pulse, 1);
delta_theta_S2 = zeros(num_pulse, 1);
for j = 1:num_pulse
    delta_theta_cplx_mat_S1 = D_S1_half(j,:) .* conj(S1_wav(j,:));
    delta_theta_cplx_mat_S2 = D_S2_half(j,:) .* conj(S2_wav(j,:)); 
    delta_theta_cplx_S1 = sum(delta_theta_cplx_mat_S1) / (num_pn)*oversamp_IF;
    delta_theta_cplx_S2 = sum(delta_theta_cplx_mat_S2) / (num_pn)*oversamp_IF;
    delta_theta_S1(j) = atan2(-imag(delta_theta_cplx_S1), real(delta_theta_cplx_S1));
    delta_theta_S2(j) = atan2(-imag(delta_theta_cplx_S2), real(delta_theta_cplx_S2));
end

% delta_theta_S1
% delta_theta_S2
dif_phi = (delta_theta_S1+delta_theta_S2)/2;

% 频偏相偏的修正

data_store = data_store .* exp(1i*repmat(dif_w_f,[num_pulse,num_data_pulse*oversamp_IF]) * 2 * pi .* repmat(t(1+num_pn*oversamp_IF:(num_pn+num_data_pulse)*oversamp_IF), [num_pulse, 1])) .* exp(1i*dif_phi);

soft_info_pre = [];
for v = 1:num_pulse
    % 获取解调数据的初相位
    % [~, phi_init] = GMSK_mode(S1(v,:)*2-1, num_pn, oversamp_IF, 0, g);
    % viterbi译码
    [~, soft_info] = GMSK_viterbi(data_store(v,:), num_data_pulse, oversamp_IF, 0, g);
    soft_info_pre(1+(v-1)*num_data_pulse:v*num_data_pulse) = soft_info;
    % [~, soft_info] = GMSK_viterbi(data_store(v,:), num_data_pulse, 4, g);
    % soft_info_pre(1+(v-1)*num_data_pulse:v*num_data_pulse) = soft_info;
end
[Lde, ~] = Ldecoder2(soft_info_pre', num_data_frame, num_data_frame*KK);

error = I_single - Lde';

error(error~=0) = 1;
% error(end-3:end) = 0;

if sum(error)
    bulk_rate = bulk_rate + 1;
end

error_rate = error_rate + sum(error);
end

erro_record = error_rate/num_data_frame/Ne
bulk_record = bulk_rate

error_cnt(error_index) = error_rate/num_data_frame/Ne;
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


