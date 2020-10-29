function [phi, I_sig_BB, Q_sig_BB, int_phi] = GMSK(bit_5, f_c, phi_last, g)
% input: 5比特数据              bit_5 
%        当前跳频载波            fc 
%        上一bit末相位          phi_last(即bit_5的第三个bit值)
%        g                    对应BT=0.3的g函数
% output：当前bit的末相位        phi
%         GMSK调制I路基带波形    I_sig_BB
%         GMSK调制Q路基带波形    Q_sig_BB
%         对应当前bit的相位路径   int_phi             

bit_rate = 16e6;  % 数据速率
T = 1/bit_rate;  % 符号时间
fs_IF = 1024e6;  % 中频信号采样频率
oversamp_IF = T*fs_IF;  % 过采样倍数
T_s_IF = 1/fs_IF;  % 采样时间
BT = 0.3;  % 高斯滤波器 BT带宽

c_bit = bit_5;

% 相位累加
intphi = zeros(1,T*fs_IF);

for i = 1:T*fs_IF
    if i == 1
        intphi(i) = phi_last + c_bit(5)*g(i)*T_s_IF + c_bit(4)*g(i + T*fs_IF)*T_s_IF + c_bit(3)*g(i + 2*T*fs_IF)*T_s_IF + c_bit(2)*g(i + 3*T*fs_IF)*T_s_IF + c_bit(1)*g(i + 4*T*fs_IF)*T_s_IF;
    else
        intphi(i) = intphi(i-1) + c_bit(5)*g(i)*T_s_IF + c_bit(4)*g(i + T*fs_IF)*T_s_IF + c_bit(3)*g(i + 2*T*fs_IF)*T_s_IF + c_bit(2)*g(i + 3*T*fs_IF)*T_s_IF + c_bit(1)*g(i + 4*T*fs_IF)*T_s_IF;
    end
end


t = -T/2:1/fs_IF:T/2-1/fs_IF;
t = t(1:end) + (T/oversamp_IF/2);  % 中心对称
int_phi = mod(pi/(2*T).*intphi, 2*pi);

phi = int_phi(T*fs_IF)/(pi/(2*T));

% I路信号
I_sig_BB = cos(int_phi);
I_sig_IF = cos(2*pi*f_c.*t);


% Q路信号
Q_sig_BB = sin(int_phi);
Q_sig_IF = sin(2*pi*f_c.*t);

