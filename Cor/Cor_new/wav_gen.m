function [wav_S1, wav_S2] = wav_gen(pn_lib_1, pn_lib_2, pn_lib_3, pn_lib_4, mode)

%  生成同步头序列对应的本地波形

load('lib/g_1024.mat');
g = g(1:16:end);  % g函数
bit_rate = 16e6;  % 数据速率
T = 1/bit_rate;  % Tb
f_s = 64e6;  % 采样频率
oversamp = T * f_s; % 过采样倍数
f_trans = 32;  % （无用）

% GMSK调制
% 2Mbps A 速率模式
if mode == 1
    num_pulses = 12;
    num_bits_pn = 24;
    pn = [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];  % 2Mbps A 模式的同步头扰码序列
    
    % 序列S1
    for idx = 1:num_pulses

        phi_last = 0;
        preamble = double(xor(pn_lib_1(idx,:), pn));  % 原始同步头序列与扰码序列异或得到最终的同步头序列
        preamble = 2*preamble-1;  %双极性
        for i = 1:num_bits_pn

            if (i==1)
                bit_5 = [preamble(i), preamble(i), preamble(i:i+2)];
            elseif (i==2)
                bit_5 = [preamble(i-1), preamble(i-1:i+2)];
            elseif (i==num_bits_pn-1)
                bit_5 = [preamble(i-2:i+1), 0];
            elseif (i==num_bits_pn)
                bit_5 = [preamble(i-2:i), 0, 0];
            else
                bit_5 = preamble(i-2:i+2);
            end

            [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, phi_last, g);
            wav_S1_64(idx, (i-1)*oversamp+1:(i)*oversamp) = complex(I_sig, Q_sig);  %创建复数

        end
        wav_S1(idx, :) = wav_S1_64(idx, 4:4:end);

    end
    
    % 序列S2
    for idx = 1:num_pulses

        phi_last = 0;
        preamble = double(xor(pn_lib_2(idx,:), pn));
        preamble = 2*preamble-1;
        
        for i = 1:num_bits_pn

            if (i==1)
                bit_5 = [preamble(i), preamble(i), preamble(i:i+2)];
            elseif (i==2)
                bit_5 = [preamble(i-1), preamble(i-1:i+2)];
            elseif (i==num_bits_pn-1)
                bit_5 = [preamble(i-2:i+1), 0];
            elseif (i==num_bits_pn)
                bit_5 = [preamble(i-2:i), 0, 0];
            else
                bit_5 = preamble(i-2:i+2);
            end

            [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, phi_last, g);
            wav_S2_64(idx, (i-1)*oversamp+1:(i)*oversamp) = complex(I_sig, Q_sig);

        end
        wav_S2(idx, :) = wav_S2_64(idx, 4:4:end);

    end
    
% 2Mbps B速率模式
elseif mode == 2
    
    num_pulses = 12;
    num_bits_pn = 24;
    pn = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0];  % 2Mbps B 模式的同步头扰码序列
    % 序列S1
    for idx = 1:num_pulses

        phi_last = 0;
        preamble = double(xor(pn_lib_1(idx,:), pn));
        preamble = 2*preamble-1;

        for i = 1:num_bits_pn

            if (i==1)
                bit_5 = [preamble(i), preamble(i), preamble(i:i+2)];
            elseif (i==2)
                bit_5 = [preamble(i-1), preamble(i-1:i+2)];
            elseif (i==num_bits_pn-1)
                bit_5 = [preamble(i-2:i+1), 0];
            elseif (i==num_bits_pn)
                bit_5 = [preamble(i-2:i), 0, 0];
            else
                bit_5 = preamble(i-2:i+2);
            end

            [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, phi_last, g);
            wav_S1_64(idx, (i-1)*oversamp+1:(i)*oversamp) = complex(I_sig, Q_sig);

        end
        wav_S1(idx, :) = wav_S1_64(idx, 4:4:end);

    end
    
    % 序列S2
    for idx = 1:num_pulses

        phi_last = 0;
        preamble = double(xor(pn_lib_2(idx,:), pn));
        preamble = 2*preamble-1;

        for i = 1:num_bits_pn

            if (i==1)
                bit_5 = [preamble(i), preamble(i), preamble(i:i+2)];
            elseif (i==2)
                bit_5 = [preamble(i-1), preamble(i-1:i+2)];
            elseif (i==num_bits_pn-1)
                bit_5 = [preamble(i-2:i+1), 0];
            elseif (i==num_bits_pn)
                bit_5 = [preamble(i-2:i), 0, 0];
            else
                bit_5 = preamble(i-2:i+2);
            end

            [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, phi_last, g);
            wav_S2_64(idx, (i-1)*oversamp+1:(i)*oversamp) = complex(I_sig, Q_sig);

        end
        wav_S2(idx, :) = wav_S2_64(idx, 4:4:end);

    end
    
% 500K速率模式    
elseif mode == 3
    num_pulses = 48;
    num_bits_pn = 45;
    pn = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];  % 500Kbps 模式的同步头S1、S2扰码序列
    pn_2 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];  % 500Kbps 模式的同步头S3、S4扰码序列
    % 序列S1、S3
    for idx = 1:num_pulses

        phi_last = 0;
        preamble = [double(xor(pn_lib_1(idx,:), pn)), double(xor(pn_lib_3(idx,:), pn_2))];
        preamble = 2*preamble-1;

        for i = 1:num_bits_pn

            if (i==1)
                bit_5 = [preamble(i), preamble(i), preamble(i:i+2)];
            elseif (i==2)
                bit_5 = [preamble(i-1), preamble(i-1:i+2)];
            elseif (i==num_bits_pn-1)
                bit_5 = [preamble(i-2:i+1), 0];
            elseif (i==num_bits_pn)
                bit_5 = [preamble(i-2:i), 0, 0];
            else
                bit_5 = preamble(i-2:i+2);
            end

            [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, phi_last, g);
            wav_S1_64(idx, (i-1)*oversamp+1:(i)*oversamp) = complex(I_sig, Q_sig);

        end
        wav_S1(idx, :) = wav_S1_64(idx, 4:4:end);

    end
    
    % 序列S4、S2
    for idx = 1:num_pulses

        phi_last = 0;
        preamble = [double(xor(pn_lib_4(idx,:), pn_2)), double(xor(pn_lib_2(idx,:), pn))];
        preamble = 2*preamble-1;

        for i = 1:num_bits_pn

            if (i==1)
                bit_5 = [preamble(i), preamble(i), preamble(i:i+2)];
            elseif (i==2)
                bit_5 = [preamble(i-1), preamble(i-1:i+2)];
            elseif (i==num_bits_pn-1)
                bit_5 = [preamble(i-2:i+1), 0];
            elseif (i==num_bits_pn)
                bit_5 = [preamble(i-2:i), 0, 0];
            else
                bit_5 = preamble(i-2:i+2);
            end

            [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, phi_last, g);
            wav_S2_64(idx, (i-1)*oversamp+1:(i)*oversamp) = complex(I_sig, Q_sig);

        end
        wav_S2(idx, :) = wav_S2_64(idx, 4:4:end);

    end
    
% 250K速率模式    
else
    num_pulses = 48;
    num_bits_pn = 45;
    pn = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    pn_2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    % 序列S1、S3
    for idx = 1:num_pulses

        phi_last = 0;
        preamble = [double(xor(pn_lib_1(idx,:), pn)), double(xor(pn_lib_3(idx,:), pn_2))];
        preamble = 2*preamble-1;

        for i = 1:num_bits_pn

            if (i==1)
                bit_5 = [preamble(i), preamble(i), preamble(i:i+2)];
            elseif (i==2)
                bit_5 = [preamble(i-1), preamble(i-1:i+2)];
            elseif (i==num_bits_pn-1)
                bit_5 = [preamble(i-2:i+1), 0];
            elseif (i==num_bits_pn)
                bit_5 = [preamble(i-2:i), 0, 0];
            else
                bit_5 = preamble(i-2:i+2);
            end

            [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, phi_last, g);
            wav_S1_64(idx, (i-1)*oversamp+1:(i)*oversamp) = complex(I_sig, Q_sig);

        end
        wav_S1_temp(idx, :) = wav_S1_64(idx, 4:4:end); %和数据速率相同

    end
    
    % 序列S4、S2
    for idx = 1:num_pulses

        phi_last = 0;
        preamble = [double(xor(pn_lib_4(idx,:), pn_2)), double(xor(pn_lib_2(idx,:), pn))];
        preamble = 2*preamble-1;

        for i = 1:num_bits_pn

            if (i==1)
                bit_5 = [preamble(i), preamble(i), preamble(i:i+2)];
            elseif (i==2)
                bit_5 = [preamble(i-1), preamble(i-1:i+2)];
            elseif (i==num_bits_pn-1)
                bit_5 = [preamble(i-2:i+1), 0];
            elseif (i==num_bits_pn)
                bit_5 = [preamble(i-2:i), 0, 0];
            else
                bit_5 = preamble(i-2:i+2);
            end

            [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, phi_last, g);
            wav_S2_64(idx, (i-1)*oversamp+1:(i)*oversamp) = complex(I_sig, Q_sig);

        end
        wav_S2_temp(idx, :) = wav_S2_64(idx, 4:4:end); %和数据速率相同

    end
    
    wav_S1 = repmat(wav_S1_temp, [2,1]);  %重复
    wav_S2 = repmat(wav_S2_temp, [2,1]);

end
    
    
    
    
