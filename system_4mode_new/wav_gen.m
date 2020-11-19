function [wav_S1, wav_S2] = wav_gen(pn_lib_1, pn_lib_2, pn_lib_3, pn_lib_4, mode)

%  ����ͬ��ͷ���ж�Ӧ�ı��ز���

load('lib/g_1024.mat');
g = g(1:16:end);  % g����
bit_rate = 16e6;  % ��������
T = 1/bit_rate;  % Tb
f_s = 64e6;  % ����Ƶ��
oversamp = T * f_s; % ����������
f_trans = 32;  % �����ã�

% GMSK����
% 2Mbps A ����ģʽ
if mode == 1
    num_pulses = 12;
    num_bits_pn = 24;
    pn = [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];  % 2Mbps A ģʽ��ͬ��ͷ��������
    
    % ����S1
    for idx = 1:num_pulses

        phi_last = 0;
        preamble = double(xor(pn_lib_1(idx,:), pn));  % ԭʼͬ��ͷ�����������������õ����յ�ͬ��ͷ����
        preamble = 2*preamble-1;  %˫����
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
            wav_S1_64(idx, (i-1)*oversamp+1:(i)*oversamp) = complex(I_sig, Q_sig);  %��������

        end
        wav_S1(idx, :) = wav_S1_64(idx, 4:4:end);

    end
    
    % ����S2
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
    
% 2Mbps B����ģʽ
elseif mode == 2
    
    num_pulses = 12;
    num_bits_pn = 24;
    pn = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0];  % 2Mbps B ģʽ��ͬ��ͷ��������
    % ����S1
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
    
    % ����S2
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
    
% 500K����ģʽ    
elseif mode == 3
    num_pulses = 48;
    num_bits_pn = 45;
    pn = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];  % 500Kbps ģʽ��ͬ��ͷS1��S2��������
    pn_2 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];  % 500Kbps ģʽ��ͬ��ͷS3��S4��������
    % ����S1��S3
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
    
    % ����S4��S2
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
    
% 250K����ģʽ    
else
    num_pulses = 48;
    num_bits_pn = 45;
    pn = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    pn_2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    % ����S1��S3
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
        wav_S1_temp(idx, :) = wav_S1_64(idx, 4:4:end); %������������ͬ

    end
    
    % ����S4��S2
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
        wav_S2_temp(idx, :) = wav_S2_64(idx, 4:4:end); %������������ͬ

    end
    
    wav_S1 = repmat(wav_S1_temp, [2,1]);  %�ظ�
    wav_S2 = repmat(wav_S2_temp, [2,1]);

end
    
    
    
    
