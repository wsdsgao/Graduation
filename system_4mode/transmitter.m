function signal_trans = transmitter(bits, fh_pat_lib_1, th_pat_lib_1, fh_pat_lib_2, th_pat_lib_2, pn_lib_S1_1, pn_lib_S1_2,...
    pn_lib_S2_1, pn_lib_S2_2, pn_lib_S3_1, pn_lib_S3_2, pn_lib_S4_1, pn_lib_S4_2, mode)

% 载入参数
load('lib/f_trans.mat');  % 21个频点
load('lib/g_1024.mat');  % GMSK调制 g函数 

% 基本参数定义
bit_rate = 16e6;  % 符号速率
T = 1/bit_rate;  % 符号时间
fs_IF = 1024e6;  % 中频信号采样速率
oversamp_IF = T * fs_IF;
num_bits_pn = 24;  % 同步头S1\S2 长度
num_bits_pn_2 = 21;  % 同步头S3\S4 长度

switch mode
    % 2Mbps A 模式
    case 1 
        
        num_pulses = 12;
        
        % 同步头序列
        pn = [1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; % 2Mbps A 模式同步头S1\S2对应的扰码序列
        S1_lib_1 = zeros(1, num_bits_pn);
        S2_lib_1 = zeros(1, num_bits_pn);
        S1_lib_2 = zeros(1, num_bits_pn);
        S2_lib_2 = zeros(1, num_bits_pn);
        for i = 1:num_pulses
            S1_lib_1(i,:) = double(xor(pn_lib_S1_1(i,:), pn));
            S2_lib_1(i,:) = double(xor(pn_lib_S2_1(i,:), pn));
            S1_lib_2(i,:) = double(xor(pn_lib_S1_2(i,:), pn));
            S2_lib_2(i,:) = double(xor(pn_lib_S2_2(i,:), pn));
        end
        
        % 跳频图案
        fh_pat_1 = fh_pat_lib_1(1:num_pulses);  %其实都只是序号的排列
        fh_pat_2 = fh_pat_lib_2(1:num_pulses);
        % 跳时图案
        th_pat_1 = th_pat_lib_1(1:num_pulses);
        th_pat_2 = th_pat_lib_2(1:num_pulses); 
        
        
        [mat_row, mat_col] = size(bits);
        frame_idx = 0;
        for i = 1:mat_row

            % 每发送完一帧 重新选择一次跳频、跳时图案和同步头序列
            if (mod(i-1,12)==0)
                frame_idx = frame_idx + 1; 
                a = rand(1);
                if (a<0.5)
                    sel = 1;
                else
                    sel = 2;
                end

                if (sel == 1)
                    fh_pat = fh_pat_1;
                    th_pat = th_pat_1;
                    pn_lib_S1 = S1_lib_1;
                    pn_lib_S2 = S2_lib_1;
                elseif (sel == 2)
                    fh_pat = fh_pat_2;
                    th_pat = th_pat_2;
                    pn_lib_S1 = S1_lib_2;
                    pn_lib_S2 = S2_lib_2;               
                end

                th_pat_idx(frame_idx) = sel;
            end

            % 调制参数初始化
            % 初相、中频频率、同步头序列、跳时长度
            phi_last = 0;  % 起始相位
            f_idx = fh_pat(mod(i-1,12)+1);
            f_IF = f_trans(f_idx);
            preamble_S1 = 2*(pn_lib_S1(mod(i-1,12)+1,:))-1;
            preamble_S2 = 2*(pn_lib_S2(mod(i-1,12)+1,:))-1;
            t_mod = -304*T/2:1/fs_IF:304*T/2-1/fs_IF;
            t_mod = t_mod(1:end) + (T/oversamp_IF/2);  % 中心对称 
            bits_trans = [preamble_S1, bits(i,25:end-24), preamble_S2];  % 组包，双极性转换

            for j = 1:304   % 每一个脉冲: 304bit长度的波形

                % 对不同位置，取5bit数据，准备送入调制模块
                if (j == 1)  %同步头第一bit
                    bit_5 = [bits_trans(j), bits_trans(j), bits_trans(j:j+2)];
                elseif (j == 2)
                    bit_5 = [bits_trans(j-1), bits_trans(j-1:j+2)];
                elseif (j == 303)  % 尾同步头倒数第二bit
                    bit_5 = [bits_trans(j-2:j+1), 0];
                elseif (j == 304)  % 尾同步头最后一bit
                    bit_5 = [bits_trans(j-2:j), 0, 0];
                else
                    bit_5 = bits_trans(j-2:j+2);
                end

                [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, f_IF, phi_last, g);
                signal_trans_temp_BB(i, (j-1)*oversamp_IF+1:(j)*oversamp_IF) = complex(I_sig, Q_sig);
                phi_all(i, (j-1)*oversamp_IF+1:(j)*oversamp_IF) = phi_int;
            end
            I_sig_IF = cos(2*pi*f_IF*t_mod);
            Q_sig_IF = sin(2*pi*f_IF*t_mod);
            signal_trans_temp_IF(i, :) = signal_trans_temp_BB(i, :) .* complex(I_sig_IF, Q_sig_IF);   % 不连续相位调制，每一个脉冲的载波的起始相位都为0
%             signal_trans_temp_IF(i,:) = signal_trans_temp_BB(i,:);
        
            %画波形
            if(i==1)  %第1帧的第1个脉冲
                figure(1);
                plot(real(signal_trans_temp_IF(i, 1:5*oversamp_IF)));  %5bit
                legend(['载频为：', num2str(f_IF/1000000), 'MHz']);
                figure(2);
                plot(phi_all(i, 20*oversamp_IF+1:30*oversamp_IF));  %跨越前同步头和数据序列的10bit
                legend(['序列为：', num2str(bits_trans(21:30))]);
                figure(3);
                fs_BB = fs_IF/8;
                t_mod1 = t_mod(8:8:end)
                N1 = length(t_mod1);  %样点个数
                df1 = fs_BB/(N1-1) ;  %分辨率
                f1 = (0:N1-1)*df1/16e6;  %其中每点的频率
                plot(f1, 20*log10(real(fft(signal_trans_temp_BB(i, 8:8:end))/N1*2)));  %除以N1/2
                figure(4);
                N2 = length(t_mod);  %样点个数
                df2 = fs_IF/(N2-1) ;  %分辨率
                f2 = (0:N2-1)*df2;  %其中每点的频率
                plot(f2, 20*log10(real(fft(real(signal_trans_temp_IF(i, :)))/N2*2)));  %除以N2/2
            end

        end

        % 将各帧调制后的信号波形按照跳时图案组成连续的信号波形
        last = 0;
        frame_idx = 0;
        for i = 1:mat_row
            if (mod(i-1,12)==0)
                frame_idx = frame_idx + 1;
                th_idx = th_pat_idx(frame_idx);
                if th_idx==1
                    th_pat = th_pat_1;
                elseif th_idx==2
                    th_pat = th_pat_2;
                end
            end
            th = th_pat(mod(i-1,12)+1);
            temp_S = [zeros(1, th/2*oversamp_IF), signal_trans_temp_IF(i,:), zeros(1, th/2*oversamp_IF)];
            signal_trans(last+1:last+length(temp_S)) = temp_S;  %都排成一行
            last = last+length(temp_S);
        end

    % 2Mbps B 模式
    case 2
        
        num_pulses = 12;
        
        % 同步头序列
        pn = [0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0];  % 2Mbps B 模式同步头S1\S2对应的扰码序列
        S1_lib_1 = zeros(1, num_bits_pn);
        S2_lib_1 = zeros(1, num_bits_pn);
        S1_lib_2 = zeros(1, num_bits_pn);
        S2_lib_2 = zeros(1, num_bits_pn);
        for i = 1:num_pulses
            S1_lib_1(i,:) = double(xor(pn_lib_S1_1(i,:), pn));
            S2_lib_1(i,:) = double(xor(pn_lib_S2_1(i,:), pn));
            S1_lib_2(i,:) = double(xor(pn_lib_S1_2(i,:), pn));
            S2_lib_2(i,:) = double(xor(pn_lib_S2_2(i,:), pn));
        end
        
        % 跳频图案
        fh_pat_1 = fh_pat_lib_1(1:num_pulses);
        fh_pat_2 = fh_pat_lib_2(1:num_pulses);
        % 跳时图案
        th_pat_1 = th_pat_lib_1(1:num_pulses);
        th_pat_2 = th_pat_lib_2(1:num_pulses);
        
        
        [mat_row, mat_col] = size(bits);
        frame_idx = 0;
        frame_length = 304*12+512*6+103*12;
        for i = 1:mat_row

            % 每发送完一帧 重新选择一次跳频、跳时图案
            frame_idx = frame_idx + 1; 
            a = rand(1);
            if (a<0.5)
                sel = 1;
            else
                sel = 2;
            end

            if (sel == 1)
                fh_pat = fh_pat_1;
                th_pat = th_pat_1;
                pn_lib_S1 = S1_lib_1;
                pn_lib_S2 = S2_lib_1;
            elseif (sel == 2)
                fh_pat = fh_pat_2;
                th_pat = th_pat_2;
                pn_lib_S1 = S1_lib_2;
                pn_lib_S2 = S2_lib_2;               
            end

            th_pat_idx(frame_idx) = sel;
            frame_last_bit = 0;
            last = 0;
            num_bits_pulse = zeros(12);
            for j = 1:num_pulses

                % 调制参数初始化
                % 初相、中频频率、同步头序列、跳时长度
                phi_last = 0;  % 起始相位
                f_idx = fh_pat(mod(j-1,12)+1);
                f_IF = f_trans(f_idx);
                preamble_S1 = 2*(pn_lib_S1(mod(j-1,12)+1,:))-1;
                preamble_S2 = 2*(pn_lib_S2(mod(j-1,12)+1,:))-1;
                th = th_pat(mod(j-1,12)+1);
                bits_trans = [bits(i, frame_last_bit+1:frame_last_bit+th/2), preamble_S1, bits(i, frame_last_bit+th/2+1:frame_last_bit+th/2+256), preamble_S2, bits(i, frame_last_bit+th/2+257:frame_last_bit+th/2+256+th/2)]; 
                frame_last_bit = frame_last_bit + 256 + th;
                signal_trans_temp_BB = zeros(1, (304+th)*oversamp_IF);
                for jj = 1:304+th   % 每一个脉冲: (304+th)bit长度的波形

                    % 对不同位置，取5bit数据，准备送入调制模块

                    % 前数据段
                    if (jj == 1)
                        bit_5 = [bits_trans(jj), bits_trans(jj), bits_trans(jj:jj+2)];
                    elseif (jj == 2)
                        bit_5 = [bits_trans(jj-1), bits_trans(jj-1:jj+2)];
                    % elseif (jj == 3)
                    %     bit_5 = [1, 1, 1, bits_trans(1:2)];
                    % elseif (jj == 4)
                    %     bit_5 = [1, 1, bits_trans(1:3)];
                    % elseif (jj == 5)
                    %     bit_5 = [1, bits_trans(1:4)];
                    % elseif (jj < th/2+3+1)
                    %     bit_5 = bits_trans(jj-3-2:jj-3+2);
                    % elseif (jj == th/2+3-1)
                    %     bit_5 = [bits_trans(jj-3-2:jj-3+1),0];
                    % elseif (jj == th/2+3)
                    %     bit_5 = [bits_trans(jj-3-2:jj-3), 0, 0];


                    % 中间数据段    
                    % elseif (jj == th/2+3+1)  % 首同步头第1bit
                    %     phi_last = 0;
                    %     bit_5 = [bits_trans(jj-3), bits_trans(jj-3), bits_trans(jj-3:jj-3+2)];
                    % elseif (jj == th/2+3+2)  % 首同步头第2bit
                    %     bit_5 = [bits_trans(jj-3-1), bits_trans(jj-3-1:jj-3+2)];
                    % elseif (jj < th/2+3+280-1)
                    %     bit_5 = bits_trans(jj-3-2:jj-3+2);
                    % elseif (jj == th/2+3+280-1)  % 数据部分倒数第2bit
                    %     bit_5 = [bits_trans(jj-3-2:jj-3+1), 0];
                    % elseif (jj == th/2+3+280)  % 数据部分最后1bit
                    %     bit_5 = [bits_trans(jj-3-2:jj-3), 0, 0];                   


                    % 尾数据段     
                    % elseif (jj == th/2+3+280+1)  % 尾同步头第1bit 
                    %     phi_last = 0;
                    %     bit_5 = [bits_trans(jj-3),bits_trans(jj-3),bits_trans(jj-3:jj-3+2)];
                    % elseif (jj == th/2+3+280+2)  % 尾同步头第2bit 
                    %     bit_5 = [bits_trans(jj-3-1), bits_trans(jj-3-1:jj-3+2)];
                    elseif (jj < th+304-1)  
                        bit_5 = bits_trans(jj-2:jj+2);
                    elseif (jj == th+304-1)  % 数据段倒数第2bit
                        bit_5 = [bits_trans(jj-2:jj+1), 0];
                    elseif (jj == th+304)  % 数据段最后1bit
                        bit_5 = [bits_trans(jj-2:jj), 0, 0];

                    end

                    [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, f_IF, phi_last, g);
                    signal_trans_temp_BB((jj-1)*oversamp_IF+1:(jj)*oversamp_IF) = complex(I_sig, Q_sig);
                    phi_all(j, (jj-1)*oversamp_IF+1:(jj)*oversamp_IF) = phi_int;

                end
                num_bits_pulse(j) = th+304;
                t_mod = -num_bits_pulse(j)*T/2:1/fs_IF:num_bits_pulse(j)*T/2-1/fs_IF;
                t_mod = t_mod(1:end) + (T/oversamp_IF/2);  % 中心对称

                I_sig_IF = cos(2*pi*f_IF*t_mod);
                Q_sig_IF = sin(2*pi*f_IF*t_mod);
                signal_trans_temp_IF = signal_trans_temp_BB .* complex(I_sig_IF, Q_sig_IF);  % 不连续相位调制，每一个脉冲的载波的起始相位都为0
%                 signal_trans_temp_IF = signal_trans_temp_BB;

                %画波形
                if(i==1&&j==1)
                    figure(1);
                    plot(real(signal_trans_temp_IF(1:5*oversamp_IF)));
                    legend(['载频为：', num2str(f_IF/1000000), 'MHz']);
                    figure(2);
                    plot(phi_all(j, 1:10*oversamp_IF));
                    legend(['序列为：', num2str(bits_trans(1:10))]);
                    figure(3);
                    N = length(t_mod);  %样点个数
                    df = fs_IF/(N-1) ;  %分辨率
                    f = (0:N-1)*df/16e6;  %其中每点的频率
                    plot(f, 20*log10(real(fft(real(signal_trans_temp_IF(:)))/N*2)));
                end


                temp = [zeros(1, 103*oversamp_IF), signal_trans_temp_IF];  %103bit的固定间隔，但A并没加
                signal_trans_temp(last+1:last+(103+num_bits_pulse(j))*oversamp_IF) = temp;
                last = last + (103+num_bits_pulse(j))*oversamp_IF;
            end
            signal_trans((i-1)*frame_length*oversamp_IF+1:i*frame_length*oversamp_IF) = signal_trans_temp;
        end
           
    % 500K 模式    
    case 3
        
        num_pulses = 48;
        
        % 同步头序列
        pn = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];  % 500Kbps 模式 同步头S1\S2对应的扰码序列
        pn_2 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];  % 500Kbps 模式 同步头S3\S4对应的扰码序列
        S1_lib_1 = zeros(1, num_bits_pn);
        S2_lib_1 = zeros(1, num_bits_pn);
        S3_lib_1 = zeros(1, num_bits_pn_2);
        S4_lib_1 = zeros(1, num_bits_pn_2);
        S1_lib_2 = zeros(1, num_bits_pn);
        S2_lib_2 = zeros(1, num_bits_pn);
        S3_lib_2 = zeros(1, num_bits_pn_2);
        S4_lib_2 = zeros(1, num_bits_pn_2);
        for i = 1:num_pulses
            S1_lib_1(i,:) = double(xor(pn_lib_S1_1(i,:), pn));
            S2_lib_1(i,:) = double(xor(pn_lib_S2_1(i,:), pn));
            S3_lib_1(i,:) = double(xor(pn_lib_S3_1(i,:), pn_2));
            S4_lib_1(i,:) = double(xor(pn_lib_S4_1(i,:), pn_2));
            S1_lib_2(i,:) = double(xor(pn_lib_S1_2(i,:), pn));
            S2_lib_2(i,:) = double(xor(pn_lib_S2_2(i,:), pn));
            S3_lib_2(i,:) = double(xor(pn_lib_S3_2(i,:), pn_2));
            S4_lib_2(i,:) = double(xor(pn_lib_S4_2(i,:), pn_2));
        end
        
        % 跳频图案
        fh_pat_1 = fh_pat_lib_1(1:num_pulses);
        fh_pat_2 = fh_pat_lib_2(1:num_pulses);
        % 跳时图案
        th_pat_1 = th_pat_lib_1(1:num_pulses);
        th_pat_2 = th_pat_lib_2(1:num_pulses);
        
        [mat_row, mat_col] = size(bits);
        frame_idx = 0;
        for i = 1:mat_row

            % 每发送完一帧 重新选择一次跳频、跳时图案
            if (mod(i-1,num_pulses)==0)
                frame_idx = frame_idx + 1; 
                a = rand(1);
                if (a<0.5)
                    sel = 1;
                else
                    sel = 2;
                end

                if (sel == 1)
                    fh_pat = fh_pat_1;
                    th_pat = th_pat_1;
                    pn_lib_S1 = S1_lib_1;
                    pn_lib_S2 = S2_lib_1;
                    pn_lib_S3 = S3_lib_1;
                    pn_lib_S4 = S4_lib_1;
                elseif (sel == 2)
                    fh_pat = fh_pat_2;
                    th_pat = th_pat_2;
                    pn_lib_S1 = S1_lib_2;
                    pn_lib_S2 = S2_lib_2;
                    pn_lib_S3 = S3_lib_2;
                    pn_lib_S4 = S4_lib_2;
                end

                th_pat_idx(frame_idx) = sel;
            end

            % 调制参数初始化
            % 初相、中频频率、同步头序列、跳时长度
            phi_last = 0;  % 起始相位
            f_idx = fh_pat(mod(i-1,num_pulses)+1);
            f_IF = f_trans(f_idx);
            preamble_S1 = 2*(pn_lib_S1(mod(i-1,num_pulses)+1,:))-1;
            preamble_S2 = 2*(pn_lib_S2(mod(i-1,num_pulses)+1,:))-1;
            preamble_S3 = 2*(pn_lib_S3(mod(i-1,num_pulses)+1,:))-1;
            preamble_S4 = 2*(pn_lib_S4(mod(i-1,num_pulses)+1,:))-1;
            t_mod = -304*T/2:1/fs_IF:304*T/2-1/fs_IF;
            t_mod = t_mod(1:end) + (T/oversamp_IF/2);  % 中心对称
            bits_trans = [preamble_S1, preamble_S3, bits(i,46:end-45), preamble_S4, preamble_S2];  % 组包，双极性转换

            for j = 1:304   % 每一帧: 304bit长度的波形

                % 对不同位置，取5bit数据，准备送入调制模块
                if (j == 1)
                    bit_5 = [bits_trans(j), bits_trans(j), bits_trans(j:j+2)];
                elseif (j == 2)
                    bit_5 = [bits_trans(j-1), bits_trans(j-1:j+2)];
                % elseif (j == 258)  % 数据部分倒数第二比特
                %     bit_5 = [bits_trans(j-2:j+1), 0];
                % elseif (j == 259)  % 数据部分最后1比特
                %     bit_5 = [bits_trans(j-2:j), 0, 0];
                % elseif (j == 260)  % 尾同步头第1bit
                    % phi_last = 0;
                    % bit_5 = [bits_trans(j), bits_trans(j), bits_trans(j:j+2)];
                % elseif (j == 261)  % 尾同步头第2bit
                %     bit_5 = [bits_trans(j-1), bits_trans(j-1:j+2)];
                elseif (j == 303)  % 尾同步头倒数第二bit
                    bit_5 = [bits_trans(j-2:j+1), 0];
                elseif (j == 304)  % 尾同步头最后一bit
                    bit_5 = [bits_trans(j-2:j), 0, 0];
                else
                    bit_5 = bits_trans(j-2:j+2);
                end

                [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, f_IF, phi_last, g);
                signal_trans_temp_BB(i, (j-1)*oversamp_IF+1:(j)*oversamp_IF) = complex(I_sig, Q_sig);
                phi_all(i, (j-1)*oversamp_IF+1:(j)*oversamp_IF) = phi_int;

            end
            I_sig_IF = cos(2*pi*f_IF*t_mod);
            Q_sig_IF = sin(2*pi*f_IF*t_mod);
            signal_trans_temp_IF(i, :) = signal_trans_temp_BB(i, :) .* complex(I_sig_IF, Q_sig_IF);  % 不连续相位调制，每一个脉冲的载波起始相位都为0
%             signal_trans_temp_IF(i,:) = signal_trans_temp_BB(i,:);

            %画波形
            if(i==1)
                figure(1);
                plot(real(signal_trans_temp_IF(i, 1:5*oversamp_IF)));
                legend(['载频为：', num2str(f_IF/1000000), 'MHz']);
                figure(2);
                plot(phi_all(i, 1:10*oversamp_IF));
                legend(['序列为：', num2str(bits_trans(1:10))]);
                figure(3);
                N = length(t_mod);  %样点个数
                df = fs_IF/(N-1) ;  %分辨率
                f = (0:N-1)*df/16e6;  %其中每点的频率
                plot(f, 20*log10(real(fft(real(signal_trans_temp_IF(i, :)))/N*2)));
            end

        end


        % 将各帧调制后的信号波形按照跳时图案组成连续的信号波形
        last = 0;
        frame_idx = 0;
        for i = 1:mat_row
            if (mod(i-1,num_pulses)==0)
                frame_idx = frame_idx + 1;
                th_idx = th_pat_idx(frame_idx);
                if th_idx==1
                    th_pat = th_pat_1;
                elseif th_idx==2
                    th_pat = th_pat_2;
                end
            end
            if (mod(i,num_pulses)~=0)  %一帧的末尾
                th = th_pat(mod(i-1,num_pulses)+1);
                temp_S = [zeros(1, th/2*oversamp_IF), signal_trans_temp_IF(i,:), zeros(1, th/2*oversamp_IF)];
                signal_trans(last+1:last+length(temp_S)) = temp_S;
                last = last+length(temp_S);
            else
                th = th_pat(mod(i-1,num_pulses)+1);
                temp_S = [zeros(1, th/2*oversamp_IF), signal_trans_temp_IF(i,:), zeros(1, th/2*oversamp_IF), zeros(1, 10*oversamp_IF)];  %这个10bit是哪的？
                signal_trans(last+1:last+length(temp_S)) = temp_S;
                last = last+length(temp_S);
            end
        end
        
    % 250K 模式    
    case 4
        
        num_pulses = 96;
        
        % 同步头序列
        pn = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];  % 250Kbps 模式 同步头S1\S2对应的扰码序列
        pn_2 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];  % 250Kbps 模式 同步头S3\S4对应的扰码序列
        S1_lib_1 = zeros(1, num_bits_pn);
        S2_lib_1 = zeros(1, num_bits_pn);
        S3_lib_1 = zeros(1, num_bits_pn_2);
        S4_lib_1 = zeros(1, num_bits_pn_2);
        S1_lib_2 = zeros(1, num_bits_pn);
        S2_lib_2 = zeros(1, num_bits_pn);
        S3_lib_2 = zeros(1, num_bits_pn_2);
        S4_lib_2 = zeros(1, num_bits_pn_2);
        for i = 1:num_pulses
            S1_lib_1(i,:) = double(xor(pn_lib_S1_1(i,:), pn));
            S2_lib_1(i,:) = double(xor(pn_lib_S2_1(i,:), pn));
            S3_lib_1(i,:) = double(xor(pn_lib_S3_1(i,:), pn_2));
            S4_lib_1(i,:) = double(xor(pn_lib_S4_1(i,:), pn_2));
            S1_lib_2(i,:) = double(xor(pn_lib_S1_2(i,:), pn));
            S2_lib_2(i,:) = double(xor(pn_lib_S2_2(i,:), pn));
            S3_lib_2(i,:) = double(xor(pn_lib_S3_2(i,:), pn_2));
            S4_lib_2(i,:) = double(xor(pn_lib_S4_2(i,:), pn_2));
        end

        % 跳频图案
        fh_pat_1 = fh_pat_lib_1(1:num_pulses);
        fh_pat_2 = fh_pat_lib_2(1:num_pulses);
        % 跳时图案
        th_pat_1 = th_pat_lib_1(1:num_pulses);
        th_pat_2 = th_pat_lib_2(1:num_pulses); 
        
        [mat_row, mat_col] = size(bits);
        
        frame_idx = 0;
        for i = 1:mat_row

            % 每发送完一帧 重新选择一次跳频、跳时图案
            if (mod(i-1,num_pulses)==0)
                frame_idx = frame_idx + 1; 
                a = rand(1);
                if (a<0.5)
                    sel = 1;
                else
                    sel = 2;
                end

                if (sel == 1)
                    fh_pat = fh_pat_1;
                    th_pat = th_pat_1;
                    pn_lib_S1 = S1_lib_1;
                    pn_lib_S2 = S2_lib_1;
                    pn_lib_S3 = S3_lib_1;
                    pn_lib_S4 = S4_lib_1;
                elseif (sel == 2)
                    fh_pat = fh_pat_2;
                    th_pat = th_pat_2;
                    pn_lib_S1 = S1_lib_2;
                    pn_lib_S2 = S2_lib_2;
                    pn_lib_S3 = S3_lib_2;
                    pn_lib_S4 = S4_lib_2;
                end

                th_pat_idx(frame_idx) = sel;
            end

            % 调制参数初始化
            % 初相、中频频率、同步头序列、跳时长度
            phi_last = 0;  % 起始相位
            f_idx = fh_pat(mod(i-1,num_pulses)+1);
            f_IF = f_trans(f_idx);
            preamble_S1 = 2*(pn_lib_S1(mod(i-1,num_pulses)+1,:))-1;
            preamble_S2 = 2*(pn_lib_S2(mod(i-1,num_pulses)+1,:))-1;
            preamble_S3 = 2*(pn_lib_S3(mod(i-1,num_pulses)+1,:))-1;
            preamble_S4 = 2*(pn_lib_S4(mod(i-1,num_pulses)+1,:))-1;
            t_mod = -304*T/2:1/fs_IF:304*T/2-1/fs_IF;
            t_mod = t_mod(1:end) + (T/oversamp_IF/2);  % 中心对称
            bits_trans = [preamble_S1, preamble_S3, bits(i,46:end-45), preamble_S4, preamble_S2];  % 组包，双极性转换

            for j = 1:304   % 每一帧: 304bit长度的波形

                % 对不同位置，取5bit数据，准备送入调制模块
                if (j == 1)
                    bit_5 = [bits_trans(j), bits_trans(j), bits_trans(j:j+2)];
                elseif (j == 2)
                    bit_5 = [bits_trans(j-1), bits_trans(j-1:j+2)];
                % elseif (j == 258)  % 数据部分倒数第二比特
                %     bit_5 = [bits_trans(j-2:j+1), 0];
                % elseif (j == 259)  % 数据部分最后1比特
                %     bit_5 = [bits_trans(j-2:j), 0, 0];
                % elseif (j == 260)  % 尾同步头第1bit
                %     phi_last = 0;
                %     bit_5 = [bits_trans(j), bits_trans(j), bits_trans(j:j+2)];
                % elseif (j == 261)  % 尾同步头第2bit
                %     bit_5 = [bits_trans(j-1), bits_trans(j-1:j+2)];
                elseif (j == 303)  % 尾同步头倒数第二bit
                    bit_5 = [bits_trans(j-2:j+1), 0];
                elseif (j == 304)  % 尾同步头最后一bit
                    bit_5 = [bits_trans(j-2:j), 0, 0];
                else
                    bit_5 = bits_trans(j-2:j+2);
                end

                [phi_last, I_sig, Q_sig, phi_int] = GMSK(bit_5, f_IF, phi_last, g);
                signal_trans_temp_BB(i, (j-1)*oversamp_IF+1:(j)*oversamp_IF) = complex(I_sig, Q_sig);
                phi_all(i, (j-1)*oversamp_IF+1:(j)*oversamp_IF) = phi_int;

            end
            I_sig_IF = cos(2*pi*f_IF*t_mod);
            Q_sig_IF = sin(2*pi*f_IF*t_mod);
            signal_trans_temp_IF(i, :) = signal_trans_temp_BB(i, :) .* complex(I_sig_IF, Q_sig_IF);  % 不连续相位调制，每一个脉冲的载波起始相位都为0
%             signal_trans_temp_IF(i,:) = signal_trans_temp_BB(i,:);

            %画波形
            if(i==1)
                figure(1);
                plot(real(signal_trans_temp_IF(i, 1:5*oversamp_IF)));
                legend(['载频为：', num2str(f_IF/1000000), 'MHz']);
                figure(2);
                plot(phi_all(i, 1:10*oversamp_IF));
                legend(['序列为：', num2str(bits_trans(1:10))]);
                figure(3);
                N = length(t_mod);  %样点个数
                df = fs_IF/(N-1) ;  %分辨率
                f = (0:N-1)*df/16e6;  %其中每点的频率
                plot(f, 20*log10(real(fft(real(signal_trans_temp_IF(i, :)))/N*2)));
            end

        end


        % 将各帧调制后的信号波形按照跳时图案组成连续的信号波形
        last = 0;
        frame_idx = 0;
        for i = 1:mat_row
            if (mod(i-1,num_pulses)==0)
                frame_idx = frame_idx + 1;
                th_idx = th_pat_idx(frame_idx);
                if th_idx==1
                    th_pat = th_pat_1;
                elseif th_idx==2
                    th_pat = th_pat_2;
                end
            end
            th = th_pat(mod(i-1,num_pulses)+1);
            temp_S = [zeros(1, th/2*oversamp_IF), signal_trans_temp_IF(i,:), zeros(1, th/2*oversamp_IF)];
            signal_trans(last+1:last+length(temp_S)) = temp_S;
            last = last+length(temp_S);
        end
        
        
end