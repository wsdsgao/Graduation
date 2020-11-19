function [de_out, soft_info] = GMSK_viterbi(signal_recv_BB, num_bits_pulse, oversamp_BB, phi_init, bit_init, g)
%GMSK - Description
%
% Syntax: [de_out, soft_info] = GMSK_viterbi()
%
% Long description
% signal_recv_BB    基带信号
% num_bits_pulse    一帧数据的符号数
% oversamp_BB       基带过采样倍数
% g                 高斯成型滤波器
% phi_init          初始相位
% 
% de_out            未经LDPC译码输出
% soft_info         LDPC软信息输入

de_out_index = 1;
de_out = zeros(1, num_bits_pulse);
soft_info = zeros(1, num_bits_pulse);

viterbi_deep = 40;
path_record = cell(32, 2);
path_weight = cell(32, 2);
path_phi = zeros(32, 2);

signal_recv_BB = signal_recv_BB * exp(-1i*phi_init);

signal_recv = signal_recv_BB(1:oversamp_BB);

for q = 1:2
    init_path = [bit_init, rem(q,2)];
    pos = bin2dec(num2str(init_path))+1;
    path_record{pos,1} = init_path;
    [phi_last, I_sig, Q_sig, ~] = GMSK(2*init_path-1, 0, g);
    path_phi(pos, 1) = phi_last;
    path_weight_m = path_weight{pos, 1};
    path_weight_m = [path_weight_m, real(sum(signal_recv .* complex(I_sig, -Q_sig)))]; % 取最大值
    path_weight{pos, 1} = path_weight_m;
end  

for n = 2:num_bits_pulse
    
    signal_recv = signal_recv_BB(1+(n-1)*oversamp_BB:n*oversamp_BB);

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
            
            path_weight0 = real(sum(signal_recv .* complex(I_sig0, -Q_sig0)));
            path_weight1 = real(sum(signal_recv .* complex(I_sig1, -Q_sig1)));

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

    if n == num_bits_pulse
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

% 解预编码
soft_info = soft_predecode(soft_info);
de_out = decode_pre(de_out);
    
end