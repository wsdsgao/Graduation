function [pn_lib_S1, pn_lib_S2, pn_lib_S3, pn_lib_S4] = pn_gen

%  生成对应96个脉冲的同步头序列
%  每一个脉冲的同步头序列都不相同
%  同步头序列 每10ms更新一次

num_pulses = 48; 
num_bits_pn = 24;  % 同步头S1\S2长度
num_bits_pn_2 = 21;  % 同步头S3\S4长度

pn_lib_S1_temp = zeros(num_pulses, num_bits_pn);
for i = 1:num_pulses
    pn_lib_S1_temp(i,:) = randi([0,1], [1,num_bits_pn]);
end

pn_lib_S2_temp = zeros(num_pulses,num_bits_pn);
for i = 1:num_pulses
    pn_lib_S2_temp(i,:) = randi([0,1], [1,num_bits_pn]);
end

pn_lib_S3_temp = zeros(num_pulses,num_bits_pn_2);
for i = 1:num_pulses
    pn_lib_S3_temp(i,:) = randi([0,1], [1,num_bits_pn_2]);
end

pn_lib_S4_temp = zeros(num_pulses,num_bits_pn_2);
for i = 1:num_pulses
    pn_lib_S4_temp(i,:) = randi([0,1], [1,num_bits_pn_2]);
end

pn_lib_S1 = repmat(pn_lib_S1_temp, [2,1]);  %重复：48*2=96行

pn_lib_S2 = repmat(pn_lib_S2_temp, [2,1]);

pn_lib_S3 = repmat(pn_lib_S3_temp, [2,1]);

pn_lib_S4 = repmat(pn_lib_S4_temp, [2,1]);

