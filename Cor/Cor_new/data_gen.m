function [data, data_record] = data_gen(num_frame, mode)
%data_gen - Description
%
% Syntax: data = data_gen
%
% Long description

% To Do: 加入交织

num_data_frame = 1024;

switch mode
case 1
    I_ldpc_1_cpl = [];
    for i = 1:num_frame
        % 生成数据并记录
        I_single_1(i,:) = randi([0,1],[1,num_data_frame]);
        % LDPC编码
        I_ldpc_1_cpl = [I_ldpc_1_cpl; encoder(I_single_1', num_data_frame, num_data_frame*3)];
    end
    data = I_ldpc_1_cpl';
    data_record = I_single_1;
case 2
    I_ldpc_2_cpl = [];
    for i = 1:num_frame
        I_single_2(i,:) = randi([0,1],[1,num_data_frame]);
        I_ldpc_2 = encoder(I_single_2', num_data_frame, num_data_frame*3);
        I_ldpc_2_cpl = [I_ldpc_2_cpl;I_ldpc_2;I_ldpc_2];
    end
    data = I_ldpc_2_cpl';
    data_record = I_single_2;
case 3
    I_ldpc_3_cpl = [];
    for i = 1:num_frame
        I_single_3(i,:) = randi([0,1],[1,num_data_frame]);
        I_ldpc_3 = encoder(I_single_3', num_data_frame, num_data_frame*5);
        I_ldpc_3_cpl = [I_ldpc_3_cpl;I_ldpc_3;I_ldpc_3];
    end
    data = I_ldpc_3_cpl';
    data_record = I_single_3;
case 4
    I_ldpc_4_cpl = [];
    for i = 1:num_frame
        I_single_4(i,:) = randi([0,1],[1,num_data_frame]);
        I_ldpc_4 = encoder(I_single_4', num_data_frame, num_data_frame*5);  
        I_ldpc_4_cpl = [I_ldpc_4_cpl;I_ldpc_4;I_ldpc_4;I_ldpc_4;I_ldpc_4];
    end
    data = I_ldpc_4_cpl';
    data_record = I_single_4;
end