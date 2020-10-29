function [bits_precode_out] = GMSK_precoding_para(bits_precode_in)
% input: 组帧后，单极性码

% output：输出预编码后的双极性码
           
%模块功能，进行GMSK预编码

% clear all;close all;
% sync_gen_poly=[0,0,1,1,0,0,0,0,0,0,0,0,1]; %x13+x4+x3+x+1
% sync_int_phase=[0,1,0,0,1,1,1,0,0,1,0,1,0];
% sync_m_seq=m_sequence( sync_gen_poly,sync_int_phase);
% bits_precode_in=sync_m_seq(1:20);

[mat_row, mat_col] = size(bits_precode_in);
bits_precode_pro=(bits_precode_in+1)/2;
% bits_precode_pro=bits_precode_in;
bits_precode_out=zeros(mat_row,mat_col);
precode=[zeros(1,mat_row),bits_precode_pro(1:end)];

% precode2=[0,bits_precode_in(2:2:end)];%取偶数
n=0;
for m=1:mat_row
    for i=1:length(precode)-1
        if n==0
            precode_Pro(m,i)= rem((precode(i)+precode(i+1)),2);
            n=1;
        else
            precode_Pro(m,i)=1-rem((precode(i)+precode(i+1)),2);
            n=0;
        end
    end
end
bits_precode_out=2*precode_Pro-1; 
% bits_precode_out(1:end)=2*precode_Pro-1;%单极性码转换成双极性码

% bits_precode_out(2:2:end)=2*precode2_Pro-1;%单极性码转换成双极性码
