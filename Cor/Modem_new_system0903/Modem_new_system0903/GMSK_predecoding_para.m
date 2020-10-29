function [bits_decode_out] = GMSK_predecoding_para(bits_decode_in)
% input: 组帧后，单极性码

% output：输出预编码后的双极性码
           
%模块功能，进行GMSK预编码

% clear all;close all;
% sync_gen_poly=[0,0,1,1,0,0,0,0,0,0,0,0,1]; %x13+x4+x3+x+1
% sync_int_phase=[0,1,0,0,1,1,1,0,0,1,0,1,0];
% sync_m_seq=m_sequence( sync_gen_poly,sync_int_phase);
% bits_precode_in=sync_m_seq(1:20);

[mat_row, mat_col] = size(bits);

bits_decode_out=zeros(mat_row,length(bits_decode_in)-3);

ini_phase=rem(sum(bits_decode_in(1:3)),2);
ini_bit=1-ini_phase;

decode_Pro=zeros(1,length(bits_decode_in)-2);

decode_Pro(1)=ini_bit;
decode_Pro(1)=ini_phase;
% decode_Pro(1)=0;
% precode2=[0,bits_precode_in(2:2:end)];%取偶数
n=0;
for i=1:length(bits_decode_in)-3
    if n==0
        decode_Pro(i+1)= xor(decode_Pro(i),bits_decode_in(i+3));
        n=1;
    else
        decode_Pro(i+1)=~xor(decode_Pro(i),bits_decode_in(i+3));
        n=0;
    end

end
bits_decode_out=decode_Pro(2:end); 
% bits_precode_out(1:end)=2*precode_Pro-1;%单极性码转换成双极性码

% bits_precode_out(2:2:end)=2*precode2_Pro-1;%单极性码转换成双极性码
