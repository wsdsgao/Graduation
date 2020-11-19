% close all;
% clear all;
% bits = randi(2,1,20);
% bits = bits - 1;

function predecode_out = decode_pre(precode_in)
%decode cription
%
% Syntax: predecode = decode()
%
% Long description

predecode_out = zeros(size(precode_in));

for i = 1:length(precode_in)
    if rem(i, 2)
        if i == 1
            predecode_out(i) = xor(precode_in(i), 0);
        else
            predecode_out(i) = xor(precode_in(i), predecode_out(i-1));
        end
    else
        predecode_out(i) = ~xor(precode_in(i), predecode_out(i-1)); 
    end
end
    
end