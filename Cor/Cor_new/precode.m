function precode_out = precode(bits)
%precode - Description
%
% Syntax: precode_out = precode(bits)
%
% Long description  
% close all;
% clear all;
% bits = randi(2,1,20);
% bits = bits - 1;


precode_out = zeros(size(bits));

for i = 1:length(bits)
    if rem(i, 2)
        if i == 1
            precode_out(i) = xor(bits(i), 0);
        else
            precode_out(i) = xor(bits(i), bits(i-1));
        end
    else
        precode_out(i) = ~xor(bits(i), bits(i-1)); 
    end
end

end