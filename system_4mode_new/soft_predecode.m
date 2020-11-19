% close all;
% clear all;
% bits = randi(2,1,20);
% bits = bits - 1;

function predecode_out = soft_predecode(precode_in)
    %decode cription
    %
    % Syntax: predecode = decode()
    %
    % Long description
    
    predecode_out = zeros(size(precode_in));
    
    for i = 1:length(precode_in)
        if rem(i, 2)
            if i == 1
                predecode_out(i) = precode_in(i)*2;
            else
                if precode_in(i).*predecode_out(i-1) > 0
                    predecode_out(i) = abs(precode_in(i)) + abs(precode_in(i-1));
                else
                    predecode_out(i) = -abs(precode_in(i)) - abs(precode_in(i-1));
                end
            end
        else
            if precode_in(i).*predecode_out(i-1) > 0
                predecode_out(i) = -abs(precode_in(i)) - abs(precode_in(i-1));
            else
                predecode_out(i) = abs(precode_in(i)) + abs(precode_in(i-1));
            end
        end
    end
        
end