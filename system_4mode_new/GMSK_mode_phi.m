function phi_out = GMSK_mode_phi(I, num, oversamp, g)
%myFun - Description
%
% Syntax: [signal_trans_BB] = myFun(input)
%
% Long description

% signal_trans_BB = zeros(1, oversamp*num);
bit_5 = zeros(1,5);
phi_last = 0;

for i = 1:num
    if i == 1
        bit_5 = [-1,-1,I(i:i+2)];
    elseif i == 2
        bit_5 = [-1,I(i-1:i+2)];
    elseif i == num-1
        bit_5 = [I(i-2:i+1),-1];
    elseif i == num
        bit_5 = [I(i-2:i),-1,-1];
    else
        bit_5 = I(i-2:i+2);
    end

    [phi_last, I_sig, Q_sig, phi_out] = GMSK(bit_5, phi_last, g);
    % signal_trans_BB((i-1)*oversamp+1:(i)*oversamp) = complex(I_sig, Q_sig);
    % phi_all((i-1)*oversamp_IF+1:(i)*oversamp_IF) = phi_int;
end
phi_out = phi_out(end);
    
end