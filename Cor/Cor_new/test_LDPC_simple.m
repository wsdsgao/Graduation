function [suctt,bitErrorRate]=test_LDPC_simple(matrix,base_matrix,K,L,N,Z,encodeform,snr)
%file_name = 'ram_ring3.txt';
%[degree,q,K,decode_list,parity_list] = convert_file(file_name);

fpga=0;
qbit=6;
rtimes = 15;


verilator=K * 2 + 1;

suctt=0;

hls=0;

suct = 0;
biterror = 0;
for kt = 1:100000
    info_block=randi([0,1],K,1);

    encode2 = encoder(info_block,K,N);
    %parity_bits2 = [info_block(2*L+1:K,1);parity_bits];


    encode = [info_block(1:K,1);encode2];
    %encode2 = [info_block(2*L+1:K,1);parity_bits];

   

    sigma = 10^(-snr/20);
    mode = 2;
    n0 = (randn(1,N/mode)+randn(1,N/mode)*1j) * sigma/sqrt(2);
    constel = mod_APSK(mode,reshape(encode2,mode,N/mode));
    signal = constel + n0;

    soft_bits = demod_APSK_VHDL(mode,real(signal),imag(signal));
    logProb = reshape(soft_bits,N,1);

if fpga==1
    
    soft_bits = soft_bits * 4/3;
    
    %logProb = round((realpow(2,qbit-3)-1)*reshape(soft_bits,N,1));
    logProb = floor((realpow(2,qbit-2))*reshape(soft_bits,N,1));
    %logProb(find(logProb>realpow(2,qbit-1)-1)) = realpow(2,qbit-1)-1;
    logProb(logProb>realpow(2,qbit-1)-1) = realpow(2,qbit-1)-1;
    %logProb(find(logProb< -realpow(2,qbit-1)+1)) = -realpow(2,qbit-1)+1;
    logProb(logProb< -realpow(2,qbit-1)) = -realpow(2,qbit-1);
    
    logProb=[zeros(2*L,1);logProb(1:K-2*L,:);ones(K0-K,1)*31;logProb(K-2*L+1:N,:)];
    
else
    %logProb=[zeros(2*L,1);logProb(1:K-2*L,:);ones(K0-K,1);logProb(K-2*L+1:N,:)];
    %logProb=logProb/sigma/sigma*2;
end;
%N=N+2*L;



[decode_bits,suc] = Ldecoder2(logProb,K,N);

info_bits = decode_bits(1:K);

for i=1:K
    %fprintf('%d',decode_bits(i));
    if decode_bits(i)~=encode(i)
        biterror = biterror + 1;
        fprintf('e:%d,%d; ',i,encode(i));
    end;
end;
%fprintf('\n');


%N=N-2*L;

err2 = sum(abs(info_bits-info_block(1:K,1)))/K;
if err2==0
    suctt=suctt+1;
end;

bitErrorRate = biterror/kt/K;

st = 'error rate of iteration of '+string(kt)+' is '+string(err2)+' nowerror '+string(kt-suctt)+' bit error rate '+string(bitErrorRate);
fprintf('%d,%d,%d',kt,err2,suc);
disp(st);

suct = suct + suc;
if kt > suctt+1000
    break;
end;
end
%imshow(reshape((decode_bits-encode'),L,kmax)');
%imshow(reshape(encode,L,kmax)');
disp(suct);
disp(suctt);
return;


