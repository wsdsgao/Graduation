function encode_out = encoder(input_block,K,N)
%M = N * rate;
lift=[2, 4, 8, 16, 32, 64, 128, 256,3, 6, 12, 24, 48, 96, 192, 384,5, 10, 20, 40, 80, 160, 320,7, 14, 28, 56, 112, 224,9, 18, 36, 72, 144, 288,11, 22, 44, 88, 176, 352,13, 26, 52, 104, 208,15, 30, 60, 120, 240];
%L=7;
%N=280;

if K > 3840 || K > (N*4/5)
    [~,matrix,bg]=bg1();
    encodeform=0;
else
    [~,matrix,bg]=bg2();
    encodeform=3;
end;
%if K<3840
%    [~,matrix,bg]=bg2();
%    encodeform=3;
%else
%    [~,matrix,bg]=bg1();
%    encodeform=0;
%end;

if K > 3840 || K > (N*4/5)
      Z = 22;
elseif K > 640
      Z = 10;
elseif K > 560
      Z = 9;
      matrix=[matrix(1:42,1:9),matrix(1:42,11:52)];
      bg=[bg(1:42,1:9),bg(1:42,11:52)];  
elseif K > 192
      Z = 9;%8;
      matrix=[matrix(1:42,1:8),matrix(1:42,11:52)];
      bg=[bg(1:42,1:8),bg(1:42,11:52)];  
else 
      Z = 9;%6;
      matrix=[matrix(1:42,1:9),matrix(1:42,11:52)];
      bg=[bg(1:42,1:9),bg(1:42,11:52)];  
end;
 
liftingTemp = K/Z;

minDistance = 384;
minJ = 0;
for j = 1 : size(lift,2)
    if liftingTemp <= lift(j)
        distance = lift(j) - liftingTemp;
        if distance < minDistance
          minDistance = distance;
          minJ = j;
        end;
    end;
end;
    L = lift(minJ);

[degree,~,K0,list,~,~] = convert_matrix(matrix,bg,L,K,N,Z);


    input_block=[input_block;zeros(K0-K,1)];

M = N-K+2*L;
parity_bits = zeros(M,1);
parity_info=zeros(M,1);

for x=1:M
    for m=1:degree
        %if list(x,m)<=112
        if list(x,m)<=K
            parity_info(x)=xor(parity_info(x),input_block(list(x,m)));
        end;
    end;
end;

% temp=xor(xor(parity_info(1:L),parity_info(L+1:2*L)),xor(parity_info(2*L+1:3*L),parity_info(3*L+1:4*L)));
% temp1=[temp(2:L);temp(1)];
% parity_bits(1:L)=temp;
% %temp=[temp(L);temp(1:L-1)];
% parity_bits(L+1:2*L)=xor(temp1(1:L),parity_info(1:L));
% parity_bits(2*L+1:3*L)=xor(parity_bits(L+1:2*L),parity_info(L+1:2*L));
% parity_bits(3*L+1:4*L)=xor(temp1(1:L),parity_info(3*L+1:4*L));

temp=xor(xor(parity_info(1:L),parity_info(L+1:2*L)),xor(parity_info(2*L+1:3*L),parity_info(3*L+1:4*L)));
temp1=[temp(2:L);temp(1)];
temp2=[temp(L);temp(1:L-1)];
if encodeform==0
    parity_bits(1:L)=temp;
    parity_bits(L+1:2*L)=xor(temp1(1:L),parity_info(1:L));
    parity_bits(2*L+1:3*L)=xor(xor(temp1(1:L),parity_info(2*L+1:3*L)),parity_info(3*L+1:4*L));
    %parity_bits(2*L+1:3*L)=xor(xor(temp1(1:L),parity_info(1:L)),parity_info(1*L+1:2*L));
    parity_bits(3*L+1:4*L)=xor(temp1(1:L),parity_info(3*L+1:4*L));  
elseif encodeform==1
    parity_bits(1:L)=temp;
    parity_bits(L+1:2*L)=xor(temp1(1:L),parity_info(1:L));
    parity_bits(2*L+1:3*L)=xor(xor(temp1(1:L),parity_info(2*L+1:3*L)),parity_info(3*L+1:4*L));
    %parity_bits(2*L+1:3*L)=xor(xor(temp1(1:L),parity_info(1:L)),parity_info(1*L+1:2*L));
    parity_bits(3*L+1:4*L)=xor([temp1(2:L),temp1(1)],parity_info(3*L+1:4*L));
elseif encodeform==2
    parity_bits(1:L)=temp2;
    parity_bits(L+1:2*L)=xor(temp2(1:L),parity_info(1:L));
    %parity_bits(2*L+1:3*L)=xor(xor(temp2(1:L),parity_info(2*L+1:3*L)),parity_info(3*L+1:4*L));
    parity_bits(2*L+1:3*L)=xor(xor(temp2(1:L),parity_info(1:L)),parity_info(1*L+1:2*L));
    parity_bits(3*L+1:4*L)=xor(temp2(1:L),parity_info(3*L+1:4*L));
elseif encodeform==3
    parity_bits(1:L)=temp;
    parity_bits(L+1:2*L)=xor(temp1(1:L),parity_info(1:L));
    %parity_bits(2*L+1:3*L)=xor(xor(temp1(1:L),parity_info(2*L+1:3*L)),parity_info(3*L+1:4*L));
    parity_bits(2*L+1:3*L)=xor(xor(temp1(1:L),parity_info(1:L)),parity_info(1*L+1:2*L));
    parity_bits(3*L+1:4*L)=xor(temp1(1:L),parity_info(3*L+1:4*L));
else
    parity_bits(1:L)=temp2;
    parity_bits(L+1:2*L)=xor(temp2(1:L),parity_info(1:L));
    %parity_bits(2*L+1:3*L)=xor(xor(temp2(1:L),parity_info(2*L+1:3*L)),parity_info(3*L+1:4*L));
    parity_bits(2*L+1:3*L)=xor(xor(temp2(1:L),parity_info(1:L)),parity_info(1*L+1:2*L));
    parity_bits(3*L+1:4*L)=xor(temp2(1:L),parity_info(3*L+1:4*L));
end;

for x=4*L+1:M
    for m=1:degree
        if list(x,m)>K0 && list(x,m)<=K0+4*L
            parity_info(x)=xor(parity_info(x),parity_bits(list(x,m)-K0));
        end;
    end;
end;
parity_bits=[parity_bits(1:4*L);parity_info(4*L+1:M)];

% parity_bits(1)=parity_info(1);
% for mmm=0:1
% parity_bits(L+1)=mmm;
% parity_bits(2*L+1)=xor(parity_info(L+1),parity_bits(L+1));
% parity_bits(3*L+1)=xor(parity_info(3*L+1),parity_bits(1));
% parity_bits(L)=xor(xor(parity_info(2*L+1),parity_bits(3*L+1)),parity_bits(2*L+1));
% for i=L:-1:2
%     parity_bits(L+i)=xor(parity_info(i),parity_bits(i));
%     parity_bits(2*L+i)=xor(parity_info(L+i),parity_bits(L+i));
%     parity_bits(3*L+i)=xor(parity_info(3*L+i),parity_bits(i));
%     parity_bits(i-1)=xor(xor(parity_info(2*L+i),parity_bits(3*L+i)),parity_bits(2*L+i));
% end;
% input_bloc=[input_block;parity_bits(1:4*L)];
% for x=4*L+1:M
%     for m=1:degree
%         if list(x,m)<=K+4*L
%             parity_bits(x)=xor(parity_bits(x),input_bloc(list(x,m)));
%         end;
%     end;
% end;
% if parity_bits(1)==parity_info(1)
%     break;
% end;
% end;

parity_info=zeros(M,1);

check=0;
for x=1:M
    for m=1:degree
        if list(x,m)<=K0
            parity_info(x)=xor(parity_info(x),input_block(list(x,m)));
        elseif list(x,m)<=N+2*L+K0-K
            parity_info(x)=xor(parity_info(x),parity_bits(list(x,m)-K0));
        end;
    end;
    if parity_info(x)==1
        a=intersect(find(list(x,:)>K),find(list(x,:)<=N+2*L));
        %parity_bits(list(x,a(randi(size(a))))-K)=abs(1-parity_bits(list(x,a(randi(size(a))))-K));
        check=1;
        %disp(x);
    end;
end;


    encode_out = [input_block(2*L+1:K,1);parity_bits];

%fprintf('%d',check);

function [degree,q,K,list,decode_list,parity_list] = convert_matrix(matrix,base_matrix,L,Ko,N,Z)
    q=length(matrix(:,1));
    M=N-Ko+2*L;
    K=Z*L;
    
    matrix=matrix';
    parity_list = cell(length(matrix(:,1))-length(matrix(1,:)),1);
    for x=1:length(matrix(:,1))-length(matrix(1,:))
        pa=0;
        for y=1:length(matrix(1,:))
            if matrix(x,y)~=0
                pa=[pa,(matrix(x,y))*q+y];
            end;
        end;
        parity_list(x)={pa(1,2:length(pa))};
    end;
    
    decode_list = cell(length(matrix(1,:)),1);
    for y=1:length(matrix(1,:))
        pa=0;
        for x=1:length(matrix(:,1))-length(matrix(1,:))
            if matrix(x,y)~=0
                %pa=[pa,(x-1)*L+L-matrix(x,y)];
                pa=[pa,(x)*L-matrix(x,y)];
            end;
        end;
        decode_list(y)={pa(1,2:length(pa))};
    end;
    
    decode_list{1} = [decode_list{1},K,N-1];
    for i = 2:q
        decode_list{i} = [decode_list{i},K+(i-2)*L,K+(i-1)*L];
    end
    
    degree=0;
    for i=1:length(matrix(1,:))
        degree=max(degree,sum(base_matrix(i,:)));
    end;
    
    matrix=matrix';
    list=99999*ones(M,degree);
    for x=1:length(matrix(:,1))
        for y=1:length(matrix(1,:))
            for m=1:L
                if matrix(x,y)~=-1
                    xi=(x-1)*L+m;
                    %out=(y-1)*L+mod(L-matrix(x,y)+m-1,L)+1;
                    out=(y-1)*L+mod(matrix(x,y)+m-1,L)+1;
                    yi=sum(matrix(x,1:y)~=-1);
                    if xi<=M
                        list(xi,yi)=out;
                    end;
                end;
            end;
        end;
    end;