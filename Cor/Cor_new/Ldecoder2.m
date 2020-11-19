function[decode,suc] =  Ldecoder2(input_block,K,N)

    iteration = 15;

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

input_block=[zeros(2*L,1);input_block(1:K-2*L,:);ones(K0-K,1)*Inf;input_block(K-2*L+1:N,:)];

N = N + 2*L;

M = N-K;
N=size(input_block,1);

fpga=0;
qbit=8;

if fpga==1
    maxm=realpow(2,qbit-1)-1;
    max31=realpow(2,qbit-2)-1;
else
    maxm=Inf;
    max31=Inf;
    %maxm=9999;
end;
%mvc(i) = ones(degree)* maxm;
mvc = zeros(M,degree);
mcv = zeros(M,degree);
%mcvt=zeros(M,degree);
suc=1;

for it=1:iteration
    
    if fpga==1
        maxmcv=0;
        for maxCount = 1:M
            if mod(maxCount-1,L) == 0
                if max(abs(mcv(maxCount,:)))>=realpow(2,qbit-3)+realpow(2,qbit-4)
                %if max(abs(mcv(maxCount,:)))>=realpow(2,qbit-4)+realpow(2,qbit-5)
                    maxmcv=maxmcv + 1;
                end;
            end;
        end;
        %if max(abs(mcv(:,:)))>=realpow(2,qbit-3)+realpow(2,qbit-4)
        if maxmcv>3%3-1.8-944-7-5            
            
            input_block=sign(input_block).*floor(sign(input_block).*input_block/2);
            mcv=sign(mcv).*floor(sign(mcv).*mcv/2);
        
        end;
        
    end;
    
if it<=20
    decl=M;
else
    decl=4*L;
end;
    
for x=1:decl
    
        %if max(abs(mcv(:)))>=24%min(abs(input_block))>8%it==4 || it==6 || it==8 || it==13
            %fprintf('sub');
            %input_block=input_block-sign(input_block).*floor(sign(input_block).*input_block/4);
            %mcv=mcv-sign(mcv).*floor(sign(mcv).*mcv/4);
        
        %end;
    
%for x=M:-1:1
    for m=1:degree
        if (list(x,m)<K+1||(list(x,m)<N+K0-K+1&&list(x,m)>K0)) && list(x,m)~=0
            mvc(x,m)=input_block(list(x,m))-mcv(x,m);
            %mvc(x,m)=max(min(mvc(x,m),63),-63);
        else
            mvc(x,m)=maxm;
        end;
    end;
    
    %mcvt(x,:)=sumproduct(mvc(x,:),maxm);
    mcv(x,:)=minsum(mvc(x,:),max31);
    
    %mcv(x,:)=adjMS(mvc(x,:),maxm);

    if fpga==1
        peng=0;
        if peng==1
            for d=1:degree
                if mcv(x,d)>=0
                    if it<3
                        mcv(x,d)=mcv(x,d)-floor(mcv(x,d)/2);
                    elseif it<10
                        mcv(x,d)=mcv(x,d)-floor(mcv(x,d)/4);
                    else
                        mcv(x,d)=mcv(x,d)-floor(mcv(x,d)/8);
                    end;
                else
                    if it<3
                        mcv(x,d)=mcv(x,d)+floor(-mcv(x,d)/2);
                    elseif it<10
                        mcv(x,d)=mcv(x,d)+floor(-mcv(x,d)/4);
                    else
                        mcv(x,d)=mcv(x,d)+floor(-mcv(x,d)/8);
                    end;
                end;
            end;
        else
            for d=1:degree
                if mcv(x,d)>=0
                    if it<54
                        mcv(x,d)=mcv(x,d)-floor(mcv(x,d)/4);
                    elseif it<110
                        mcv(x,d)=mcv(x,d)-floor(mcv(x,d)/8);
                    else
                        mcv(x,d)=mcv(x,d)-floor(mcv(x,d)/16);
                    end;
                else
                    if it<54
                        mcv(x,d)=mcv(x,d)+floor(-mcv(x,d)/4);
                    elseif it<110
                        mcv(x,d)=mcv(x,d)+floor(-mcv(x,d)/8);
                    else
                        mcv(x,d)=mcv(x,d)+floor(-mcv(x,d)/16);
                    end;
                end;
            end;
        end;
    else
        %%%%%%%%%%%%%%%%%%%
       %if it<2
            %mcv(x,:)=mcv(x,:)*0.5;
       %elseif it<4
           mcv(x,:)=mcv(x,:)*0.75;
       %elseif it<6
           %mcv(x,:)=mcv(x,:)*0.85;
       %elseif it<9
           %mcv(x,:)=mcv(x,:)*0.95;
       %else
           %mcv(x,:)=mcv(x,:)*0.95;
       %end;
       %if it<4
            %mcv(x,:)=round(mcv(x,:)*0.75);
       %elseif it<10
           %mcv(x,:)=round(mcv(x,:)*0.875);
       %else
           %mcv(x,:)=round(mcv(x,:)*(1-1/16));
       %end;
       %%%%%%%%%%%%%%%%%%%%%
       %mcv(x,:)=mcv(x,:)*(0.75);
   end;

   for m=1:degree
        if list(x,m)<=N+K0-K && list(x,m)~=0
            input_block(list(x,m))=mcv(x,m)+mvc(x,m);
        end;
    end;
    if fpga==1
        for m=1:degree
            if list(x,m)<=N && list(x,m)~=0
                
                    if input_block(list(x,m))>maxm
                        input_block(list(x,m))=maxm;
                    elseif input_block(list(x,m))<-maxm
                        input_block(list(x,m))=-maxm;
                    end;
                
            end;
        end;
    end;
    debugp=0;
    if debugp==1
    if mod(x-1,L) == 7 && (x-8)/L == 0
        %if x==1 || x == 209 || x == 208*2+1 || x == 208*3 + 1 || x == 208*4 + 1 || x == 208*5 + 1 || x == 208*6 + 1 || x == 208*7 + 1 || x == 208*8 + 1 || x == 208*9 + 1 || x == 208*10 + 1 || x == 208*11 + 1 %M/22*4+1
        fprintf('mvc:');
        for m=1:degree
            fprintf('%d,',mvc(x,m));
        end;
        fprintf('\n');
        %fprintf('inp:');
        for m=1:degree
            if list(x,m)<=N && list(x,m)~=0
                %fprintf('%d,%d; ',input_block(list(x,m)),list(x,m));
            end;
        end;
        %fprintf('\n');
    end;
    end;
    
end;

%success
suc=1;
a=mvc<0;
for i=1:M
    if mod(sum(a(i,:)),2)~=0
        suc=0;
        if it==iteration
            %disp(i);
        end;
    end;
end;
if suc==1
    %fprintf('%d',it);
    
    %break;
end;
   
   %input_block=inputo;


end;

decode=input_block<0;
decode = decode(1:K);


function z = minsum(x,maxm)
degree=length(x);
signa=sign(1);
for i=1:degree
    if x(1,i)~=0
        signa=sign(x(1,i))*signa;
    end;
end;
sgn = zeros(1,degree);
for i=1:degree
    if sign(x(1,i))==0
        sgn(i)=signa;
    else
        sgn(i) = sign(x(1,i))*signa;
    end;
end;
zt=ones(1,degree)*maxm;
for i=1:degree
    for j=1:degree
        if i~=j
            zt(1,i)=min(abs(zt(1,i)),abs(x(1,j)));
        end;
    end;
end;
%zt(1,:)=zt(1,:)-floor(zt(1,:)/4);
z = zt(1,:) .* sgn;


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