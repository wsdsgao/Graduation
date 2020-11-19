function [Z,matrix,bg]=bg2()

Z=10;

ils=4;

bg=zeros(42,52);
matrix=zeros(42,52)-1;
data=importdata('bg2.txt');
row=1;
last=0;
for i=1:length(data(:,1))
    col=data(i,1);
    n=data(i,1+ils);
    if col<last
        row=row+1;
    end;
    last=col;
    bg(row,col+1)=1;
    matrix(row,col+1)=n;
end;

show=1;
if show==1
    matix=matrix;
for i=1:42
    for j=1:52
        if matrix(i,j)<0
            matix(i,j)=0;
        else
            matix(i,j)=matrix(i,j)+100;
        end;
    end;
end;

ma=1/max(max(matix));
% imshow(1-ma*(matix+1));
end;