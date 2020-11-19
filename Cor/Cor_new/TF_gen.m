function [th_pat, fh_pat] = TF_gen
% ������96�������Ӧ����Ƶ\��ʱͼ��������48��Ϊһ���ظ�2��

% 2Mbps A ȡͼ����ǰ12��
% 2Mbps B ȡͼ����ǰ12��
% 500K ȡͼ����ǰ48��
% 250K ȡ����ͼ��

% ��ʱ (48������)
for i = 1:4
    th_pat_part1(i,:) = 2*randi([15, 240], [1, 6]);  
    th_pat_part2(i,:) = 512 - th_pat_part1(i,:);  
    th_pat_temp((i-1)*12+1:i*12) = [th_pat_part1(i,:), th_pat_part2(i,:)];
end

% ��Ƶ�� Ĭ��21��Ƶ��ȫ������
% ���˴��ɸı����Ƶ��������
for i = 1:3
    temp = randperm(21);
    fh_pat_temp1((i-1)*21+1:i*21) = temp;
end
fh_pat_temp = fh_pat_temp1(1:48);

th_pat = repmat(th_pat_temp, [1,2]);
fh_pat = repmat(fh_pat_temp, [1,2]);
