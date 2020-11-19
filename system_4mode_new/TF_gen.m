function [th_pat, fh_pat] = TF_gen
% 共生成96个脉冲对应的跳频\跳时图案，是以48个为一组重复2遍

% 2Mbps A 取图案的前12个
% 2Mbps B 取图案的前12个
% 500K 取图案的前48个
% 250K 取完整图案

% 跳时 (48个脉冲)
for i = 1:4
    th_pat_part1(i,:) = 2*randi([15, 240], [1, 6]);  
    th_pat_part2(i,:) = 512 - th_pat_part1(i,:);  
    th_pat_temp((i-1)*12+1:i*12) = [th_pat_part1(i,:), th_pat_part2(i,:)];
end

% 跳频， 默认21个频点全部可用
% （此处可改变可用频点总数）
for i = 1:3
    temp = randperm(21);
    fh_pat_temp1((i-1)*21+1:i*21) = temp;
end
fh_pat_temp = fh_pat_temp1(1:48);

th_pat = repmat(th_pat_temp, [1,2]);
fh_pat = repmat(fh_pat_temp, [1,2]);
