clear all;
close all;

open('./fig/res_lpf_7_9_deep_30_win.fig')
lh = findall(gca,'type','line');
bit1 = get(lh, 'ydata');
x = get(lh, 'xdata');
close all;

open('./fig/res2_lpf_7_9_deep_40_win.fig')
lh = findall(gca,'type','line');
bit2 = get(lh, 'ydata');
close all;

open('./fig/res3_lpf_7_9_deep_80_win.fig')
lh = findall(gca,'type','line');
bit3 = get(lh, 'ydata');
close all;

open('./fig/res_lpf_7_9_deep_30_win_freq.fig')
lh = findall(gca,'type','line');
bit1_fre = get(lh, 'ydata');
close all;

open('./fig/res2_lpf_7_9_deep_40_win_freq2.fig')
lh = findall(gca,'type','line');
bit2_fre = get(lh, 'ydata');
close all;

open('./fig/res3_lpf_7_9_deep_80_win_freq.fig')
lh = findall(gca,'type','line');
bit3_fre = get(lh, 'ydata');
close all;

f = figure;
f.PaperUnits = 'centimeters';
f.PaperSize = [16, 12];
f.Units = 'centimeters';
f.Position = [0, 0, 16, 12];
semilogy(x, bit1, '-rs', 'LineWidth', 2);
hold on;
semilogy(x, bit2, '-gs', 'LineWidth', 2);
hold on;
semilogy(x, bit3, '-bs', 'LineWidth', 2);
hold on;
semilogy(x, bit1_fre, '-s', 'LineWidth', 2);
hold on;
semilogy(x, bit2_fre, '-s', 'LineWidth', 2);
hold on;
semilogy(x, bit3_fre, '-s', 'LineWidth', 2);
xlim([min(x)-1, max(x)+1]);
grid on;
legend('1bit','2bit','3bit','1bit\_fre','2bit\_fre','3bit\_fre');
print(f,'-dpng','-r300','./image/BER_fre.png');

f = figure;
f.PaperUnits = 'centimeters';
f.PaperSize = [16, 12];
f.Units = 'centimeters';
f.Position = [0, 0, 16, 12];
semilogy(x, bit1, '-rs', 'LineWidth', 2);
hold on;
semilogy(x, bit2, '-gs', 'LineWidth', 2);
hold on;
semilogy(x, bit3, '-bs', 'LineWidth', 2);
xlim([min(x)-1, max(x)+1]);
grid on;
legend('1bit','2bit','3bit');
print(f,'-dpng','-r300','./image/BER.png');
