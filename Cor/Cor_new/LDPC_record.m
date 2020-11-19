close all;
clear all;

Eb_N0 = -1.6:0.2:-0.8;
% Eb_N0 = [Eb_N0, -0.8, -0.4, 0, 0.2, 0.6, 1];
coherent_ber_record = [0.00971679687500000,0.00605468750000000,0.00578125000000000,0.00238281250000000,0]
coherent_bulk_record = [7 6 8 4 0];
coherent_bulk_record = coherent_bulk_record/100;

f = figure;
f.PaperUnits = 'centimeters';
f.PaperSize = [16, 12];
f.Units = 'centimeters';
f.Position = [0, 0, 16, 12];
semilogy(Eb_N0, coherent_ber_record, '-s', 'LineWidth', 2);
grid on;
set(get(gca, 'XLabel'), 'String', 'Es/N0');

figure;
semilogy(Eb_N0, coherent_bulk_record, '-s', 'LineWidth', 2);
grid on;

% [-1.60000000000000,-1.40000000000000,-1.20000000000000,-1,-0.8]

% [0.00971679687500000,0.00605468750000000,0.00578125000000000,0.00238281250000000,0]
% 7 6 8 4