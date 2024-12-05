clc; 
clear all;
close all;

x = 0 : 0.001 : 1;
n = length(x);

y = 1 * (x <= 0.5) + 0.846 * (0.5 < x & x <= 0.55) + 0.444 * (0.55 < x & x <= 0.6) + ...
    0 * (x > 0.6);

plot(x, y, 'r', 'LineWidth', 3);
xlabel('Recall');
ylabel('Precision');
%title('P-R curve');
legend('P-R curve');