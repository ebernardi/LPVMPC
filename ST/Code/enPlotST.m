%% plot qLPV-MPC on ST
clc; clear; close all;
load('run.mat')

red = [0.7; 0; 0]; 
black = [.1; .1; .1]; 
gray = [.5; .7; .5];
blue = [0 0.4470 0.7410];
orange = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];
violet = [0.4940 0.1840 0.5560];
green = [0.4660 0.6740 0.1880];
lightblue = [0.3010 0.7450 0.9330];
purple = [0.6350 0.0780 0.1840];

orange_red = [255 69 0]/255;
forest_green = [34 139 34]/255;
royal_blue = [65 105 225]/255;
dark_blue = [0 0 139]/255;
gold = [255 215 0]/255;
chocolate = [210 105 30]/255;
arrow = [212 55 144]/255;

error = abs(X(2, :) - Xsp(2));
IAE = sum(error)/length(error);
msg = ['IAE = ', num2str(IAE)];
disp(msg)

time_avg = mean(time_MHE+time_MPC) ;
msg = ['Mean time = ', num2str(time_avg)];
disp(msg)

% Perfonmance indices
error = abs(X(2, :) - Xsp(2));
IAE = trapz(t, abs(error));
ISE = trapz(t, error.^2);
ITAE = trapz(t, t.*abs(error));
msg = ['IAE = ', num2str(IAE)];
disp(msg)
msg = ['TV = ', num2str(sum(deltaU(1, 2:end)))];
disp(msg)
msg = ['ITAE = ', num2str(ITAE)];
disp(msg)
msg = ['ISE = ', num2str(ISE)];
disp(msg)

error = abs(X(2, 1:7000) - Xsp(2));
IAEtrack = sum(error)/length(error);
error = abs(X(2, 7001:end) - Xsp(2));
IAEpert = sum(error)/length(error);

msg = ['IAE tracking = ', num2str(IAEtrack)];
disp(msg)
msg = ['ITAE disturbance = ', num2str(IAEpert)];
disp(msg)

time_avg = mean(time_MHE) ;
msg = ['Mean time = ', num2str(time_avg)];
disp(msg)
time_avg = max(time_MHE) ;
msg = ['Max time = ', num2str(time_avg)];
disp(msg)
time_avg = min(time_MHE) ;
msg = ['Min time = ', num2str(time_avg)];
disp(msg)

%% Outputs
fig = figure(3);
subplot(211)
plot(t, X(1, :), 'Color', chocolate, 'LineWidth', 1.5);
axis([0 Time 100 300]); grid on
xlabel('Time [s]'); ylabel('T_p [°C]');
leg = legend('x_1(t)');
set(leg, 'Location', 'NorthEast');
leg.ItemTokenSize = [10, 10];
subplot(212)
plot(t, Xsp(2)*ones(1, length(t)), ':', 'Color', orange_red, 'LineWidth', 1.5);
hold on
plot(t, X(2, :), '-', 'Color', blue, 'LineWidth', 1.5);
hold off; axis([0 Time 92 98]); grid on
xlabel('Time [s]'); ylabel('T_f [°C]');
leg = legend('x_2sp', 'x_2(t)');
set(leg, 'Location', 'East');
leg.ItemTokenSize = [10, 10];
% Create new axes
axes_1 = axes('Parent', fig, 'Position', [0.307 0.167 0.426 0.212], 'FontSize', 6);
hold(axes_1, 'on');
plot(t, Xsp(2)*ones(1, length(t)), ':', 'Color', orange_red, 'LineWidth', 1.5); hold on; grid on
plot(t, X(2, :), '-', 'Color', blue, 'LineWidth', 1.5); hold off;
box(axes_1, 'on'); grid(axes_1, 'on');
xlim(axes_1, [0 80]); ylim(axes_1, [92 98]);
% Create arrow
annotation(fig, 'arrow', [0.629 0.707], [0.726 0.63]);
annotation(fig, 'textarrow', [0.491 0.421], [0.723 0.635], ...
                    'String', {'Disturbance', 'effect'}, 'HorizontalAlignment', 'center');
print -dsvg ../Figs/state.svg

%% Input
figure(4)
stairs(Tsim, input, 'Color', orange_red, 'LineWidth', 1.5);
xlabel('Time [s]'); ylabel('u [m^3/s]'); grid on
xlim([0 Time])
pbaspect([2 1 1]);
print -dsvg ../Figs/input.svg

%% Objective function
figure(6)
plot(Tsim, Obj(:))
xlim([0 Time])
xlabel('Muestra'); ylabel('objective');

%% Disturbance
figure(7)
subplot(211)
plot(Tsim, W(1, :), '-.', 'Color', blue, 'LineWidth', 1.5);
xlim([0 Time])
xlabel('Time [s]'); ylabel('I [W/m^2]'); grid on
subplot(212)
plot(Tsim, W(2, :), 'g-', 'LineWidth', 1.5);
xlim([0 Time])
xlabel('Time [s]'); ylabel('T_e [°C]'); grid on
print -dsvg ../Figs/disturbance.svg

figure(5);
hold on
plot(Tsim, mu_fuzzy(1, :)); hold on; grid on;
plot(Tsim, mu_mhe(1,:));
plot(Tsim, mu_fuzzy(2, :));
plot(Tsim, mu_mhe(2,:));
plot(Tsim, mu_fuzzy(3, :));
plot(Tsim, mu_mhe(3,:));
plot(Tsim, mu_fuzzy(4, :));
plot(Tsim, mu_mhe(4,:));
xlim([0 Time]); hold off; grid on
legend('real1', 'MHE 1', 'real2', 'MHE 2', 'real3', 'MHE 3', 'real4', 'MHE 4')

figure(10);
area((time_MHE+time_MPC)*100/3, 'FaceColor', green);
hold on
area(time_MPC*100/Ts, 'FaceColor', yellow, 'FaceAlpha', 0.5);
area(time_MHE*100/Ts, 'FaceColor', blue, 'FaceAlpha', 1);
xlabel('Iteration [k]'); ylabel('Computation effort [%]'); grid on
axis([0 Nsim 0 3]); grid on
pbaspect([2 1 1]);

%% State space
figure(6)
% plot3(x0(1), x0(2), x0(3), 'o', 'Color', verde, 'LineWidth', 1.5);
% plot3(Xsp(1, :), Xsp(2, :), Xsp(3, :), '*', 'Color', bordo, 'LineWidth', 1.5);
% plot3(X(1, :), X(2, :), X(3, :), '-.', 'Color', azul, 'LineWidth', 1.5);
% plot3(X(1, :), X(2, :), X(3, :), 'Color', chocolate, 'LineWidth', 1.5);
% plot(X, 'color', 'b', 'alpha', 0.05, 'edgecolor',  blue)
ellipse(Wbmi, Xsp, 20, 'black', '-')
hold on
plot(X(1, :), X(2, :), '-g')
plot(Xsp(1), Xsp(2), '*', 'Color', purple, 'LineWidth', 1.5);
xlabel('T_p [K]'), ylabel('T_f [K]')
% leg = legend('x_0', 'x_s', 'nominal-MPC', 'LPV-MPC', 'X_f_{LPV}', 'Location', 'Northeast');
% leg.ItemTokenSize = [20, 15];

% print -dsvg ../Figs/space_state.svg

% %% Contractive Terminal Sets
% fig = figure(8); hold on; box on
% xlabel('T_p [°C]'); ylabel('T_f [°C]'); grid on
% plot(X(1, :), X(2, :), 'b-.', 'LineWidth', 1.5);
% plot(xsp(1), xsp(2), '*', 'Color', bordo, 'LineWidth', 1.5);
% plot(X(1, 1), X(2, 1), 'o', 'Color', forest_green, 'LineWidth', 1.5); 
% for j = 1:Nr
%     vect_Color = [j/Nr; 1-j/Nr; 1-j/Nr];
%     select_Upsilon = max([Nr-j+1 1]);
%     plot(Upsilon_j(select_Upsilon), 'Color', vect_Color, 'Alpha', 0.01, 'edgecolor', vect_Color);
% end
% plot(X(1, :), X(2, :), 'b-.', 'LineWidth', 1.5);
% plot(xsp(1), xsp(2), '*', 'Color', bordo, 'LineWidth', 1.5);
% plot(X(1, 1), X(2, 1), 'o', 'Color', forest_green, 'LineWidth', 1.5); 
% axis([70 620 90 135]); hold off
% leg = legend('x(t)', 'p_t', 'x_0');
% set(leg, 'Position', [0.786 0.642 0.075 0.147]);
% leg.ItemTokenSize = [10, 10];
% legend boxoff
% annotation(fig, 'rectangle', [0.786 0.642 0.075 0.147]);
% % Create new axes
% axes_1 = axes('Parent', fig, 'Position', [0.394 0.337 0.306 0.291], 'FontSize', 6);
% hold(axes_1, 'on');
% for j = 1:Nr
%     vect_Color = [j/Nr; 1-j/Nr; 1-j/Nr];
%     select_Upsilon = max([Nr-j+1 1]);
%     plot(Upsilon_j(select_Upsilon), 'Color', vect_Color, 'Alpha', 0.01, 'edgecolor', vect_Color);
% end     
% plot(X(1, :), X(2, :), 'b-.', 'LineWidth', 1.5);
% plot(xsp(1), xsp(2), '*', 'Color', bordo, 'LineWidth', 1.5);
% plot(X(1, 1), X(2, 1), 'o', 'Color', forest_green, 'LineWidth', 1.5);
% box(axes_1, 'on'); grid(axes_1, 'on');
% xlim(axes_1, [80 160]); ylim(axes_1, [95 99]);
% 
% % Create arrow
% annotation(fig, 'arrow', [0.444 0.434], [0.878 0.799], ...
%     'HeadWidth', 6, 'HeadStyle', 'vback1', 'HeadLength', 6);
% annotation(fig, 'textbox', [0.445 0.821 0.135 0.0799], ...
%     'String', {'Contractive', 'Sets'}, 'LineStyle', 'none', ...
%     'HorizontalAlignment', 'center', 'FontSize', 8, 'FitBoxToText', 'off');
% annotation(fig,'textarrow', [0.214 0.193], [0.307 0.259], 'String', {'Target', 'State'}, 'FontSize', 8, ...
%     'HorizontalAlignment', 'center', 'HeadWidth', 6, 'HeadStyle', 'vback1', 'HeadLength', 6);
% annotation(fig, 'textarrow', [0.478 0.44], [0.181 0.156], 'String', {'Initial State'}, 'FontSize', 8, ...
%     'HorizontalAlignment','center', 'HeadWidth', 6, 'HeadStyle', 'vback1', 'HeadLength', 6);
% print -dsvg sets.svg