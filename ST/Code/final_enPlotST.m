%% plot qLPV-MPC on ST
clc; clear; close all;
% load('run-qLPV-MPC.mat')
% load('run-linear.mat')
% load('run-MHE-MPC.mat')

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

% Perfonmance
for i = 1:3
    switch i
        case 1
            load('run-sMPC.mat')
            disp('<strong>sMPC</strong>')
        case 2
            load('run-qLPV-MPC.mat')
            disp('<strong>qLPV-MPC</strong>')
        case 3
            load('run-MHE-MPC.mat')
            disp('<strong>MHE-MPC</strong>')
            time_avg = mean(time_MHE+time_MPC);
            msg = ['Total mean percent = ', num2str(time_avg*100/Ts)];
            disp(msg)
            msg = ['Total mean time = ', num2str(time_avg)];
            disp(msg)
            time_avg = max(time_MHE+time_MPC);
            msg = ['Total max percent = ', num2str(time_avg*100/Ts)];
            disp(msg)
            msg = ['Total max time = ', num2str(time_avg)];
            disp(msg)
            time_avg = min(time_MHE+time_MPC);
            msg = ['Total min percent = ', num2str(time_avg*100/Ts)];
            disp(msg)
            msg = ['Total min time = ', num2str(time_avg)];
            disp(msg)
        otherwise
            break;
    end
    
    % Times
    time_avg = mean(time_MPC);
    msg = ['MPC mean percent = ', num2str(time_avg*100/Ts)];
    disp(msg)
    msg = ['MPC mean time = ', num2str(time_avg)];
    disp(msg)
    time_avg = max(time_MPC);
    msg = ['MPC max percent = ', num2str(time_avg*100/Ts)];
    disp(msg)
    msg = ['MPC max time = ', num2str(time_avg)];
    disp(msg)
    time_avg = min(time_MPC);
    msg = ['MPC min percent = ', num2str(time_avg*100/Ts)];
    disp(msg)
    msg = ['MPC min time = ', num2str(time_avg)];
    disp(msg)
    
    % Indices
    error = abs(X(2, :) - Xsp(2));
    IAE = trapz(t, abs(error));
    ISE = trapz(t, error.^2);
    ITAE = trapz(t, t.*abs(error));
    ITSE = trapz(t, t.*error.^2);
    msg = ['IAE = ', num2str(IAE)];
    disp(msg)
    msg = ['ITAE = ', num2str(ITAE)];
    disp(msg)
    msg = ['ISE = ', num2str(ISE)];
    disp(msg)
    msg = ['ITSE = ', num2str(ITSE)];
    disp(msg)
    
    input_shifted = [0; input(1:end-1)];
    msg = ['TV = ', num2str(sum(abs(input_shifted-input)))];
    disp(msg)

    error = abs(X(2, 1:5000) - Xsp(2));
    IAEtrack = trapz(t(1:5000), abs(error));
    ITAEtrack = trapz(t(1:5000), t(1:5000).*abs(error));
    error = abs(X(2, 5001:end) - Xsp(2));
    IAEpert = trapz(t(5001:end), abs(error));
    ITAEpert = trapz(t(5001:end), t(5001:end).*abs(error));

    msg = ['IAE tracking = ', num2str(IAEtrack)];
    disp(msg)
    msg = ['ITAE tracking = ', num2str(ITAEtrack)];
    disp(msg)
    msg = ['IAE disturbance = ', num2str(IAEpert)];
    disp(msg)
    msg = ['ITAE disturbance = ', num2str(ITAEpert)];
    disp(msg)
    disp(' ')
end

%% Outputs
fig = figure(3);
axes_p1 = axes('Parent', fig, 'Position', [0.13 0.72 0.775 0.2]);
plot(t, Xsp(1)*ones(1, length(t)), ':', 'Color', black, 'LineWidth', 1.5); hold(axes_p1, 'on');
load('run-sMPC.mat')
plot(t, X(1, :), '-', 'Color', orange_red, 'LineWidth', 1);
load('run-qLPV-MPC.mat')
plot(t, X(1, :), '--', 'Color', forest_green, 'LineWidth', 1);
load('run-MHE-MPC.mat')
plot(t, X(1, :), '-.', 'Color', royal_blue, 'LineWidth', 1); hold(axes_p1, 'off');
axis([0 Time 100 140]); grid on
xlabel('Time [s]'); ylabel('T_p [°C]');
leg = legend('set-point', 'sMPC', 'qLPV-MPC', 'MHE-MPC');
set(leg, 'Position', [0.728 0.741 0.164 0.156]);
leg.ItemTokenSize = [20, 10];
set(axes_p1, 'FontSize', 8)

axes_p2 = axes('Parent', fig, 'Position', [0.13 0.117 0.775 0.487]);
plot(t, Xsp(2)*ones(1, length(t)), ':', 'Color', black, 'LineWidth', 1.5); hold(axes_p2, 'on');
load('run-sMPC.mat')
plot(t, X(2, :), '-', 'Color', orange_red, 'LineWidth', 1);
load('run-qLPV-MPC.mat')
plot(t, X(2, :), '--', 'Color', forest_green, 'LineWidth', 1);
load('run-MHE-MPC.mat')
plot(t, X(2, :), '-.', 'Color', royal_blue, 'LineWidth', 1); hold(axes_p2, 'off');
axis([0 Time 95 97.5]); grid on
xlabel('Time [s]'); ylabel('T_f [°C]');
set(axes_p2, 'FontSize', 8)

% Create new axes
axes_1 = axes('Parent', fig, 'Position', [0.407 0.154 0.47 0.2912]);
plot(t, Xsp(2)*ones(1, length(t)), ':', 'Color', black, 'LineWidth', 1.5); hold(axes_1, 'on');
load('run-sMPC.mat')
plot(t, X(2, :), '-', 'Color', orange_red, 'LineWidth', 1);
load('run-qLPV-MPC.mat')
plot(t, X(2, :), '--', 'Color', forest_green, 'LineWidth', 1);
load('run-MHE-MPC.mat')
plot(t, X(2, :), '-.', 'Color', royal_blue, 'LineWidth', 1); hold(axes_1, 'off');
box(axes_1, 'on'); grid(axes_1, 'on');
set(axes_1, 'FontSize', 6)
xlim(axes_1, [500 2500]); ylim(axes_1, [96.96 97.01]);

% Create new axes
axes_2 = axes('Parent', fig, 'Position', [0.198 0.154 0.147 0.2912]);
plot(t, Xsp(2)*ones(1, length(t)), ':', 'Color', black, 'LineWidth', 1.5); hold(axes_2, 'on');
load('run-sMPC.mat')
plot(t, X(2, :), '-', 'Color', orange_red, 'LineWidth', 1);
load('run-qLPV-MPC.mat')
plot(t, X(2, :), '--', 'Color', forest_green, 'LineWidth', 1);
load('run-MHE-MPC.mat')
plot(t, X(2, :), '-.', 'Color', royal_blue, 'LineWidth', 1); hold(axes_2, 'off');
box(axes_2, 'on'); grid(axes_2, 'on');
set(axes_2, 'FontSize', 6)
xlim(axes_2, [0 100]); ylim(axes_2, [96 97.01]);

% Create annotations
annotation(fig, 'arrow', [0.569 0.607], [0.823 0.771], 'HeadWidth', 6, 'HeadLength', 6);
annotation(fig, 'textarrow', [0.375 0.339], [0.82 0.771], 'FontSize', 8, 'HeadWidth', 6, ...
                    'HeadLength', 6, 'String', {'Disturbance effect'}, 'HorizontalAlignment', 'center');
annotation(fig, 'doublearrow', [0.585 0.585], [0.37 0.223], 'Head2Width', 6, ...
                    'Head2Length', 6, 'Head1Width', 6, 'Head1Length', 6);
annotation(fig, 'textbox', [0.578 0.285 0.223 0.039], 'String', {'sMPC tracking error'}, ...
                    'FontSize', 8, 'FitBoxToText', 'off', 'EdgeColor', 'none');
            
print -dsvg ../Figs/state.svg

%% Input
fig4 = figure(4);
axes1 = axes('Parent', fig4, 'Position', [0.121 0.383 0.775 0.325]);
load('run-sMPC.mat')
stairs(Tsim, input, '-', 'Color', orange_red, 'LineWidth', 1); hold(axes1, 'on');
load('run-qLPV-MPC.mat')
stairs(Tsim, input, '--', 'Color', forest_green, 'LineWidth', 1);
load('run-MHE-MPC.mat')
stairs(Tsim, input, '-.', 'Color', royal_blue, 'LineWidth', 1);
xlabel('Time [s]'); ylabel('u [m^3/s]');
xlim(axes1, [0 Time]); ylim(axes1, [0 1e-3]); grid(axes1, 'on');
set(axes1, 'FontSize', 8);
leg = legend('sMPC', 'qLPV-MPC', 'MHE-MPC');
set(leg, 'Location', 'NorthEast'); %'Position', [0.728 0.741 0.164 0.156]);
leg.ItemTokenSize = [20, 10];

% Create axes
axes2 = axes('Parent', fig4, 'Position', [0.2526 0.5047 0.4071 0.1738]);
load('run-sMPC.mat')
stairs(Tsim, input, '-', 'Color', orange_red, 'LineWidth', 1); hold(axes2, 'on');
load('run-qLPV-MPC.mat')
stairs(Tsim, input, '--', 'Color', forest_green, 'LineWidth', 1);
load('run-MHE-MPC.mat')
stairs(Tsim, input, '-.', 'Color', royal_blue, 'LineWidth', 1);
xlim(axes2, [700 1200]); ylim(axes2, [1.1e-4 2e-4]); grid(axes2, 'on');
set(axes2, 'FontSize', 6);

print -dsvg ../Figs/input.svg

%% Objective function
figure(8)
hold on
plot(Tsim, Obj(:))
xlim([0 Time])
xlabel('Muestra'); ylabel('objective');

%% Disturbance
figure(7)
sub1 = subplot(211, 'FontSize', 4);
plot(Tsim, W(1, :), '-.', 'Color', orange, 'LineWidth', 1.5);
axis([0 Time 0 805])
xlabel('Time [s]'); ylabel('I [W/m^2]'); grid on
set(sub1, 'FontSize', 8, 'Position', [0.13 0.583 0.775 0.266]);
sub2 = subplot(212);
plot(Tsim, W(2, :), '-', 'Color', dark_blue, 'LineWidth', 1.5);
axis([0 Time 27.9 28.8])
xlabel('Time [s]'); ylabel('T_e [°C]'); grid on
set(sub2, 'FontSize', 8, 'Position', [0.13 0.192 0.775 0.258]);

print -dsvg ../Figs/disturbance.svg

%% Membership
fig5 = figure(5);
axes1 = axes('Parent', fig5, 'Position', [0.121 0.311 0.775 0.447]);
plot(Tsim, mu_fuzzy(1, :), '-', 'Color', lightblue); hold on; grid on;
plot(Tsim, mu_mhe(1, :), '-.', 'Color', blue);
plot(Tsim, mu_fuzzy(2, :), '-', 'Color', forest_green);
plot(Tsim, mu_mhe(2, :), '-.', 'Color', green);
plot(Tsim, mu_fuzzy(3, :), '-', 'Color', violet);
plot(Tsim, mu_mhe(3, :), '-.', 'Color', red);
plot(Tsim, mu_fuzzy(4, :), '-', 'Color', chocolate);
plot(Tsim, mu_mhe(4, :), '-.', 'Color', gold);
xlim([0 Time]); hold off;
xlabel('Time [s]'); ylabel('\mu_i'); grid on
leg = legend('Real \mu_1', 'MHE \mu_1', 'Real \mu_2', 'MHE \mu_2', 'Real \mu_3', 'MHE \mu_3', 'Real \mu_4', 'MHE \mu_4');
set(leg, 'Orientation', 'horizontal', 'Location', 'North', 'FontSize', 6);
leg.ItemTokenSize = [12, 10];
set(axes1, 'FontSize', 8);
print -dsvg ../Figs/membership.svg

%% Period consumption
fig = figure(10);
axes1 = axes('Parent', fig, 'Position', [0.13 0.314 0.775 0.442]);
area(Tsim, (time_MHE+time_MPC)*100/Ts, 'FaceColor', green); hold(axes1, 'on');
area(Tsim, time_MPC*100/Ts, 'FaceColor', yellow, 'FaceAlpha', 0.5);
area(Tsim, time_MHE*100/Ts, 'FaceColor', blue, 'FaceAlpha', 1);
xlabel('Time [s]'); ylabel('Period [%]'); grid on
axis([0 Time 0 2]); grid(axes1, 'on');
leg = legend('MHE+MPC', 'MPC', 'MHE');
set(leg, 'Location', 'NorthEast');
set(axes1, 'FontSize', 8);
leg.ItemTokenSize = [20, 10];

print -dsvg ../Figs/periodLoad.svg

%% State space
load('run-MHE-MPC.mat')
figure(6)
ellipse(Wbmi, Xsp, 20, 'black', '-')
hold on
plot(X(1, :), X(2, :), '-g')
plot(Xsp(1), Xsp(2), '*', 'Color', purple, 'LineWidth', 1.5);
xlabel('T_p [K]'), ylabel('T_f [K]')

% print -dsvg ../Figs/space_state.svg