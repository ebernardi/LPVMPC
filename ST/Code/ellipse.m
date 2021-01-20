function [xMax, yMax] =  ellipse(W, center, N, color, style)
% Plot the ellipse asociated to W matrix (x'Wx<=1)
%
% Inputs:
% W => W matrix (x'Wx<=1)
% center => Ellipse center
% N => Number of points
% color => Plot color line (only for 2D plot)
% style => Plot line style (only for 2D plot)
%
% Outputs:
% xMax = maximun value of x
% yMax = maximun value of y

    if nargin == 3  % Color and style by default
        color = 'blue';
        style = '-';
    end


    if length(center) == 3
        Type = '3D';
    elseif length(center) == 2
        Type = '2D';
    else
        disp('Cannot plot an ellipse with more than 3 dimensions!');
        return
    end

    [U, D, V] = svd(W);
    if strcmp(Type, '2D')
        % Get the major and minor axes
        a = 1/sqrt(D(1, 1));
        b = 1/sqrt(D(2, 2));
        theta = 0:1/N:2*pi+1/N;
        % Parametric equation of the ellipse
        state(1, :) = a*cos(theta); 
        state(2, :) = b*sin(theta);
        % Coordinate transform 
        X = V * state;
        X(1, :) = X(1,:) + center(1);
        X(2, :) = X(2,:) + center(2);

    elseif strcmp(Type, '3D')
        % Generate the ellipsoid at (0, 0, 0)
        a = 1/sqrt(D(1, 1));
        b = 1/sqrt(D(2, 2));
        c = 1/sqrt(D(3, 3));
        [X, Y, Z] = ellipsoid(0, 0, 0, a, b, c, N);

        % Rotate and center the ellipsoid to the actual center point
        XX = zeros(N+1, N+1);
        YY = zeros(N+1, N+1);
        ZZ = zeros(N+1, N+1);
        for k = 1:length(X)
            for j = 1:length(X)
                point = [X(k, j) Y(k, j) Z(k, j)]';
                P = V * point;
                XX(k, j) = P(1)+center(1);
                YY(k, j) = P(2)+center(2);
                ZZ(k, j) = P(3)+center(3);
            end
        end
    end

    % Plot the ellipse
    if strcmp(Type, '2D')
        plot(X(1, :), X(2, :), 'Color', color, 'LineWidth', 1, 'LineStyle', style);
        grid on
    else
        mesh(XX, YY, ZZ);
        hidden off
    end
    xMax = max(X(1, :));
    yMax = max(X(2, :));
end