function [time, x] = solvesystem_yookjoeu(f, g, t0, tN, x0, h)
    len = round((tN - t0) / h, 0);
    time = linspace(t0, tN, len);
    x = zeros(2, len);
    
    % Initial values:
    time(1) = t0;
    x(:,1) = x0; % x(:,1) is the first column 
    
    for i = 2:len
        time(i) = time(i - 1) + (i - 1) * h;
        
        % First step
        x(1,i) = x(1, i - 1) + h * f(time(i - 1), x(1, i - 1), x(2, i - 1));
        x(2,i) = x(2, i - 1) + h * g(time(i - 1), x(1, i - 1), x(2, i - 1));
        
        % Slope of left
        S_L = [f(time(i - 1), x(1, i - 1), x(2, i - 1));
               g(time(i - 1), x(1, i - 1), x(2, i - 1))];
        
        % Slope of right
        S_R = [f(time(i), x(1, i), x(2, i));
               g(time(i), x(1, i), x(2, i))];
        
        % Improved slope
        improv_slope = (S_L + S_R) ./ 2;
        x(:,i) = x(:, i - 1) + h * improv_slope;
    end 
end
