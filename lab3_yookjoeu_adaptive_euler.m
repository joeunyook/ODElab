
function [t, y] = lab3_yookjoeu_adaptive_euler(f, t0, tN, y0, h)
    % Initialize
    t = t0;
    y = y0;
    tol = 1e-8;  % tolerance for the error
    
    % Helper function to compute Y, Z, and D 
    function [Y, Z, D] = helper(y, h, t)
        Y = y + h * f(t, y);
        Z1 = y + 0.5 * h * f(t, y);  % First half-step
        Z = Z1 + 0.5 * h * f(t + 0.5 * h, Z1);  % Second half-step

        %  difference D for estimate of the error
        D = Z - Y;
    end
    
    % Main loop: Perform the Adaptive Euler method
    while t(end) < tN
        [Y, Z, D] = helper(y(end), h, t(end));
        
        while abs(D) >= tol
            % If the error is too large, update the step size
            h = 0.9 * h * min(max(tol / abs(D), 0.3), 2);
            [Y, Z, D] = helper(y(end), h, t(end)); 
        end
        
        % Once the error is acceptable, update the solution
        y = [y, Z];  % Use the more accurate Z value
        t = [t, t(end) + h];  % Move to the next time step
    end
end


