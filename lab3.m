%% ODE Lab: Creating your own ODE solver in MATLAB
% In this lab, you will write your own ODE solver for the Improved Euler method 
% (also known as the Heun method), and compare its results to those of |ode45|.
% 
% You will also learn how to write a function in a separate m-file and execute 
% it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each part using 
% cell mode to see the results. Compare the output with the PDF, which was generated 
% from this m-file.
% 
% There are six (6) exercises in this lab that are to be handed in on the due 
% date. Write your solutions in the template, including appropriate descriptions 
% in each step. Save the .m files and submit them online on Quercus.
%% Student Information
% Student Name: Joeun Yook
% 
% Student Number: 1010101462
%% Creating new functions using m-files.
% Create a new function in a separate m-file:
% 
% Specifics: Create a text file with the file name f.m with the following lines 
% of code (text):
%%
% 
%  function y = f(a,b,c)
%  y = a+b+c;
%
%% 
% Now MATLAB can call the new function f (which simply accepts 3 numbers and 
% adds them together). To see how this works, type the following in the matlab 
% command window: sum = f(1,2,3)
%% Exercise 1
% Objective: Write your own ODE solver (using the Heun/Improved Euler Method).
% 
% Details: This m-file should be a function which accepts as variables (t0,tN,y0,h), 
% where t0 and tN are the start and end points of the interval on which to solve 
% the ODE, y0 is the initial condition of the ODE, and h is the stepsize. You 
% may also want to pass the function into the ODE the way |ode45| does (check 
% lab 2).
% 
% Note: you will need to use a loop to do this exercise. You will also need 
% to recall the Heun/Improved Euler algorithm learned in lectures.
% 
% 

% solution written in lab3_yookjoeu_improved_euler.m file
%% Exercise 2
% Objective: Compare Heun with |ode45|.
% 
% Specifics: For the following initial-value problems (from lab 2, exercises 
% 1, 4-6), approximate the solutions with your function from exercise 1 (Improved 
% Euler Method). Plot the graphs of your Improved Euler Approximation with the 
% |ode45| approximation.
% 
% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|
% 
% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|
% 
% (c) |y' = 1 - t y / 2, y(0) = -1| from |t=0| to |t=10|
% 
% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|
% 
% Comment on any major differences, or the lack thereof. You do not need to 
% reproduce all the code here. Simply make note of any differences for each of 
% the four IVPs.
% 
% 

% (a)
% for ODE45
f_a = @(t,y) y.*tan(t) + sin(t);
% initial condition
t0 = 0;
tN = pi;
y0 = -1/2;
sol_by_45 = ode45(f_a, [t0, tN], y0);

% ODE solver (Improved Euler)
[Imp_x, Imp_y] = lab3_yookjoeu_improved_euler(f_a, t0,tN,y0,0.001);
subplot(2,2,1);
plot(sol_by_45.x, sol_by_45.y(1,:), 'b', Imp_x, Imp_y, 'r'); 
legend('ODE45', 'My ODE', 'Location', 'Best');
title("a) y' = y*(tan(t))+(sin(t))");
ylabel('y');
% No recognizable difference 


% (b)
f_b = @(t,y) 1/(y^2);
t0 = 1;
y0 = 1;
tN = 10;
sol_45 = ode45(f_b, [t0, tN], y0);
[Imp_x, Imp_y] = lab3_yookjoeu_improved_euler(f_b, t0,tN,y0,0.001);
subplot(2,2,2);
plot(sol_45.x, sol_45.y(1,:), 'b', Imp_x, Imp_y, 'r'); 
legend('ODE 45', 'My ODE','Location','Best');
title("b) y' = 1 / y^2 ");
xlabel('t');
ylabel('y');
% No recognizable difference 

% (c)
f_c = @(t,y) 1 - t*y/2;
t0 = 0;
y0 = -1;
tN = 10;
sol45 = ode45(f_c, [t0, tN], y0);
[Imp_x, Imp_y] = lab3_yookjoeu_improved_euler(f_c, t0,tN,y0,0.001);
subplot(2,2,3);
plot(sol45.x, sol45.y(1,:), 'b', Imp_x, Imp_y, 'r'); 
xlabel('t');
ylabel('y');
title("c) y' = 1 - t y / 2");
legend("ODE45", "My ODE", "Location", "Best");
%My ode has smaller step size than ode45, resulting in smoother curve than ode45


% (d)
f_d = @(t,y) y.^3 - t.^2;
t0 = 0;
y0 = 1;
tN = 1;
sol45 = ode45(f_d, [t0, tN], y0);
[Imp_x, Imp_y] = lab3_yookjoeu_improved_euler(f_d, t0,tN,y0,0.001);
subplot(2,2,4);
plot(sol45.x, sol45.y(1,:), 'b', Imp_x, Imp_y, 'r'); 
xlabel('t');
ylabel('y');
title("d) y' = y^3 - t^2");
legend("ODE45", "My ODE", "Location", "Best");

% ode45 does not show the plot(it goes to infinity). Instead, results in the Warning message
% saying : Failure at t=5.066046e-01. Unable to meet integration tolerances without 
% reducing the step size below the smallest value allowed (1.776357e-15) at time t.
% At appraximately 0.5, My ode plot soars to infinity to the factor of 10.^22


%% Exercise 3
% Objective: Use Euler's method and verify an estimate for the global error.
% 
% Details:
% 
% (a) Use Euler's method (you can use euler.m from iode) to solve the IVP
% 
% |y' = 2 t sqrt( 1 - y^2 ) , y(0) = 0|
% 
% from |t=0| to |t=0.5|.
% 
% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.
% 
% (c) Read the attached derivation of an estimate of the global error for Euler's 
% method. Type out the resulting bound for En here in a comment. Define each variable.
% 
% (d) Compute the error estimate for |t=0.5| and compare with the actual error.
% 
% (e) Change the time step and compare the new error estimate with the actual 
% error. Comment on how it confirms the order of Euler's method.

% (a) Use Euler's method to solve the IVP y' = 2 t sqrt(1 - y^2), y(0) = 0 from t=0 to t=0.5
f = @(t,y) 2.*t.*sqrt(1 - y.^2);
dt = 0.01;
t = 0:dt:0.5;
xc = euler(f, 0, t); % Euler's method

% Result at t=0.5
fprintf('Eulers method result at t = 0.5: %f\n', xc(end));




% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.

% Solving ode gives y = sin(t^2 + C)
% Then, C = 0 given initial value is (0,0) 
exact_sol = @(t) sin(t.^2);
fprintf('Solution of IVP at t = 0.5 y(0.5): %f. \n', sin(0.5^2));


% (c) Read the attached derivation of an estimate of the global error for Euler's 
% method. Type out the resulting bound for En here in a comment. Define each variable.

% E_n = (1+M)/2 * dt * M * dt * n
% E_n = error at n
% M = [t0, tN]: there exists an M > 0 s.t |f | ≤ M , |∂tf | ≤ M , and |∂yf | ≤ M .
% dt = h = time = step size = 0.01 here
% n = number of steps, = 100 here
% f = 2 t sqrt( 1 - y^2 )
% ∂tf = 2*sqrt( 1 - y^2 )
% ∂yf = t*(1-y^2)^(-1/2)*(-2y)

t = 0:0.01:0.5;
f = 2.*t.*sqrt(1-xc.^2);
d_tf = 2.*sqrt(1-xc.^2);
d_yf = t.*(1-xc.^2).^(-1/2).*(-2*xc);

max(abs(f));
max(abs(d_tf));
max(abs(d_yf));

M = max(abs(d_tf)) ;
dt = 0.01;
n = 100;
E_n = (1+M)/2 * dt * M * dt * n 
% resulting bound of error is  ±0.0300
%% 
% (d) Compute the error estimate for |t=0.5| and compare with the actual error.

actual_error = abs(xc(end) - exact_sol(0.5));
fprintf('With a step size of 0.01 and n = 100, the error bound is %f.\n', E_n);
fprintf('Actual error at t = 0.5: %f\n', actual_error);
error_difference = abs(actual_error - E_n);
fprintf('Difference between estimated and actual error: %f\n', error_difference);


% (e) Change the time step and compare the new error estimate with the actual 
% error. Comment on how it confirms the order of Euler's method.

xc=euler(@(t, y) 2*t*sqrt(1-y.^2), 0, 0:0.001:0.5);
t = 0:0.001:0.5;
f = 2.*t.*sqrt(1-xc.^2);
dtf = 2.*sqrt(1-xc.^2);
dyf = t.*(1-xc.^2).^(-1/2).*(-2*xc);

max(abs(f));
max(abs(dtf));
max(abs(dyf));

M = max(abs(dtf));
dt = 0.001;
n = 100;
E_n = (1+M)/2 * dt * M * dt * n
fprintf('With a step size of 0.001 and n = 100, the error bound is %f.\n', E_n);
fprintf('The actual error computed as %f.\n', abs(xc(end)-exact_sol(0.5)));

%Euler's method is first order because the error decreases 
%linearly with the step size. With the step size decrease by a factor of 10
%(0.01 to 0.001), the actual error decreases by the factor of 10; (0.004732
%to 0.000472)
%% Adaptive Step Size
% As mentioned in lab 2, the step size in |ode45| is adapted to a specific error 
% tolerance.
% 
% The idea of adaptive step size is to change the step size |h| to a smaller 
% number whenever the derivative of the solution changes quickly. This is done 
% by evaluating f(t,y) and checking how it changes from one iteration to the next.
%% Exercise 4
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
% 
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as in 
% exercise 1, where |h| is an initial step size. You may also want to pass the 
% function into the ODE the way |ode45| does.
% 
% Create an implementation of Euler's method by modifying your solution to exercise 
% 1. Change it to include the following:
% 
% (a) On each timestep, make two estimates of the value of the solution at the 
% end of the timestep: |Y| from one Euler step of size |h| and |Z| from two successive 
% Euler steps of size |h/2|. The difference in these two values is an estimate 
% for the error.
% 
% (b) Let |tol=1e-8| and |D=Z-Y|. If |abs(D)<tol|, declare the step to be successful 
% and set the new solution value to be |Z+D|. This value has local error |O(h^3)|. 
% If |abs(D)>=tol|, reject this step and repeat it with a new step size, from 
% (c).
% 
% (c) Update the step size as |h = 0.9*h*min(max(tol/abs(D),0.3),2)|.
% 
% Comment on what the formula for updating the step size is attempting to achieve.
% 
% 

% solution written in lab3_yookjoeu_adaptive_euler.m file
% Test function from the previous exercise to verify the adptive euler code:
f = @(t, y) 2 * t * sqrt(1 - y.^2);  % Define the function f(t, y)

t0 = 0;   % Initial time
tN = 0.5; % Final time
y0 = 0;   % Initial condition
h = 0.01; % Initial step size
[t, y] = lab3_yookjoeu_adaptive_euler(f, t0, tN, y0, h);
% Plot the solution
plot(t, y);
xlabel('t');
ylabel('y');
title('Adaptive Euler Solution');

%The formula reduces the step size (h) whenever the error (|D|) exceeds the tolerance (tol),
% ensuring a more accurate result by taking smaller steps.
%It uses a safety factor of 0.9 to avoid large jumps in step size .
%The formula also allows the step size to increase if the error is significantly smaller
% than the tolerance, speeding up computation.
%% Exercise 5
% Objective: Compare Euler to your Adaptive Euler method.
% 
% Details: Consider the IVP from exercise 3.
% 
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75| with 
% |h=0.025|.
% 
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.
% 
% (c) Plot both approximations together with the exact solution.

% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75| with 
% |h=0.025|.

xc=euler(@(t, y) 2*t*sqrt(1-y.^2), 0, 0:0.025:0.75);
xc(end)
%% 
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.

% My adavtive euler
[t_adapted, y_adapted] = lab3_yookjoeu_adaptive_euler(@(t, y) 2*t*sqrt(1-y.^2), 0, 0.75, 0, 0.025);
results_table = table(t_adapted', y_adapted', 'VariableNames', {'Time', 'Solution'});

% (c) Plot both approximations together with the exact solution.

t = 0:0.025:0.75;
figure; 
plot(t, xc , t_adapted, y_adapted);
legend('euler', 'My Adaptive', 'Location', 'Best');
title("Exercise 5: Comparing euler and adapted");
ylabel('y');
xlabel('t');
%% Exercise 6
% Objective: Problems with Numerical Methods.
% 
% Details: Consider the IVP from exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is closer 
% to the actual solution (done in 3.b)? Explain why.
% 
% (b) Plot the exact solution (from exercise 3.b), the Euler's approximation 
% (from exercise 3.a) and the adaptive Euler's approximation (from exercise 5) 
% from |t=0| to |t=1.5|.
% 
% (c) Notice how the exact solution and the approximations become very different. 
% Why is that? Write your answer as a comment.
% 
% 

% (a) From the two approximations calculated in exercise 5, which one is closer 
% to the actual solution (done in 3.b)? Explain why.

% The value from 3.b was ans = 0.2427.
% The Adaptive method is closer to the actual solution than Euler's.
% This is because Adaptive method adjusts the step size based on error tolerance: 
% it decreases the step size when the error exceeds the tolerance, 
% leading to higher accuracy.

%% 
% (b) Plot the exact solution (from exercise 3.b), Euler's approximation 
% (from exercise 3.a), and the adaptive Euler's approximation (from exercise 5) 
% from t=0 to t=1.5.

xc = euler(@(t, y) 2*t*sqrt(1-y.^2), 0, 0:0.025:1.5);
[t_adapt, y_adapt] = lab3_yookjoeu_adaptive_euler(@(t, y) 2*t*sqrt(1-y.^2), 0, 1.5, 0, 0.025);
t = 0:0.025:1.5;
y = @(t) sin(t.^2); % Exact solution
figure;
plot(t, xc, t_adapt, y_adapt, t, y(t))
title("Exercise 6: Comparing Adaptive vs Euler vs Exact");
legend("Euler", "Adaptive", "Exact", "Location", "best");
ylabel("y");
xlabel("t");

%% 
% (c) Why do the exact solution and the approximations become very different?

% Around t=1.5, the Euler, adaptive, and exact solutions deviate significantly.
% This happens because the differential equation becomes problematic for y > 1:
% the term sqrt(1 - y^2) becomes imaginary.
% Additionally, when y approaches 1, the slope f(t, y) approaches 0, 
% causing both Euler and Adaptive methods to predict y(n+1) = y(n).
% This results in a horizontal slope and the breakdown of further approximation.