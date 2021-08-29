clear
clc

% This is a MATLAB file explaining:
    % 1. How to use function handles
    % 2. How to use MATLAB's optimization function fminsearch
    
% The mathematical function used in this example is:
    % y=x^2-4x
    
% Our goal is to find the minimum of this function using fminsearch


%% Call the function: y=x^2-4x for a given series of x-values and plot it

% create the series of x values
x_val=(-5:0.5:5)';

% call the function computeSquare.m
y_val=mini_ex(x_val);

% Plot the function to eye-ball the minimum
figure;
plot(x_val,y_val);
xlabel('x values');
ylabel('y values');
title('Evaluate the minimum of $y=x^2-4x$');
close all


% Note: Eyeballing the graph (and also with very simple math), the minimum
% of the graph is at x=2
% Now, let's see whether we get the same result using fminsearch.

%% Creating function handles

% We could also call the function mini_ex.m in a different way using
% function handles (for more details on function handles, take a look at:
% https://de.mathworks.com/help/matlab/matlab_prog/creating-a-function-handle.html) 

func=@mini_ex;

% Take a look at the workspace, now we have a new variable called "func", which is our function handle for mini_ex.m
% If we give func any value as an input, it will return the x^2-4x of that value as follows:

a=3;
b=func(a);

% This is very useful, as we can now hand in our function handle, func, to MATLAB's optimization function fminsearch.

%% Using fminsearch

% fminsearch has the following syntax:
%   [x,fval, exitflag, output]=fminsearch(fun,x0)
%   Input arguments:    fun      --> This is the function for which you
%                                    want to find the local minimum. Here we can put our function handle
%                                    func. Alternatively, we could have put @mini_ex directly.
%                       x0       --> This is the starting value at which fminsearch will start searching for the local minimum of fun.
%   Output arguments:   x        --> This is the local minimum of the function fun. 
%                       fval     --> This is the value of the function fun at the local minimum x.
%                       exitflag --> returns a value that describes the exit condition. 1: means the function did converge
%                                    and find a solution x. 0: means the maximum number of iterations was reached before a solution could
%                                    be found. -1: means the algorithm was terminated by the output function.
%                       output   --> returns a structure with information about the optimization process  

% This information and many other details can be found in the excellent
% documentation on fminsearch in MATLAB

% suggest a starting value for fminsearch
start_val=5;

% use fminsearch to find the minimum of y=x^2
[min,fval,exitflag,output]=fminsearch(func,start_val);

% Note that we would get the same result, if we change the syntax to:
[min2,fval2,exitflag2,output2]=fminsearch(@(start_val)mini_ex(start_val),start_val);

% The syntax in line 76 is especially useful in our assignment. The syntax
% for the function handle generally is then: @(variableInput)funName(variableInput,fixedInput_1,fixedInput2,...)
% This means that the optimizer, fminsearch, will try to find the minimum
% of the function funName with respect to the variableInput, while keeping
% all other inputs to the function fixed, i.e. fixedInput_1,
% fixedInput2,...,etc.

% In our mini_ex.m function, there is only 1 input which is, of course then, the variable input. 
% The optimizer fminsearch, then tries to minimize the function mini_ex
% with respect to that variable input.









