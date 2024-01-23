%% Array Synthesis
%  Copyright 2017  The MathWorks, Inc.
%% Optimize pattern using only weights 
N = 8;                                            % 8 elements in a linear array
azimuth = -90:90;                                 % Define azimuth field of view
weights_d = hamming(N);                           % Create weights for desired pattern 
weights_d = weights_d/norm(weights_d);            % Normalize weights for pattern
stvmat = steervec((0:N-1)/2,azimuth);             % Generate a pattern for this example
Beam_d = weights_d'*stvmat;                       % Apply weights to desired pattern
%% Set up optimization 
objfun = @(w)norm(w'*stvmat-Beam_d);     % Define objective function used in fmincon
                                              % Goal is to minimize the norms between
                                              % the desired pattern and
                                              % resulting pattern
weights_i = ones(N,1);                        % Initial setting for array amplitudes
                                              % Serves as starting point to
                                              % optimization
weights_o = fmincon(objfun,weights_i,[],[],[],[],zeros(N,1),ones(N,1));
                                             % fmincon takes in the objfun,
                                             % the initial weights, and 
                                             % upper and lower bounds of the weights   
                                             % In this example,
                                             % 0 <= weights_o <= 1
                                             % weights_o holds the weights
                                             % which can be used to create
                                             % a beam that matches our
                                             % desired pattern
 %%                                           
figure;
plot(azimuth,mag2db([abs(Beam_d);abs(weights_o'*stvmat)]'))
legend('desired','synthesized')
ylim([-100, 50]);
%% Optimize Both Weights and Element Position
% Build Desired Pattern
N = 8;                                        % 8 elements in array
azimuth = -90:90;
weights_d = hamming(N).*steervec((0:N-1)/2,10);
stvmat = steervec((0:N-1)/2,azimuth);
Beam_d = abs(weights_d'*stvmat);
% Set up optimization
% break complex weights to real and imaginary
w_i_re = ones(N,1);                          % initial real values
w_i_im = zeros(N,1);                         % initial imaginary values 
x_i = (0:N-1)'/2;                            % initial locations - uniform
objfun = @(x)norm(abs((x(1:N)+1i*x(N+1:2*N))'*steervec(x(2*N+1:end).',azimuth))-Beam_d);
x_ini = [w_i_re;w_i_im;x_i];
[x_o, fval, ef, output] = ...
    fmincon(objfun,x_ini,[],[],[],[],[-ones(2*N,1);(0:N-1)'/2-0.2],...
     [ones(2*N,1);(0:N-1)'/2+0.2],...
    [], optimoptions(@fmincon,'Display','iter'));
%%
figure;
plot(azimuth,mag2db([Beam_d;abs((x_o(1:N)+1i*x_o(N+1:2*N))'*steervec(x_o(2*N+1:end).',azimuth))]'))
legend('desired','synthesized')
%%
scatter(x_o(2*N+1:end), zeros(1,8))
