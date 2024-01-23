%% Array Synthesis
%  Copyright 2017  The MathWorks, Inc.
%% Optimize pattern using only weights 
format long g
N_x = 64;                                            % 8 elements in a linear array
N_y = 64;
azimuth = -90:1:90;                                 % Define azimuth field of view
elevation = -90:1:90;
%elevation = 0;
azrep = repmat(azimuth, 1, size(elevation,2));
elrep = repelem(elevation, 1, size(azimuth,2));
angs = [azrep; elrep];
%disp(angs);
azrand = 180 * rand(1,1000) - 90;
elrand = 180 * rand(1,1000) - 90;
angsrand = [azrand;elrand];
angstrain = angsrand;
angstest = angs;
%weights_d = hamming(N);                           % Create weights for desired pattern 
%weights_d = ones(N,1);3j6ilBYnhcjl
%weights_d = weights_d/norm(weights_d);            % Normalize weights for pattern
locs_x = (0:N_x-1)/4 - ((N_x-1)/8);
locs_y = (N_y-1:-1:0)/4 - ((N_y-1)/8);
locs_x = imresize(locs_x, [1 N_x*N_y], 'nearest');
locs_y = repmat(locs_y,1, N_x);
locs = [locs_x; locs_y];
stvmat = steervec(locs,angstrain);             % Generate a pattern for this example
%Beam_d = weights_d'*stvmat;                       % Apply weights to desired pattern
%disp(size(Beam_d));
Beam_d = ones(1,size(angstrain,2));
%% Set up optimization
objfun = @(ph)max(abs(exp(1i*ph)'*stvmat))- min(abs(exp(1i*ph)'*stvmat)); 
%objfun = @(ph)(1.0/max(abs(exp(1i*ph)'*stvmat))); 

weights_i = 1*ones(N_x*N_y,1);                        % Initial setting for array amplitudes
                                              % Serves as starting point to
                                              % optimization
options = optimoptions(@fmincon,'Algorithm','active-set','Display','off');

tic;

weights_o = fmincon(objfun,weights_i,[],[],[],[],-1*pi*ones(N_x*N_y,1),1*pi*ones(N_x*N_y,1), [], options);
                                             % fmincon takes in the objfun,
                                             % the initial weights, and 
                                             % upper and lower bounds of the weights   
                                             % In this example,
                                             % 0 <= weights_o <= 1
                                             % weights_o holds the weights
                                             % which can be used to create
                                             % a beam that matches our
                                             % desired pattern
                                             
after = toc;
%%                 
stvmattest = steervec(locs,angstest);
disp(weights_o);
disp(after);
disp(mean(abs(exp(1i*weights_o)'*stvmattest)))
return;

figure;
t = tiledlayout(2,1);
nexttile;
plot(1:size(angstrain,2),mag2db([abs(Beam_d);abs(exp(1i*weights_o)'*stvmat)]'))
%ylim([-50, 30]);
legend('desired','synthesized')
nexttile;
plot(1:size(angstest,2),mag2db([abs(exp(1i*weights_o)'*stvmattest)]'))
%ylim([-50, 30]);
return;
