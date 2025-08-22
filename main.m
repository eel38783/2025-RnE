clear;close all;clc;

Nx = 4;
Ny = 4;

K = [3,6,6,3]; %vehicular set [K1, K2, K3, K4]
rho = randi([2,8], 15, 1); %allocated resource
Ds = randi([2,10], 15, 1); %amount of data
C = randi([2,10], 15, 1); %required computation
t_tole = 10 * ones(1,15);
t_hold = 7 * ones(15, 15, 2);

action_space = linspace(-5, 5, 16) * (pi / 10);
theta = action_space(randi(numel(action_space), 1, 16)); %steering vector
veh_pos = [rand(20, 2) * 200, zeros(20,1)]; % vechile location
UAV = [10,10,20]; % fix UAV Height as 20m <-- use optimized coordinate of UAV by DDQN


psi2 = V_RIS(Nx, Ny, theta, UAV, veh_pos); %V2UAV Mode SNR
psi1 = V2V_SNR(veh_pos,K);

psi = Psi(veh_pos, psi1, psi2); % transmission rate

T = Opt_func(K, rho, Ds, C, psi, t_tole, t_hold); %Optimized Delay Time

disp(T);