clear; clc; rng(0);

K = 20;
vehicles = [rand(K,2)*200, 1.5*ones(K,1)];  % Kx3 [x,y,H]

% cfg (psi 계산용)
cfg = struct();
cfg.N = 32; cfg.z_uav = 60; cfg.lambda_c = 0.01; cfg.d_elem = 0.5*cfg.lambda_c;
cfg.rho0_db = -20; cfg.uav_xy = [100,100];
cfg.PT_dbm = 0; cfg.GT_db = 10; cfg.GR_db = 10;
cfg.noise_dbm = -80; cfg.BW_hz = 100e6; cfg.N0_WHz = 4e-21;

cfg.phy  = struct('Pt',20,'Gt',0,'Gr',0,'Pn',-90);
cfg.pl   = struct('mu_PL0',100,'dPL_per_block',6,'sigma_PL0',3,'beta_sigma',0.5);
cfg.geom = struct('rho',0.05,'lv',4.5,'wv',2.0,'ds',2.0,'lane_w',3.5,'M',4);
cfg.B_v2v = 6; cfg.gammas = linspace(-70,70,2001);

j0 = 0;                          % 수식의 j=0
Ppsi = psi(vehicles, cfg, j0);   % KxKx2 (psi(:,:,1)=†, psi(:,:,2)=0)

% params
params.K_TV = 1:15;
params.K_SV = 0:14;              % 0부터 시작 → solve 함수가 내부에서 +1 매핑
params.K1   = 1:5;
params.j0   = j0;

params.c = randi([2,10], K, 1);
params.d = randi([2,10], K, 1);
params.rho_max = randi([2,8], K, 1);

params.psi = Ppsi;               % 또는 psid/psi0로 별도 전달 가능

% 시간 제약(예시 값)
params.t_tol    = 0.50 + 0.10*rand(K,1);
params.t_hold0  = 0.40 + 0.10*rand(K,1);
params.t_holdd  = 0.45 + 0.10*rand(K,1);
params.t_holdij = 0.40 + 0.20*rand(K,K);
params.t_holdij(1:K+1:end) = 1e6;

% SINR 행렬(선형). 예시로 psi의 큰 값 사용
Glin = max(Ppsi(:,:,1), Ppsi(:,:,2));
Glin(1:K+1:end) = 1e6;
params.gamma = Glin;

params.gamma_VUE = 10^(5/10);    % 5 dB

% CVX 경로 세팅 (처음 1회)
% addpath('C:\cvx'); cvx_setup;

out = solve_minT_eps_cvx(params);
disp(out.status);
disp(out.T_value);
