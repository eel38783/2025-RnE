function P = psi(vehicles, cfg, j0, v, quantBits)
% psi 함수 — vehicles, cfg, j0 등이 주어지지 않으면 parameter.m에서 자동 로딩

% ================== parameter.m 로딩 ==================
if nargin < 1 || isempty(vehicles) || nargin < 2 || isempty(cfg)
    % parameter.m 또는 parameter.mat 로드
    if exist('parameter.mat','file')
        S = load('parameter.mat');
        if isfield(S,'params') && isstruct(S.params), S = S.params; end
    elseif exist('parameter.m','file')
        run('parameter.m');
        if exist('params','var') && isstruct(params)
            S = params;
        else
            vv = whos;
            S  = struct();
            skip = {'S','vv'};
            for k = 1:numel(vv)
                n = vv(k).name;
                if any(strcmp(n, skip)), continue; end
                S.(n) = eval(n);
            end
        end
    else
        error('psi: parameter.m 또는 parameter.mat을 찾을 수 없습니다.');
    end

    % 필요한 필드 꺼내기
    if (~exist('vehicles','var') || isempty(vehicles))
        if isfield(S,'vehicles'), vehicles = S.vehicles;
        else, error('psi: vehicles 데이터가 없습니다.'); end
    end
    if (~exist('cfg','var') || isempty(cfg))
        if isfield(S,'cfg'), cfg = S.cfg;
        else, error('psi: cfg 데이터가 없습니다.'); end
    end
    if nargin < 3 || isempty(j0)
        if isfield(S,'j0'), j0 = S.j0; else, j0 = 0; end
    end
end

% -------------------- defaults --------------------
if nargin < 3 || isempty(j0),      j0 = 0;                 end
if nargin < 4,                     v = [];                 end
if nargin < 5,                     quantBits = [];         end

if ~isfield(cfg, 'N'),          cfg.N = 32;                     end
if ~isfield(cfg, 'z_uav'),      cfg.z_uav = 60;                 end
if ~isfield(cfg, 'lambda_c'),   cfg.lambda_c = 0.01;            end
if ~isfield(cfg, 'd_elem'),     cfg.d_elem = 0.5*cfg.lambda_c;  end
if ~isfield(cfg, 'rho0_db'),    cfg.rho0_db = -20;              end
if ~isfield(cfg, 'uav_xy'),     cfg.uav_xy = [100,100];         end
if ~isfield(cfg, 'PT_dbm'),     cfg.PT_dbm = 0;                 end
% 호환 필드명 처리
if ~isfield(cfg, 'GT_db'),      cfg.GT_db = getfield_or(cfg,'GT_dB',10); end
if ~isfield(cfg, 'GR_db'),      cfg.GR_db = getfield_or(cfg,'GR_dB',10); end
if ~isfield(cfg, 'noise_dbm'),  cfg.noise_dbm = -80;            end
if ~isfield(cfg, 'BW_hz'),      cfg.BW_hz = 100e6;              end
if ~isfield(cfg, 'N0_WHz'),     cfg.N0_WHz = 4e-21;             end

% ---- V2V PDF 기본 파라미터 ----
if ~isfield(cfg,'phy'),   cfg.phy  = struct('Pt',20,'Gt',0,'Gr',0,'Pn',-90); end
if ~isfield(cfg,'pl'),    cfg.pl   = struct('mu_PL0',100,'dPL_per_block',6,'sigma_PL0',3,'beta_sigma',0.5); end
if ~isfield(cfg,'geom'),  cfg.geom = struct('rho',0.05,'lv',4.5,'wv',2.0,'ds',2.0,'lane_w',3.5,'M',4); end
if ~isfield(cfg,'B_v2v'), cfg.B_v2v = 6; end
if ~isfield(cfg,'gammas'), cfg.gammas = linspace(-70,70,2001); end

% -------------------- 준비 --------------------
K = size(vehicles,1);
P = nan(K,K,2);
PT_W = 10.^((cfg.PT_dbm - 30)/10);
j0_col = j0 + 1;

% 상수 (후처리용)
eps_small = 1e-9;  % 너무 작은 값 하한
big_val  = 1e9;    % 대각/불사용 자리용 큰 값

% -------------------- psi(:,:,1) = psi_{i,j,†} --------------------
for i = 1:K
    for j_col = 1:K
        if i == j_col
            % dagger 대각은 사용 안 함: 0 또는 big_val 아무거나 가능
            P(i,j_col,1) = 0;
            continue;
        end

        [h_eff, gain_lin, snr_db] = effective_channel_via_ris( ...
            vehicles, i, j_col, cfg, v, quantBits);

        if j_col == j0_col
            % RIS 링크 SNR (선형)
            P(i,j_col,1) = PT_W * gain_lin / (cfg.BW_hz * cfg.N0_WHz);
        else
            P(i,j_col,1) = 10.^(snr_db/10);
        end
    end
end

% -------------------- psi(:,:,2) = psi_{i,j,0} (V2V & j0 처리) --------------------
for i = 1:K
    for j_col = 1:K
        if i == j_col
            % 대각은 분모로 쓰지 않게 매우 큰 값으로
            P(i,j_col,2) = big_val;
            continue;
        end

        if j_col == j0_col
            % j0 컬럼도 반드시 유한한 양수로 채워야 함 (RIS 경로 기반)
            [h_eff, gain_lin, ~] = effective_channel_via_ris( ...
                vehicles, i, j_col, cfg, v, quantBits);
            P(i,j_col,2) = max(PT_W * gain_lin / (cfg.BW_hz * cfg.N0_WHz), eps_small);
        else
            out = v2v_snr_pdf(vehicles, i, j_col, cfg.phy, cfg.pl, cfg.geom, cfg.B_v2v, cfg.gammas);
            % E[10^(γ/10)] = 선형 SNR 기대값
            P(i,j_col,2) = trapz(out.gammas, (10.^(out.gammas/10)) .* out.pdf);
        end
    end
end

% --------- sanitize P(:,:,1) / P(:,:,2) ---------
P(:,:,1) = real(P(:,:,1));
P(:,:,2) = real(P(:,:,2));

% NaN/Inf 치환 (임시 변수 사용)
tmp1 = P(:,:,1);
tmp2 = P(:,:,2);

mask1 = ~isfinite(tmp1);
mask2 = ~isfinite(tmp2);

tmp1(mask1) = eps_small;
tmp2(mask2) = eps_small;

% 하한 적용 (0 또는 음수 방지)
tmp1(tmp1 <= 0) = eps_small;
tmp2(tmp2 <= 0) = eps_small;

% 대각은 큰 값으로 잠그기 (psi(:,:,2)만)
for ii = 1:K
    tmp2(ii,ii) = big_val;
end

% 반영
P(:,:,1) = tmp1;
P(:,:,2) = tmp2;

end % <-- psi 함수 'end' 유지

% ================== 내부 헬퍼들 ==================
function val = getfield_or(S, name, defaultVal)
if isfield(S, name), val = S.(name); else, val = defaultVal; end
end

function [h_eff, gain_lin, snr_db] = effective_channel_via_ris( ...
    vehicles, tx_idx, rx_idx, cfg, v, quantBits)

h_VR = channel_vehicle_to_ris(vehicles, tx_idx, cfg);
h_RV = channel_ris_to_vehicle(vehicles, rx_idx, cfg);

if isempty(v)
    theta = -angle(h_VR) - angle(h_RV);
    v = exp(1j*theta);
else
    v = v(:);
    if numel(v) ~= cfg.N
        error('RIS phase length mismatch (N=%d).', cfg.N);
    end
end

if ~isempty(quantBits)
    L = 2^quantBits;
    step = 2*pi/L;
    ang  = angle(v);
    angq = round(ang/step)*step;
    v    = exp(1j*angq);
end

h_eff    = sum(h_RV .* v .* h_VR);
gain_lin = abs(h_eff)^2;

Pr_dbm = cfg.PT_dbm + cfg.GT_db + cfg.GR_db + 20*log10(abs(h_eff) + 1e-30);
snr_db  = Pr_dbm - cfg.noise_dbm;
end

function h = channel_vehicle_to_ris(vehicles, idx, cfg)
x = vehicles(idx,1);  y = vehicles(idx,2);  H = vehicles(idx,3);
w_xy = [x, y];  q = cfg.uav_xy(:).';  Z = cfg.z_uav;

vec_xy = q - w_xy;
d = sqrt(sum(vec_xy.^2) + (Z - H)^2) + 1e-12;
path   = path_gain_los(d, cfg.rho0_db);
cosphi = cos_angle_from_xy(vec_xy);
a = steering_vec(cosphi, cfg);
h = path .* a(:);
end

function h = channel_ris_to_vehicle(vehicles, idx, cfg)
x = vehicles(idx,1);  y = vehicles(idx,2);  H = vehicles(idx,3);
w_xy = [x, y];  q = cfg.uav_xy(:).';  Z = cfg.z_uav;

d = sqrt(sum((q - w_xy).^2) + (Z - H)^2) + 1e-12;
path   = path_gain_los(d, cfg.rho0_db);
cosphi = cos_angle_from_xy(w_xy - q);
a = steering_vec(cosphi, cfg);
h = path .* a(:);
end

function path = path_gain_los(d, rho0_db)
rho0_lin = 10.^(rho0_db/10);
path = sqrt(rho0_lin) / d;
end

function a = steering_vec(cosphi, cfg)
n = (0:cfg.N-1).';
a = exp(-1j * 2*pi * cfg.d_elem / cfg.lambda_c .* n * cosphi);
end

function c = cos_angle_from_xy(vec2)
nrm = hypot(vec2(1), vec2(2));
if nrm < 1e-12
    c = 1.0;
else
    c = vec2(1) / nrm;
end
end
