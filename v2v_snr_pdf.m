function out = v2v_snr_pdf(veh_pos, i, j, phy, pl, geom, varargin)
% out = v2v_snr_pdf(veh_pos, i, j, phy, pl, geom, [B], [gammas], [gamma_th])

% ==== parameter.m / .mat 자동 로딩 ====
if nargin < 1 || isempty(veh_pos) || nargin < 4 || isempty(phy) || isempty(pl) || isempty(geom)
    S = load_or_run_parameter();
    if (~exist('veh_pos','var') || isempty(veh_pos)) && isfield(S,'vehicles'), veh_pos = S.vehicles; end
    if (~exist('phy','var')     || isempty(phy))     && isfield(S,'cfg') && isfield(S.cfg,'phy'),   phy  = S.cfg.phy;  end
    if (~exist('pl','var')      || isempty(pl))      && isfield(S,'cfg') && isfield(S.cfg,'pl'),    pl   = S.cfg.pl;   end
    if (~exist('geom','var')    || isempty(geom))    && isfield(S,'cfg') && isfield(S.cfg,'geom'),  geom = S.cfg.geom; end
end

% ==== 기본값/보정 ====
if ~isfield(phy,'Pt'), phy.Pt = getfield_or(phy,'Pt',20); end
if ~isfield(phy,'Gt'), phy.Gt = getfield_or(phy,'Gt',0);  end
if ~isfield(phy,'Gr'), phy.Gr = getfield_or(phy,'Gr',0);  end
if ~isfield(phy,'Pn'), phy.Pn = getfield_or(phy,'Pn',-90);end

if ~isfield(pl,'mu_PL0'),       pl.mu_PL0       = 100;  end
if ~isfield(pl,'dPL_per_block'),pl.dPL_per_block= 6;    end
if ~isfield(pl,'sigma_PL0'),    pl.sigma_PL0    = 3;    end
if ~isfield(pl,'beta_sigma'),   pl.beta_sigma   = 0.5;  end

if ~isfield(geom,'rho'),     geom.rho     = 0.05; end
if ~isfield(geom,'lv'),      geom.lv      = 4.5;  end
if ~isfield(geom,'wv'),      geom.wv      = 2.0;  end
if ~isfield(geom,'ds'),      geom.ds      = 2.0;  end
if ~isfield(geom,'lane_w'),  geom.lane_w  = 3.5;  end
if ~isfield(geom,'M'),       geom.M       = 4;    end

B = 6; gammas = linspace(-70,70,2001); gamma_th = [];
if ~isempty(varargin)
    if numel(varargin) >= 1 && ~isempty(varargin{1}), B = varargin{1}; end
    if numel(varargin) >= 2 && ~isempty(varargin{2}), gammas = varargin{2}; end
    if numel(varargin) >= 3 && ~isempty(varargin{3}), gamma_th = varargin{3}; end
end

% ==== 거리 ====
dij = veh_pos(i,1:2) - veh_pos(j,1:2);
d_tr = hypot(dij(1), dij(2));

% ==== PDF ====
[pdf_vals, p_los, p_nlos_vec] = snr_pdf(gammas, d_tr, phy, pl, geom, B);

out.gammas = gammas(:).';
out.pdf    = pdf_vals(:).';
out.d_tr   = d_tr;
out.p_los  = p_los;
out.p_nlos_vec = p_nlos_vec;

if ~isempty(gamma_th)
    mask = gammas >= gamma_th;
    out.P_success = trapz(gammas(mask), pdf_vals(mask));
end
end

% ---------- 로컬 헬퍼 ----------
function S = load_or_run_parameter()
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
        error('v2v_snr_pdf: parameter.m / parameter.mat not found.');
    end
end

function val = getfield_or(S, name, defaultVal)
    if isfield(S,name), val = S.(name); else, val = defaultVal; end
end

function [mu_gamma_k, sigma_gamma_k] = snr_gauss_params(k, phy, pl)
    mu_PL_k = pl.mu_PL0 + pl.dPL_per_block * k;
    sigma_PL_k = pl.sigma_PL0 * sqrt(1.0 + pl.beta_sigma * k);
    mu_gamma_k = (phy.Pt + phy.Gt + phy.Gr) - mu_PL_k - phy.Pn;
    sigma_gamma_k = sigma_PL_k;
end

function P = P_occ_slot(length_m, rho)
    P = 1.0 - exp(-rho * length_m);
end

function Ns = Ns_same_lane(d_tr, geom)
    d_a = geom.lv + geom.ds;
    d_eff = max(d_tr - geom.lv, 0.0);
    Ns = max(floor(d_eff / d_a), 0);
end

function p = same_lane_p_k(d_tr, geom, k)
    Ns = Ns_same_lane(d_tr, geom);
    if k > Ns, p = 0.0; return; end
    d_a = geom.lv + geom.ds;
    P_a = P_occ_slot(d_a, geom.rho);
    p = nchoosek(Ns, k) * (P_a^k) * ((1.0 - P_a)^(Ns - k));
end

function [d_b, d_c] = diff_lane_slot_lengths(d_tr, geom)
    Dy = geom.lane_w;
    root = max(d_tr^2 - Dy^2, 0.0);
    if Dy > 0
        d_b = geom.wv * sqrt(root) / (2.0 * Dy);
        d_c = geom.wv * sqrt(root) / Dy + geom.lv;
    else
        d_b = 0.0;
        d_c = geom.lv;
    end
end

function p = prob_K_in_nplus1_slots(k, Ps)
    n1 = numel(Ps);
    if k > n1, p = 0.0; return; end
    idx = 1:n1;
    C = nchoosek(idx, k);
    total = 0.0;
    for r = 1:size(C,1)
        A = C(r,:);
        pA = 1.0;
        for t = 1:n1
            if any(A==t), pA = pA * Ps(t); else, pA = pA * (1.0 - Ps(t)); end
        end
        total = total + pA;
    end
    p = total;
end

function p = different_lane_p_k(d_tr, geom, k)
    M = geom.M;
    if M < 2, p = 0.0; return; end
    [d_b, d_c] = diff_lane_slot_lengths(d_tr, geom);
    P_b = P_occ_slot(d_b, geom.rho);
    P_c = P_occ_slot(d_c, geom.rho);
    term_neigh = 0.0;
    if k <= 2
        term_neigh = (2*(M-1)/M^2) * (nchoosek(2, k) * (P_b^k) * ((1-P_b)^(2-k)));
    end
    term_far = 0.0;
    for n = 2:M-1
        weight = 2*(M - n)/M^2;
        Ps = [P_b, repmat(P_c, 1, n-1), P_b];
        term_far = term_far + weight * prob_K_in_nplus1_slots(k, Ps);
    end
    p = term_neigh + term_far;
end

function p = P_NLoS_k(d_tr, geom, k)
    p = different_lane_p_k(d_tr, geom, k) + (1.0/geom.M) * same_lane_p_k(d_tr, geom, k);
end

function p = P_LoS(d_tr, geom, B)
    acc = 0.0;
    for k = 1:B, acc = acc + P_NLoS_k(d_tr, geom, k); end
    p = max(1.0 - acc, 0.0);
end

function [pdf_vals, p_los, p_nlos_vec] = snr_pdf(gamma_vals, d_tr, phy, pl, geom, B)
    gamma_vals = gamma_vals(:).';
    pdf_vals = zeros(size(gamma_vals));
    p_nlos_vec = zeros(1,B);
    [mu0, sig0] = snr_gauss_params(0, phy, pl);
    p_los = P_LoS(d_tr, geom, B);
    pdf_vals = pdf_vals + p_los * normal_pdf(gamma_vals, mu0, sig0);
    for k = 1:B
        pk = P_NLoS_k(d_tr, geom, k);
        p_nlos_vec(k) = pk;
        if pk <= 0, continue; end
        [muk, sigk] = snr_gauss_params(k, phy, pl);
        pdf_vals = pdf_vals + pk * normal_pdf(gamma_vals, muk, sigk);
    end
    area = trapz(gamma_vals, pdf_vals);
    if area > 0, pdf_vals = pdf_vals / area; end
end

function y = normal_pdf(x, mu, sigma)
    y = (1./(sigma*sqrt(2*pi))) .* exp(-0.5*((x - mu)./sigma).^2);
end
