function pdf_vals = snr_pdf(d_tr)
% SNR_PDF  V2V 링크의 SNR(dB) 확률밀도함수 계산
% 사용법:
%   pdf = snr_pdf(gammas, d_tr, phy, pl, geom)         % B 기본값=6
%   [pdf, p_los, p_nlos_vec] = snr_pdf(gammas, d_tr, phy, pl, geom, B)
%
% 입력:
%   gamma_vals : 1xN 또는 Nx1, SNR(dB) 축
%   d_tr       : Tx-Rx 거리 [m]
%   phy        : struct with fields {Pt, Gt, Gr, Pn}  (dB 단위)
%   pl         : struct with fields {mu_PL0, dPL_per_block, sigma_PL0, beta_sigma}
%   geom       : struct with fields {rho, lv, wv, ds, lane_w, M}
%   B          : (옵션) 최대 차폐 블록 수, 기본 6
%
% 출력:
%   pdf_vals    : gamma_vals와 같은 길이의 PDF 값 (정규화되어 적분=1)
%   p_los       : LoS 확률
%   p_nlos_vec  : k=1..B 각 NLoS 확률 벡터
    gamma_vals = linspace(-70, 70, 2001);
    phy = ;
    pl = ;
    geom = ;
    B = 6;

    if nargin < 6 || isempty(B), B = 6; end

    % 행 벡터로 정리
    gamma_vals = gamma_vals(:).';
    pdf_vals   = zeros(size(gamma_vals));
    p_nlos_vec = zeros(1, B);

    % k=0 (LoS) 가우시안 파라미터
    [mu0, sig0] = snr_gauss_params(0, phy, pl);

    % LoS probability
    p_los = P_LoS(d_tr, geom, B);
    pdf_vals = pdf_vals + p_los * normal_pdf(gamma_vals, mu0, sig0);

    % k=1..B 혼합
    for k = 1:B
        pk = P_NLoS_k(d_tr, geom, k);
        p_nlos_vec(k) = pk;
        if pk <= 0, continue; end
        [muk, sigk] = snr_gauss_params(k, phy, pl);
        pdf_vals = pdf_vals + pk * normal_pdf(gamma_vals, muk, sigk);
    end

    % 정규화 (수치 오차 보정)
    area = trapz(gamma_vals, pdf_vals);
    if area > 0
        pdf_vals = pdf_vals / area;
    end
end

% sub functions 

function [mu_gamma_k, sigma_gamma_k] = snr_gauss_params(k, phy, pl)
% k차 차폐 시의 SNR(dB) 가우시안 근사 파라미터
    mu_PL_k      = pl.mu_PL0 + pl.dPL_per_block * k;
    sigma_PL_k   = pl.sigma_PL0 * sqrt(1.0 + pl.beta_sigma * k);
    mu_gamma_k   = (phy.Pt + phy.i + phy.Gr) - mu_PL_k - phy.Pn;
    sigma_gamma_k= sigma_PL_k;
end

function P = P_occ_slot(length_m, rho)
% 길이 length_m 슬롯 점유 확률 (포아송 가정)
    P = 1.0 - exp(-rho * length_m);
end

function Ns = Ns_same_lane(d_tr, geom)
% 동일 차선에서 가능한 차폐 차량 슬롯 수
    d_a = geom.lv + geom.ds;
    d_eff = max(d_tr - geom.lv, 0.0);
    Ns = max(floor(d_eff / d_a), 0);
end

function p = same_lane_p_k(d_tr, geom, k)
% 동일 차선에서 k대가 차폐할 확률
    Ns = Ns_same_lane(d_tr, geom);
    if k > Ns, p = 0.0; return; end
    d_a = geom.lv + geom.ds;
    P_a = P_occ_slot(d_a, geom.rho);
    p = nchoosek(Ns, k) * (P_a^k) * ((1.0 - P_a)^(Ns - k));
end

function [d_b, d_c] = diff_lane_slot_lengths(d_tr, geom)
% 이웃/원거리 차선에서의 유효 슬롯 길이
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
% n+1개의 독립 슬롯 확률 Ps에서 정확히 k개 점유될 확률
    n1 = numel(Ps);
    if k > n1, p = 0.0; return; end
    if k == 0
        p = prod(1 - Ps);
        return;
    end
    idx = 1:n1;
    C = nchoosek(idx, k);
    total = 0.0;
    for r = 1:size(C,1)
        A = C(r,:);
        sel = false(1,n1); sel(A) = true;
        pA = prod(Ps(sel)) * prod(1 - Ps(~sel));
        total = total + pA;
    end
    p = total;
end

function p = different_lane_p_k(d_tr, geom, k)
% 다른 차선들로 인한 k대 차폐 확률
    M = geom.M;
    if M < 2, p = 0.0; return; end

    [d_b, d_c] = diff_lane_slot_lengths(d_tr, geom);
    P_b = P_occ_slot(d_b, geom.rho);
    P_c = P_occ_slot(d_c, geom.rho);

    % 바로 이웃한 양쪽 차선 (총 2개 슬롯)
    term_neigh = 0.0;
    if k <= 2
        term_neigh = (2*(M-1)/M^2) * (nchoosek(2, k) * (P_b^k) * ((1-P_b)^(2-k)));
    end

    % 더 먼 차선들
    term_far = 0.0;
    for n = 2:(M-1)   % Python range(2,M)에 대응
        weight = 2*(M - n)/M^2;
        Ps = [P_b, repmat(P_c, 1, n-1), P_b];  % 길이 n+1
        term_far = term_far + weight * prob_K_in_nplus1_slots(k, Ps);
    end

    p = term_neigh + term_far;
end

function p = P_NLoS_k(d_tr, geom, k)
% k대 차폐에 의한 NLoS 확률(동일 차선 + 다른 차선)
    p = different_lane_p_k(d_tr, geom, k) + (1.0/geom.M) * same_lane_p_k(d_tr, geom, k);
end

function p = P_LoS(d_tr, geom, B)
% LoS 확률 = 1 - sum_{k=1..B} P_NLoS_k
    acc = 0.0;
    for k = 1:B
        acc = acc + P_NLoS_k(d_tr, geom, k);
    end
    p = max(1.0 - acc, 0.0);
end

function y = normal_pdf(x, mu, sigma)
% 정규분포 PDF (도구상자 없이)
    y = (1./(sigma*sqrt(2*pi))) .* exp(-0.5*((x - mu)./sigma).^2);
end

function psi1 = V2V_SNR(veh_pos, i, j)
    pdf_vals = snr_pdf(distance(veh_pos(i),veh_pos(j)));
    cdf_vals = cumtrapz(gamma_vals, pdf_vals);
    cdf_vals = cdf_vals / cdf_vals(end);

    N = 1000;
    u = rand(1,N);
    psi1 = interp1(cdf_vals, gamma_vals, u, 'linear', 'extrap');
end

