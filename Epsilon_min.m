function out = Epsilon_min(params)
% (상단 주석 생략) — parameter.m 자동 로딩 지원

% ---------- 입력 로딩/정규화 ----------
if nargin==0 || isempty(params)
    params = 'parameter';
end

if ischar(params) || (isstring(params) && isscalar(params))
    name = char(params);
    [pth,base,ext] = fileparts(name);
    if isempty(ext)
        cand_mat = fullfile(pth,[base '.mat']);
        cand_m   = fullfile(pth,[base '.m']);
        if exist(cand_mat,'file')
            S = load(cand_mat);
        elseif exist(cand_m,'file')
            run(cand_m);
            if exist('params','var') && isstruct(params)
                S = params;
            else
                vv = whos;
                S  = struct();
                skip = {'name','pth','base','ext','cand_mat','cand_m','S','vv','skip'};
                for k = 1:numel(vv)
                    n = vv(k).name;
                    if any(strcmp(n, skip)), continue; end
                    S.(n) = eval(n);
                end
            end
        else
            error('Epsilon_min: ''%s(.m/.mat)''을 찾을 수 없습니다.', name);
        end
    elseif strcmpi(ext,'.mat')
        S = load(name);
    elseif strcmpi(ext,'.m')
        run(name);
        if exist('params','var') && isstruct(params)
            S = params;
        else
            vv = whos;
            S  = struct();
            skip = {'name','pth','base','ext','S','vv','skip'};
            for k = 1:numel(vv)
                n = vv(k).name;
                if any(strcmp(n, skip)), continue; end
                S.(n) = eval(n);
            end
        end
    else
        error('Epsilon_min: 지원하지 않는 파일 확장자입니다: %s', ext);
    end

    if isstruct(S) && ~isfield(S,'params')
        params = S;
    elseif isfield(S,'params') && isstruct(S.params)
        params = S.params;
    else
        params = S;
    end
elseif ~isstruct(params)
    error('Epsilon_min: params는 struct 또는 .m/.mat 경로여야 합니다.');
end

% ---------- psi 마련(없으면 생성/결합) ----------
if ~isfield(params,'psi') || isempty(params.psi)
    if isfield(params,'psid') && isfield(params,'psi0') && ~isempty(params.psid) && ~isempty(params.psi0)
        K = size(params.psid,1);
        params.psi = nan(K,K,2);
        params.psi(:,:,1) = params.psid;
        params.psi(:,:,2) = params.psi0;
    elseif isfield(params,'vehicles') && isfield(params,'cfg') && ~isempty(params.vehicles)
        if ~isfield(params,'j0') || isempty(params.j0), params.j0 = 0; end
        params.psi = psi(params.vehicles, params.cfg, params.j0);
    else
        error('psi가 없고 (psid,psi0) 또는 (vehicles,cfg)도 없어 psi를 만들 수 없습니다.');
    end
end

% ---------- 기본값 채우기 ----------
K    = size(params.psi,1);
BigM = 1e9;

if ~isfield(params,'K_TV') || isempty(params.K_TV), params.K_TV = 1:K; end
if ~isfield(params,'K_SV') || isempty(params.K_SV), params.K_SV = 0:K-1; end
if ~isfield(params,'K1')   || isempty(params.K1),   params.K1   = [];  end
if ~isfield(params,'j0')   || isempty(params.j0),   params.j0   = 0;   end

if ~isfield(params,'c')       || isempty(params.c),       params.c       = ones(K,1); end
if ~isfield(params,'d')       || isempty(params.d),       params.d       = ones(K,1); end
if ~isfield(params,'rho_max') || isempty(params.rho_max), params.rho_max = ones(K,1); end

if ~isfield(params,'t_tol')   || isempty(params.t_tol),   params.t_tol   = BigM*ones(K,1); end
if ~isfield(params,'t_hold0') || isempty(params.t_hold0), params.t_hold0 = BigM*ones(K,1); end
if ~isfield(params,'t_holdd') || isempty(params.t_holdd), params.t_holdd = BigM*ones(K,1); end
if ~isfield(params,'t_holdij')|| isempty(params.t_holdij),params.t_holdij= BigM*ones(K,K); end
params.t_tol   = params.t_tol(:);   params.t_tol   = params.t_tol(1:K);
params.t_hold0 = params.t_hold0(:); params.t_hold0 = params.t_hold0(1:K);
params.t_holdd = params.t_holdd(:); params.t_holdd = params.t_holdd(1:K);
if ~isequal(size(params.t_holdij),[K K])
    Th = BigM*ones(K,K);  sz = size(params.t_holdij);
    Th(1:min(K,sz(1)),1:min(K,sz(2))) = params.t_holdij(1:min(K,sz(1)),1:min(K,sz(2)));
    params.t_holdij = Th;
end

if ~isfield(params,'gamma')     || isempty(params.gamma),     params.gamma     = zeros(K,K); end
if ~isfield(params,'gamma_VUE') || isempty(params.gamma_VUE), params.gamma_VUE = 0;         end

% ---------- psid/psi0 언팩 ----------
P = params.psi;
psid = P(:,:,1);
if size(P,3) >= 2 && any(~isnan(P(:,:,2)),'all')
    psi0 = P(:,:,2);
elseif isfield(params,'psi0') && ~isempty(params.psi0)
    psi0 = params.psi0;
else
    error('psi(:,:,2) 또는 params.psi0가 필요합니다 (psi_{i,j,0}).');
end

[nK1, nK2] = size(psid);
assert(nK1==nK2, 'psid must be square KxK.');
nK = nK1;

% ---------- 집합/라벨 ----------
K_TV = params.K_TV(:)';     
K1   = params.K1(:)';

K_SV_lbl = params.K_SV(:)';  
if any(K_SV_lbl==0)
    K_SV = K_SV_lbl + 1;     
else
    K_SV = K_SV_lbl;
end

has_j0 = isfield(params,'j0') && ~isempty(params.j0);
if has_j0
    j0_lbl = params.j0;
    if j0_lbl==0, j0 = 1; else, j0 = j0_lbl; end
else
    j0 = [];
end

assert(all(K_TV>=1 & K_TV<=nK), 'K_TV indices out of range.');
assert(isempty(K1) || all(K1>=1 & K1<=nK), 'K1 indices out of range.');
assert(all(K_SV>=1 & K_SV<=nK), 'K_SV indices out of range.');

% ---------- 파라미터 벡터 ----------
c        = params.c(:);       c = c(1:nK);
d        = params.d(:);       d = d(1:nK);
rho_max  = params.rho_max(:); rho_max = rho_max(1:nK);

t_tol    = params.t_tol(:);    t_tol    = t_tol(1:nK);
t_hold0  = params.t_hold0(:);  t_hold0  = t_hold0(1:nK);
t_holdd  = params.t_holdd(:);  t_holdd  = t_holdd(1:nK);
t_holdij = params.t_holdij;    
gamma     = params.gamma;      
gamma_VUE = params.gamma_VUE;

% ---------- 유효성 검사 / 자동 보정 옵션 ----------
allow_auto_clip = isfield(params,'allow_zero_den') && params.allow_zero_den;
eps_den = 1e-12;

bad_rho = ~isfinite(rho_max) | rho_max<=0;
bad_psi0= ~isfinite(psi0)    | psi0==0;
bad_psid= ~isfinite(psid)    | psid==0;

if any(bad_rho)
    if allow_auto_clip
        rho_max(bad_rho) = eps_den;
    else
        error('rho_max에 0/NaN/Inf 또는 비양수 값이 있습니다. (분모가 됨)');
    end
end
if any(bad_psi0,'all')
    if allow_auto_clip
        z = abs(psi0)<eps_den | ~isfinite(psi0);
        psi0(z) = eps_den;
    else
        error('psi0에 0/NaN/Inf가 있습니다. (분모가 됨)');
    end
end
if any(bad_psid,'all')
    if allow_auto_clip
        z = abs(psid)<eps_den | ~isfinite(psid);
        psid(z) = eps_den;
    else
        error('psid에 0/NaN/Inf가 있습니다. (분모가 됨)');
    end
end

% ---------- CVX 시작 ----------
cvx_clear
cvx_begin quiet
    variables eps0(nK,nK) epsd(nK,nK)
    expression T

    T = 0;
    K_TV_minus_K1 = setdiff(K_TV, K1);

    for ii = K_TV_minus_K1
        ci = c(ii);  di = d(ii);
        for jj = 1:nK
            if jj == ii, continue; end
            T = T ...
              + eps0(ii,jj) * ( ci/rho_max(jj) + 2*di/psi0(ii,jj) ) ...
              + epsd(ii,jj) * ( ci/rho_max(jj) + 2*di/psid(ii,jj) );
        end
        T = T + eps0(ii,ii) * ( ci / rho_max(ii) );
    end

    for ii = K1
        T = T + c(ii)/rho_max(ii);
    end

    minimize( T )
    subject to
    eps0 >= 0; eps0 <= 1;
    epsd >= 0; epsd <= 1;

    % (1)
    for ii = K_TV
        ci = c(ii);
        eps0(ii,ii) * ( ci/rho_max(ii) ) <= t_tol(ii);
    end

    % (2)
    if ~isempty(j0)
        for ii = K_TV
            ci = c(ii);  di = d(ii);
            ub = min( t_tol(ii), t_hold0(ii) );
            eps0(ii,j0) * ( ci/rho_max(j0) + 2*di/psi0(ii,j0) ) <= ub;
        end
    end

    % (3)(4)
    for ii = K_TV
        ci = c(ii);  di = d(ii);
        for jj = 1:nK
            if jj == ii, continue; end
            ub_d = min( t_tol(ii), t_holdd(ii) );
            epsd(ii,jj) * ( ci/rho_max(jj) + 2*di/psid(ii,jj) ) <= ub_d;

            ub_0 = min( t_tol(ii), t_holdij(ii,jj) );
            eps0(ii,jj) * ( ci/rho_max(jj) + 2*di/psi0(ii,jj) ) <= ub_0;
        end
    end

    % (6)
    for jj = K_SV
        sum( eps0(K_TV, jj) ) == 1;
        sum( epsd(K_TV, jj) ) == 1;
    end

    % (7)
    for ii = K_TV
        for jj = K_SV
            eps0(ii,jj) + epsd(ii,jj) == 1;
        end
    end

    % (8)
    for ii = K_TV
        for jj = 1:nK
            if jj == ii, continue; end
            gamma(ii,jj) >= gamma_VUE * eps0(ii,jj);
            gamma(ii,jj) >= gamma_VUE * epsd(ii,jj);
        end
    end
cvx_end

% -------- one-hot 라운딩 --------
eps0_val = full(eps0);
epsd_val = full(epsd);

eps0_onehot = zeros(nK,nK);
epsd_onehot = zeros(nK,nK);
choice = repmat(struct('j',[],'k',[],'val',[]), nK, 1);

for ii = K_TV
    [v0, idx0] = max(eps0_val(ii, K_SV));
    [vd, idxd] = max(epsd_val(ii, K_SV));

    if vd > v0
        jj = K_SV(idxd);
        epsd_onehot(ii, jj) = 1;
        choice(ii).j = jj; choice(ii).k = 'dagger'; choice(ii).val = vd;
    else
        jj = K_SV(idx0);
        eps0_onehot(ii, jj) = 1;
        choice(ii).j = jj; choice(ii).k = '0';      choice(ii).val = v0;
    end
end

% -------- 결과 --------
out.status        = cvx_status;
out.T_value       = full(T);
out.eps0          = eps0_val;
out.epsd          = epsd_val;
out.eps0_onehot   = eps0_onehot;
out.epsd_onehot   = epsd_onehot;
out.choice        = choice;
end
