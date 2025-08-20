action_space = [-5:5] * (pi/10);
theta = action_space(randi(length(action_space)));


function [psi2] = V_RIS(Nx, Ny, theta, UAV, veh_pos)
%UAV-Mounter RIS channel
    N = Nx*Ny;
    kappa = 40.05; %40.05dB
    lambda = 10^(-2); %10(mm)

    cos_phi_RIS = zeros(1:N);

    for i=1:N 
        cos_phi_RIS(i) = (UAV(1)-veh_pos(i,1)/norm(UAV-veh_pos(i)));
    end
    g_RIS = zeros(1:N);
    %AoA direction cosine angle (ith vehicle to RIS)

    for i=1:N 
        g_RIS(i) = sqrt(kappa/distance(veh_pos(i),UAV)^2)*exp(-2*(i-1)*pi*cos_phi_RIS(i)*(i)/lambda);
    end
    g_RIS = g_RIS'; %Hermitian

    
    cos_phi_UAV = zeros(1:N);
    %AoA direction cosine angle (RIS to jth vehicle)
    for i=1:N 
        cos_phi_UAV(i) = (veh_pos(i,1)-UAV(1)/norm(UAV-veh_pos(i)));
    end

    g_UAV = zeros(1:N);

    for i=1:N 
        g_UAV(i) = sqrt(kappa/distance(veh_pos(i),UAV)^2)*exp(-2*(i-1)*pi*cos_phi_UAV(i)*(i)/lambda);
    end

    g_UAV = g_UAV'; %Hermitian

    arr_str = zeros(1:N);
    for i=1:N 
        arr_str(i) = exp(i*theta(i)); 
    end
    Theta = diag(theta);
    h = g_RIS + g_UAV*Theta;
    psi2 = h;
end
