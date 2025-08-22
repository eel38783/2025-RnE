
function psi2 = V_RIS(Nx, Ny, theta, UAV, veh_pos)
%UAV-Mounter RIS channel
    N = Nx*Ny;
    kappa = 40.05; %40.05dB
    lambda = 10^(-2); %10(mm)

    cos_phi_RIS = zeros(1,N);

    for i=1:N
        for j=1:15
            cos_phi_RIS(i,j) = (UAV(1)-veh_pos(j,1))/norm(UAV-veh_pos(j));
        end
    end
    g_RIS = zeros(15,N);
    %AoA direction cosine angle (ith vehicle to RIS)

    for i=1:N 
        for j=1:15
            g_RIS(i,j) = sqrt(kappa/(norm(veh_pos(j)-UAV))^2)*exp((-2*(i-1)*pi*cos_phi_RIS(i,j)*norm(veh_pos(j)-UAV))/lambda);
        end
    end
    g_RIS = conj(g_RIS); %Hermitian

    
    cos_phi_UAV = zeros(15,N);
    %AoA direction cosine angle (RIS to jth vehicle)
    for i=1:N 
        for j=1:15
            cos_phi_UAV(i,j) = (veh_pos(j,1)-UAV(1))/norm(UAV-veh_pos(j));
        end
    end

    g_UAV = zeros(15,N);

    for i=1:N 
        for j=1:15
            g_UAV(i,j) = sqrt(kappa/(norm(veh_pos(j)-UAV))^2)*exp((-2*(i-1)*pi*cos_phi_UAV(i,j)*norm(veh_pos(j)-UAV))/lambda);
        end
    end

    g_UAV = conj(g_UAV); %Hermitian

    arr_str = zeros(1,N);
    for j=1:N 
        arr_str(j) = exp(j*theta(j)); 
    end

    Theta = diag(theta);
    h = g_RIS + g_UAV * Theta;
    Pt = 20;
    N0 = 20;
    
    psi2 = zeros(15,N);
    for i=1:N 
        for j=1:15 
            psi2(i,j) = (abs(h(i,j))^(2) * Pt)/N0;
        end
    end
end