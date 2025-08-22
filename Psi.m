
function psi = Psi(veh_pos, psi1, psi2)   % SNR로 Transmission Rate 계산
    psi = zeros(15, 15, 2);
    
    for i=1:15 
        d0 = norm(veh_pos(i,1:3) - [0,0,5]); %RSU Height <-- 10m
        psi0 = V_RSU(d0);
        psi(i,1,1) = log2(1 + psi0); 
        psi(i,1,2) = log2(1 + psi0); 
        for j=2:15
            psi(i,j,1) = log2(1 + psi1(i,j));
            psi(i,j,2) = log2(1 + psi2(i,j));
        end
    end 
                  
end 
