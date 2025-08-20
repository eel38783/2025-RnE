function T = opt_func( K, epsilon, rho, Ds, C)
% Ds (task data saize list), C (required computation list)

veh_pos = [rand(20,2)*200, 0]; % vechile location
psi = SNR(V2V, V2I, V2UAV, veh_pos);
UAV = [10,10,20]; % fix UAV Height as 20m
psi = Rate(veh_pos); % transmission rate

T = 0;
for i=1:(K(2)+K(3))
    for j=1:(K(2)+K(3)+K(4))
        if j==i; continue; end
            for k=1:2 
                psi = Rate(veh_pos);
                T = T + epsilon(i,j,k)*(C(i)/rho(j) ...
                    + (2*Ds(i))/psi(i,j,k))+epsilon(i,i,1)*(C(i)/rho(j));
            end
    end
end

for i=1:K(1) 
    T = T + C(i)/rho(i);
end

end
K = [3,6,6,3]; %vehicular set
rho = randi([2,8], length(K), 1); %allocated resource
Ds = randi([2,10], length(K), 1); %amount of data
C = randi([2,10], length(K), 1);

cvx_begin
    variable epsilon
    minimize Opt_func(K, epsilon, rho, Ds, C, psi)
    subject to
        for i=1:(K(2)+K(3)) 
            for j=1:(K(2)+K(3)+K(4))
                for k=1:2 
                    eps(i,j,k) >= 0;
                    eps(i,j,k) <= 1;
                end
            end
        end
        % relaxation factor
        for i=1:K(1)
            psi = Rate(veh_pos);
            epsilon(i,1,1)*(C(i)/rho(j)+(2*Ds(i)/psi(i,1,1))) ...
                   <= min(t_tole(i),t_hold(i,1,1));
        end
        
        for i=1:(K(2)+K(3)) 
            for j=2:(K(2)+K(3)+K(4)) 
                epsilon(i,j,1)*(C(i)/rho(j)+(2*Ds(i)/psi(i,j,1))) ...
                    <= min(t_tole(i), t_hold(i,j,1));
            end
        end

        for i=1:(K(2)+K(3))
            for j=2:(K(2)+K(3)+K(4)) 
                epsilon(i,j,2)*(C(i)/rho(j)+(2*Ds(i)/psi(i,j,2))) ...
                    <= min(t_tole(i), t_hold(i,j,2));
            end
        end
        % time constraints
        e = 0;
        for j=1:(K(2)+K(3)+K(4)) 
            for k=1:2 
                for i=(K(1)+K(2)+K(3)) 
                    e = e + epsilon(i,j,k);
                end
                e == 1;
            end
        end
        % to ensure 1:1 match

        for i=1:(K(2)+K(3)) 
            for j=2:(K(2)+K(3)+K(4)) 
                epsilon(i,j,1) + eps(i,j,2) == 1; 
            end 
        end
        % to ensure V2V Mode is either on the ground or through UAV-RIS
cvx_end