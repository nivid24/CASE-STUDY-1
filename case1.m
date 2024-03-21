% Data
Lx = 25; 
Ly = 20;
Nx = 75; 
Ny = 75; 
dx = Lx / (Nx - 1);
dy = Ly / (Ny - 1);
tol = 1e-6; % Tolerance for convergence
maxIter = 100000; % Maximum number of iterations

% Initialization
T = zeros(Nx, Ny); %Temperature matrix
T(:, 1) = 100; % Dirichlet 
T(:, end) = 0; % Dirichlet
T(1, :) = T(2, :); % Neumann
T(end, :) = T(end-1, :); %Neumann
Tnew = T; % New temperature matrix


for iter = 1:maxIter
    T = Tnew;
    for i = 2:Nx-1
        for j = 2:Ny-1
            Tnew(i, j) = (T(i+1, j) + T(i-1, j) + T(i, j+1) + T(i, j-1)) / 4; 
        end
    end
    % Neumann boundary conditions
    Tnew(:, end) = Tnew(:, end-1); 
    Tnew(end, :) = Tnew(end-1, :); 
    
    %convergence point
    if max(abs(Tnew(:) - T(:))) < tol
        disp(['Converged at iteration ', num2str(iter)]);
        break;
    end
end

% Plot
[X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));
contourf(X, Y, Tnew', 75 , 'LineColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
title('D-N temp distribution');