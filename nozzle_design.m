%% =========================================================================
%  DIVERGENT PART OF A CONVERGENT-DIVERGENT NOZZLE - METHOD OF CHARACTERISTICS
%  Minimum-length nozzle (sharp-edged throat), 2D planar, calorically perfect gas
%  -------------------------------------------------------------------------
%  IDEA OF THE METHOD
%  The supersonic flow is expanded through a sharp corner at the throat
%  (Prandtl-Meyer fan). The expansion waves travel to the centerline, where they
%  are reflected, and then travel towards the wall. The wall is shaped so that
%  these reflected waves are CANCELLED (no reflection) -> the flow becomes
%  uniform and parallel at the exit with the desired Mach number Me.
%
%  Only the upper half of the nozzle is computed (the flow is symmetric
%  with respect to the centerline y = 0).
%
%  CONVENTIONS (angles in radians inside the code)
%    theta = flow angle w.r.t. the x axis
%    nu    = Prandtl-Meyer function
%    mu    = Mach angle = asin(1/M)
%    C- characteristic : slope dy/dx = tan(theta - mu),  theta + nu = K-  (const)
%    C+ characteristic : slope dy/dx = tan(theta + mu),  theta - nu = K+  (const)
%
%  NAMING OF THE GRID POINTS
%    P(i,j) : point where the i-th C- line (coming from the throat corner)
%             meets the j-th C+ line (coming from the centerline).
%             j = 0 is the throat corner itself (start of the C- line),
%             j = i is the point on the centerline (end of the C- line).
%    In MATLAB the index j starts from 1, so the column used is  j+1.
%    (column 1 = corner, column 2 = first C+ line, ..., column i+1 = centerline)
%
%  Reference: J.D. Anderson, "Modern Compressible Flow", Chapter 11.
% =========================================================================
clear; clc; close all;

%% 1) INPUT DATA --------------------------------------------------------------
gamma = 1.4;      % ratio of specific heats [-]
Me    = 2.4;      % desired exit Mach number [-]
N     = 2;       % number of characteristic lines of the expansion fan
                  % (more lines = more accurate contour)
h     = 1.0;      % throat half-height [any length unit]

%% 2) MAXIMUM WALL ANGLE ------------------------------------------------------
% For a minimum-length nozzle the wall leaves the throat corner with the
% maximum expansion angle, which is half of the Prandtl-Meyer angle at the exit
nu_e      = PM(Me, gamma);      % Prandtl-Meyer function at the exit [rad]
theta_max = nu_e/2;             % maximum flow/wall angle [rad]

fprintf('Exit Prandtl-Meyer angle  nu_e     = %8.4f deg\n', rad2deg(nu_e));
fprintf('Maximum wall angle        theta_max = %8.4f deg\n\n', rad2deg(theta_max));

%% 3) PRE-ALLOCATION OF THE ARRAYS --------------------------------------------
% Rows    : i-th C- line (i = 1..N)
% Columns : j+1 (see the note at the top), so N+1 columns
x     = nan(N, N+1);   % x coordinate of the points
y     = nan(N, N+1);   % y coordinate of the points
theta = nan(N, N+1);   % flow angle
nu    = nan(N, N+1);   % Prandtl-Meyer function
M     = nan(N, N+1);   % Mach number
mu    = nan(N, N+1);   % Mach angle

%% 4) EXPANSION FAN AT THE THROAT CORNER (first column, j = 0) -----------------
% The corner is at (0,h). Each C- line leaves it with a different flow angle,
% going from theta_max/N (weakest wave) up to theta_max (last wave).
% Across a simple Prandtl-Meyer expansion starting from M = 1 (nu = 0,
% theta = 0):   nu = theta.
for i = 1:N
    theta(i,1) = i*theta_max/N;
    nu(i,1)    = theta(i,1);
    M(i,1)     = invPM(nu(i,1), gamma);   % Mach number from nu
    mu(i,1)    = asin(1/M(i,1));          % Mach angle
    x(i,1)     = 0;                       % all lines start at the corner
    y(i,1)     = h;
end

%% 5) INTERIOR POINTS AND CENTERLINE POINTS -----------------------------------
% We march along each C- line (index i). On the i-th C- line we compute
% the intersections with the C+ lines j = 1 ... i-1 and finally the point
% on the centerline (j = i).
for i = 1:N
    for j = 1:i
        c = j + 1;          % MATLAB column of the point being computed
        % Point B = previous point on the SAME C- line  -> column c-1
        % Point A = point on the C+ line j of the PREVIOUS C- line -> row i-1

        if j < i
            % ---------- INTERIOR POINT (unit process) ------------------------
            % Invariants carried by the two characteristics arriving at P:
            Kminus = theta(i,c-1)   + nu(i,c-1);     % from B (along C-)
            Kplus  = theta(i-1,c)   - nu(i-1,c);     % from A (along C+)
            % Flow properties in P (obtained by adding / subtracting):
            theta(i,c) = 0.5*(Kminus + Kplus);
            nu(i,c)    = 0.5*(Kminus - Kplus);
        else
            % ---------- POINT ON THE CENTERLINE ------------------------------
            % By symmetry the flow angle on the axis is zero, so the C-
            % invariant gives directly nu.
            theta(i,c) = 0;
            nu(i,c)    = theta(i,c-1) + nu(i,c-1);
        end

        % Mach number and Mach angle in P from nu
        M(i,c)  = invPM(nu(i,c), gamma);
        mu(i,c) = asin(1/M(i,c));

        % ---------- POSITION OF P ---------------------------------------------
        % Slope of the C- line B->P (average of the values in B and P)
        mB = tan(0.5*((theta(i,c-1) - mu(i,c-1)) + (theta(i,c) - mu(i,c))));

        if j < i
            % Slope of the C+ line A->P (average of the values in A and P)
            mA = tan(0.5*((theta(i-1,c) + mu(i-1,c)) + (theta(i,c) + mu(i,c))));
            % Intersection of the two straight lines:
            %   y = yB + mB*(x - xB)   and   y = yA + mA*(x - xA)
            x(i,c) = (y(i-1,c) - y(i,c-1) + mB*x(i,c-1) - mA*x(i-1,c)) / (mB - mA);
            y(i,c) = y(i,c-1) + mB*(x(i,c) - x(i,c-1));
        else
            % On the centerline y = 0: intersection of the C- line with y = 0
            y(i,c) = 0;
            x(i,c) = x(i,c-1) - y(i,c-1)/mB;
        end
    end
end

%% 6) WALL POINTS -------------------------------------------------------------
% Each C+ line j (j = 1..N) starting from the centerline reaches the wall.
% The wall is designed so that the wave is cancelled: the flow properties
% at the wall point W are the same as those of the last point A of that C+
% line (i.e. the point on the last C- line, row N). The wall point W_j
% is the intersection between:
%   - the C+ line coming from A (slope tan(theta_A + mu_A))
%   - the wall segment starting from the previous wall point, with slope equal
%     to the mean of the previous and the new wall angles.
xw     = nan(1, N+1);      % wall x coordinates
yw     = nan(1, N+1);      % wall y coordinates
thetaw = nan(1, N+1);      % wall angle
Mw     = nan(1, N+1);      % Mach number at the wall

% first wall point = throat corner, with maximum wall angle
xw(1) = 0;  yw(1) = h;  thetaw(1) = theta_max;  Mw(1) = 1;

for j = 1:N
    c = j + 1;                          % column of point A in row N
    thetaw(j+1) = theta(N,c);           % wall angle = flow angle in A
    Mw(j+1)     = M(N,c);               % same Mach number as A

    mC = tan(theta(N,c) + mu(N,c));                 % slope of the C+ line from A
    mW = tan(0.5*(thetaw(j) + thetaw(j+1)));        % mean wall slope

    % Intersection C+ line / wall segment
    xw(j+1) = (y(N,c) - yw(j) + mW*xw(j) - mC*x(N,c)) / (mW - mC);
    yw(j+1) = yw(j) + mW*(xw(j+1) - xw(j));
end

%% 7) CHECK WITH THE ISENTROPIC AREA RATIO ------------------------------------
% For 2D planar nozzle: A/A* = y_exit / h. It must match the isentropic
% area-Mach relation.
AR_isen = (1/Me) * ((2/(gamma+1))*(1 + (gamma-1)/2*Me^2))^((gamma+1)/(2*(gamma-1)));
AR_moc  = yw(end)/h;

fprintf('Exit half-height (MOC)      y_e/h = %8.4f\n', AR_moc);
fprintf('Isentropic area ratio       A/A*  = %8.4f\n', AR_isen);
fprintf('Relative error                    = %8.4f %%\n', 100*abs(AR_moc-AR_isen)/AR_isen);
fprintf('Nozzle length (divergent)   L/h   = %8.4f\n\n', xw(end)/h);

%% 8) PRINT / SAVE THE WALL CONTOUR -------------------------------------------
fprintf('   Point      x/h        y/h     theta[deg]     M\n');
for k = 1:N+1
    fprintf('  %4d   %9.4f  %9.4f  %9.4f  %9.4f\n', ...
            k-1, xw(k)/h, yw(k)/h, rad2deg(thetaw(k)), Mw(k));
end

fid = fopen('nozzle_contour.txt', 'w');          % simple text file: x  y
fprintf(fid, '%% x  y   (upper wall, same units as h)\n');
fprintf(fid, '%12.6f  %12.6f\n', [xw; yw]);
fclose(fid);
fprintf('\nContour saved in nozzle_contour.txt\n');

%% 9) PLOT 1: NOZZLE CONTOUR AND CHARACTERISTIC NET ---------------------------
figure('Name','Characteristic net','Color','w'); hold on; axis equal; box on;

for i = 1:N
    % C- line i : corner -> interior points -> centerline
    plot(x(i,1:i+1),  y(i,1:i+1),  'b-');
    plot(x(i,1:i+1), -y(i,1:i+1),  'b-');      % mirror image (lower half)
end
for j = 1:N
    % C+ line j : centerline point P(j,j) -> P(j+1,j) ... -> P(N,j) -> wall
    xr = [x(j,j+1), x(j+1:N,j+1)', xw(j+1)];
    yr = [y(j,j+1), y(j+1:N,j+1)', yw(j+1)];
    plot(xr,  yr, 'r-');
    plot(xr, -yr, 'r-');
end
plot(xw,  yw, 'k-', 'LineWidth', 2.5);          % upper wall
plot(xw, -yw, 'k-', 'LineWidth', 2.5);          % lower wall
plot([0 xw(end)], [0 0], 'k--');                % centerline
xlabel('x'); ylabel('y');
title(sprintf('Minimum-length nozzle, M_e = %.2f, \\gamma = %.2f, N = %d', Me, gamma, N));
% blue = C- characteristics, red = C+ characteristics, black = wall

%% 10) PLOT 2: MACH NUMBER IN THE GRID POINTS ---------------------------------
figure('Name','Mach number field','Color','w'); hold on; axis equal; box on;
scatter(x(:), y(:), 25, M(:), 'filled');
scatter(x(:), -y(:), 25, M(:), 'filled');
plot(xw,  yw, 'k-', 'LineWidth', 2);
plot(xw, -yw, 'k-', 'LineWidth', 2);
colorbar; xlabel('x'); ylabel('y'); title('Mach number at the characteristic points');

%% =========================================================================
%  LOCAL FUNCTIONS  (in a MATLAB script they must be at the end of the file)
% =========================================================================

function nu = PM(M, gamma)
% Prandtl-Meyer function nu(M) [rad], valid for M >= 1
    g  = (gamma+1)/(gamma-1);
    nu = sqrt(g)*atan(sqrt((M.^2 - 1)/g)) - atan(sqrt(M.^2 - 1));
end

function M = invPM(nu, gamma)
% Inverse Prandtl-Meyer function: finds M such that PM(M) = nu.
% Solved numerically with fzero in the interval M in [1, 1000].
    if nu <= 0
        M = 1;                                  % sonic condition
        return
    end
    f = @(Mach) PM(Mach, gamma) - nu;
    M = fzero(f, [1, 1000]);
end