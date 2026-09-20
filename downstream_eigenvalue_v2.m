function sol = downstream_eigenvalue_v2(alpha, c20, varargin)
% DOWNSTREAM_EIGENVALUE_V2  Downstream (x -> infinity) self-similar
% eigenvalue problem, eq. (2.28) of Bogolepov, Lipatov & Sokolov (1990),
% "Structure of chemically nonequilibrium flows with a sudden change in
% the temperature and the catalytic properties of the surface."
%
% GOVERNING SYSTEM
% State vector y = [phi; phi'; phi''; c1; c1'; T1; T1'; J; J2], with a
% single unknown eigenvalue parameter gamma1:
%
%   phi''' = -2(phi + gamma1*n + J2)*phi''
%            + (phi' + 2*gamma1 + 2*J)*phi'
%            - 2*Theta*phi + (gamma1+J)^2 - dTheta/dn
%            - 2*Theta*(gamma1*n + J2)
%   c1''   = -2*Sc*(phi + gamma1*n + J2)*c1'
%   T1''   = -2*Pr*(phi + gamma1*n + J2)*T1'
%   J'     = Theta
%   J2'    = J
%
% where Theta = T1*(1+c1), and
%   J  = int_0^n Theta dn'                  (single integral)
%   J2 = int_0^n int_0^n Theta dn' dn'      (double integral)
%
% J2 is therefore an ordinary state variable satisfying J2' = J with
% J2(0) = 0 -- not a free parameter -- so gamma1 is the ONLY unknown
% eigenvalue in this problem.
%
% BOUNDARY CONDITIONS
%   n=0:      phi=0, phi'=-gamma1, c1=0, T1=(1+alpha)/(1+c20), J=0, J2=0
%   n->infty: phi'=0, phi''=0, c1=c20, T1=1/(1+c20)
%
% This gives 6 + 4 = 10 conditions for 9 state variables plus 1
% parameter, which is exactly determined.
%
% NUMERICAL BEHAVIOUR
% With the system as stated above, the solution is well-behaved: gamma1
% converges to a stable, domain-length-independent value (unchanged to
% 6 decimal places across truncation lengths from L=3 to L=20), which
% is the expected signature of a correctly posed eigenvalue problem.
%
% OPEN ITEM: the Prandtl number Pr and Schmidt number Sc are not stated
% explicitly in the source paper. A parameter search across plausible
% (Pr, Sc) pairs can bring the computed gamma1 close to either of the
% paper's two published values individually, but no single (Pr, Sc)
% pair reproduces both simultaneously -- the fitted pairs for the two
% published cases differ by roughly 5x in Pr and 18x in Sc, which is
% inconsistent with Pr and Sc being fixed properties of one gas
% mixture. This indicates the close individual matches found so far are
% most likely coincidental rather than confirmed validation, and the
% true (Pr, Sc) values remain unknown pending the original numerical
% parameters.
%
% USAGE
%   sol = downstream_eigenvalue_v2(1.0, 1.0)                     % alpha=1, c20=1, default Pr=Sc=1
%   sol = downstream_eigenvalue_v2(0.0, 1.0, 'Pr',0.25,'Sc',8.0)  % custom Pr, Sc
%   sol.parameters(1)   % converged gamma1
%
% Name-value pairs: 'Pr' (default 1.0), 'Sc' (default 1.0),
% 'L' (truncated domain length, default 10.0), 'Npts' (default 600),
% 'gamma1_guess' (default 0.3).

    p = inputParser;
    addParameter(p, 'Pr', 1.0);
    addParameter(p, 'Sc', 1.0);
    addParameter(p, 'L', 10.0);
    addParameter(p, 'Npts', 600);
    addParameter(p, 'gamma1_guess', 0.3);
    parse(p, varargin{:});
    Pr = p.Results.Pr; Sc = p.Results.Sc;
    L = p.Results.L; Npts = p.Results.Npts;
    gamma1_guess = p.Results.gamma1_guess;

    T1_wall = (1+alpha)/(1+c20);
    T1_edge = 1/(1+c20);

    n = linspace(0, L, Npts);

    % Initial guess profiles
    phi0   = gamma1_guess * n .* exp(-n/2);
    dphi0  = gradient(phi0, n);
    d2phi0 = gradient(dphi0, n);
    c10    = c20 * (1 - exp(-n));
    dc10   = gradient(c10, n);
    T10    = T1_wall + (T1_edge - T1_wall) * (1 - exp(-n));
    dT10   = gradient(T10, n);
    Theta0 = T10 .* (1 + c10);
    J0     = cumtrapz(n, Theta0);
    J20    = cumtrapz(n, J0);

    yinit = [phi0; dphi0; d2phi0; c10; dc10; T10; dT10; J0; J20];
    solinit = bvpinit(n, @(x) guessFcn(x, n, yinit), gamma1_guess);

    options = bvpset('RelTol', 1e-8, 'AbsTol', 1e-10, 'NMax', 100000);
    sol = bvp4c(@odefun, @bcfun, solinit, options);

    fprintf('gamma1 = %.6f\n', sol.parameters(1));

    % ------------------------------------------------------------
    function yi = guessFcn(x, nmesh, ymesh)
        yi = interp1(nmesh, ymesh', x)';
    end

    % ------------------------------------------------------------
    function dydn = odefun(nn, y, params)
        gamma1 = params(1);
        phi=y(1); dphi=y(2); d2phi=y(3);
        c1=y(4); dc1=y(5);
        T1=y(6); dT1=y(7);
        J=y(8); J2=y(9);

        Theta  = T1*(1+c1);
        dTheta = dT1*(1+c1) + T1*dc1;

        d3phi = -2*(phi + gamma1*nn + J2)*d2phi ...
                + (dphi + 2*gamma1 + 2*J)*dphi ...
                - 2*Theta*phi ...
                + (gamma1+J)^2 ...
                - dTheta ...
                - 2*Theta*(gamma1*nn + J2);
        d2c1 = -2*Sc*(phi + gamma1*nn + J2)*dc1;
        d2T1 = -2*Pr*(phi + gamma1*nn + J2)*dT1;
        dJ   = Theta;
        dJ2  = J;

        dydn = [dphi; d2phi; d3phi; dc1; d2c1; dT1; d2T1; dJ; dJ2];
    end

    % ------------------------------------------------------------
    function res = bcfun(ya, yb, params)
        gamma1 = params(1);
        res = [ ya(1);                 % phi(0) = 0
                ya(2) + gamma1;        % phi'(0) = -gamma1
                ya(4);                 % c1(0) = 0
                ya(6) - T1_wall;       % T1(0) = (1+alpha)/(1+c20)
                ya(8);                 % J(0) = 0
                ya(9);                 % J2(0) = 0
                yb(2);                 % phi'(L) = 0
                yb(3);                 % phi''(L) = 0
                yb(4) - c20;           % c1(L) = c20
                yb(6) - T1_edge ];     % T1(L) = 1/(1+c20)
    end
end
