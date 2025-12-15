function [r,a,b,mu,nu,poles,residues] = ratdisk(f,m,n,N,tol)
% RATDISK  Rational approximation on the unit circle.
%
% Input:  Function f or vector of data at z_j = exp(2i*pi*(0:N)/(N+1))
%         for some N >= m+n. If N >> m+n, it is best to choose N odd.
%         Maximal numerator, denominator degrees m,n.
%
%         An optional 5th argument specifies relative tolerance tol.
%         If omitted, tol = 1e-14. Use tol = 0 to turn off robustness.
%
% Output: function handle r of exact type (mu,nu) approximant to f
%         with coefficient vectors a and b and optional poles,residues.
%
% (Structure based on Gonnet–Pachon–Trefethen, cleaned for MATLAB.)

% ---------------------------------------------------------------
% 1. Set N if not given
% ---------------------------------------------------------------
if nargin < 4
    if isnumeric(f)
        N = length(f) - 1;      % do interpolation if f is data vector
    else
        N = m + n;              % minimal N if f is function handle
    end
end
N1 = N + 1;                      % no. of roots of unity

% ---------------------------------------------------------------
% 2. Default tolerance
% ---------------------------------------------------------------
if nargin < 5
    tol = 1e-14;                 % default relative tolerance
end

% ---------------------------------------------------------------
% 3. Sample values fj either from data vector or function handle
% ---------------------------------------------------------------
if isnumeric(f)
    fj = f(:);                   % allow for either function or data
else
    zj = exp(2i*pi*(0:N).' / N1);
    fj = f(zj);                  % handle or data vector
end

ts = tol * norm(fj,inf);         % absolute tolerance
M  = floor(N/2);                 % no. of points in upper half-plane

f1 = fj(2:M+1);                  % fj in upper half-plane
f2 = fj(N+2-M:N1);               % fj in lower half-plane

realf = norm(f1(M:-1:1) - conj(f2), inf) < ts;  % true if fj is real symmetric
oddN  = mod(N,2) == 1;                           % true if N is odd
evenf = oddN && (norm(f1 - f2, inf) < ts);      % true if fj is even
oddf  = oddN && (norm(f1 + f2, inf) < ts);      % true if fj is odd

% ---------------------------------------------------------------
% 4. First row/column of Toeplitz matrix
% ---------------------------------------------------------------
row = conj(fft(conj(fj))) / N1;  % 1st row of Toeplitz matrix
col = fft(fj) / N1;              % 1st column of Toeplitz matrix
col(1) = row(1);

if realf
    row = real(row);             % discard negligible imaginary parts
    col = real(col);
end

d = xor(evenf, mod(m,2) == 1);   % either 0 or 1

% ---------------------------------------------------------------
% 5. Main stabilization loop
% ---------------------------------------------------------------
while true
    % Toeplitz matrix
    Z = toeplitz(col, row(1:n+1));

    if ~oddf && ~evenf            % fj is neither even nor odd
        [~,S,V] = svd(Z(m+2:N1,:), 0);   % SVD
        b = V(:, n+1);                    % coeffs of q
    else                           % fj is even or odd
        [~,S,V] = svd(Z(m+2+d:2:N1, 1:2:n+1), 0);  % symmetry treatment
        b = zeros(n+1,1);
        b(1:2:end) = V(:, end);    % coeffs of q
    end

    % smallest singular value (for robustness)
    if N > m+n && n > 0
        ssv = S(end,end);          % smallest singular value
    else
        ssv = 0;                   % or 0 in case of interpolation
    end

    % values of q at z_j
    qj = ifft(b, N1);
    qj = qj(:);

    % coefficients of p-hat
    ah = fft(qj .* fj);
    a  = ah(1:m+1);               % coeffs of p

    if realf, a = real(a); end    % discard imaginary rounding errors
    if evenf, a(2:2:end) = 0; end % enforce even symmetry of coeffs
    if oddf,  a(1:2:end) = 0; end % enforce odd symmetry of coeffs

    % --- stabilization by discarding singular values ---
    if tol > 0                    % tol=0 means no stabilization
        ns = n;                   % number of singular values
        if oddf || evenf
            ns = floor(n/2);
        end
        s  = diag(S(1:ns,1:ns));  % extract singular values
        nz = sum(s - ssv <= ts);  % number of singular values to discard
        if nz == 0
            break                 % if no discards, we are done
        else
            n = n - nz;           % reduce denominator degree
        end
    else
        break                     % no iteration if tol=0
    end
end  % end main loop

% ---------------------------------------------------------------
% 6. Remove tiny trailing coefficients
% ---------------------------------------------------------------
nna = abs(a) > ts;                % nonnegligible a coeffs
nnb = abs(b) > tol;               % nonnegligible b coeffs
kk  = 1:min(m+1,n+1);             %#ok<NASGU> indices a and b have in common

if any(nna)
    a = a(1:find(nna,1,'last'));  % discard trailing zeros of a
else
    a = 0;
end

if any(nnb)
    b = b(1:find(nnb,1,'last'));  % discard trailing zeros of b
else
    b = 1;
end

% special case of zero function
if isempty(a)
    a = 0;
    b = 1;
end

mu = length(a) - 1;               % exact numerator degree
nu = length(b) - 1;               % exact denominator degree

% function handle for r
r = @(z) polyval(a(end:-1:1), z) ./ polyval(b(end:-1:1), z);

% ---------------------------------------------------------------
% 7. Poles and residues (optional)
% ---------------------------------------------------------------
if nargout > 5
    poles = roots(b(end:-1:1));           % poles
    t = max(tol, 1e-7);                   % perturbation for residue estimate
    residues = t * (r(poles + t) - r(poles - t)) / 2;  % residue estimates
else
    poles = [];
    residues = [];
end
end
