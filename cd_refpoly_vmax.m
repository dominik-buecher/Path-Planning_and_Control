function [c, te] = cd_refpoly_vmax(vmax, x0, xs)
% CD_REFPOLY_VMAX Computes coefficients for a reference polynomial with specified maximum velocity

% Calculate time to reach maximum velocity
te = (15 * xs) / (8 * vmax);

% Calculate coefficients of the reference polynomial
c5 = 6 * xs / (te.^5);
c4 = -15 * xs / (te.^4);
c3 = 10 * xs / (te.^3);
c2 = 0;
c1 = 0;
c0 = x0;

% Combine coefficients into an output vector
c = [c5, c4, c3, c2, c1, c0];
end
