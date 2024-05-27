function cff = cd_refpoly_ff(c, k, T, T_t, kr, T_i)
% CD_REFPOLY_FF Computes the feedforward coefficients for a reference polynomial

% Extract coefficients from input vector
c5 = c(1);
c4 = c(2);
c3 = c(3);
c2 = c(4);
c1 = c(5);
c0 = c(6);

% Calculate feedforward coefficients
cff5 = c5;
cff4 = (c4 + 5 * T_i * c5 * (1/(k*kr) + 1));
cff3 = (c3 + 4 * T_i * c4 * (1/(k*kr) + 1) + 20 * T_i * c5 * (T + T_t) / (k*kr));
cff2 = (c2 + (3 * T_i * c3 * (1/(k*kr) + 1)) + (12 * T_i * c4 * (T + T_t) / (k*kr)) + (60 * T * T_i * T_t * c5 / (k*kr)));
cff1 = ((c1 + 2 * T_i * c2 * (1/(k*kr) + 1)) + (6 * T_i * c3 * (T + T_t) / (k*kr)) + (24 * T * T_i * T_t* c4 / (k*kr)));
cff0 = (c0 + T_i * c1 * (1/(k*kr) + 1)) + (2 * T_i * c2 * (T + T_t) / (k*kr)) + (6 * T * T_i * T_t* c3 / (k*kr));

% Combine feedforward coefficients into an output vector
cff = [cff5, cff4, cff3, cff2, cff1, cff0];
end

