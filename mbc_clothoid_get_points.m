function points = mbc_clothoid_get_points(track, idx, xstart, dx, alpha)
% This function calculates the points on a clothoid curve.
% track:   The track structure containing the clothoid.
% idx:     The index of the current clothoid in the track structure.
% xstart:  The starting x-value along the clothoid.
% dx:      The distance between points to be calculated.
% alpha:   A parameter to adjust the position along the width of the track.

    % Default alpha to 0.5 (center line) if not enough input arguments
    if nargin < 5
        alpha = 0.5; % return center line points
    end

    % Extracting point and track information from the track structure
    p = track.points{idx};
    t = track.tracks{idx};

    % Assigning clothoid parameters from the track
    opening = t.opening; % Whether the clothoid is opening (1) or closing (0)
    a = t.a; % The parameter 'a' of the clothoid
    xe = t.xe; % Endpoint x-coordinate
    w = t.w; % Width of the track

    % Extracting start position and orientation
    s1_p = p.s1; % Start position on s1-axis
    s2_p = p.s2; % Start position on s2-axis
    x_p = p.x; % Start x-coordinate
    psi_p = p.psi; % Start orientation

    % Calculating indices and x-coordinates for points
    index = 0:floor((abs(xe) - (xstart - x_p) + mbc_cmp_eps)/dx);
    x = (xstart - x_p) + dx * index;

    % Defining functions for s1 and s2 calculations based on opening/closing
    if opening == 1
        % For opening clothoid
        s1_func = @(x) cos(-a/2 * x.^2 + a*xe.*x + psi_p);
        s2_func = @(x) sin(-a/2 * x.^2 + a*xe.*x + psi_p);
        psi = -a/2 .* (xe - x).^2 + psi_p + a/2 * xe^2;
    else
        % For closing clothoid
        s1_func = @(x) cos(a/2 * x.^2 + psi_p);
        s2_func = @(x) sin(a/2 * x.^2 + psi_p);
        psi = a/2 .* x.^2 + psi_p;
    end

    % Initializing arrays for integrated s1 and s2 values
    s1_integ = zeros(size(x));
    s2_integ = zeros(size(x));

    % Calculating integrated s1 and s2 values for each point
    for i = 1:length(x)
        s1_integ(i) = integral(s1_func, 0, x(i));
        s2_integ(i) = integral(s2_func, 0, x(i));
    end

    % Calculating final s1 and s2 positions for each point
    s1 = s1_integ + s1_p + (alpha - 0.5) * w * cos(psi + pi/2);
    s2 = s2_integ + s2_p + (alpha - 0.5) * w * sin(psi + pi/2);

    % Compiling the calculated points into a matrix
    points = [p.x+x; s1; s2];
end
