function track = mbc_clothoid_create(track, a, rad, w, opening)
    % This function creates a clothoid segment and adds it to the track structure.
    % track:   The existing track structure to which the new segment will be added.
    % a:       The parameter 'a' of the clothoid, defining its curvature.
    % rad:     The rotation angle or radial change for the clothoid.
    % w:       The width of the clothoid track.
    % opening: Specifies if the clothoid is opening (1) or closing (0).

    % Get the current count of segments in the track
    cnt = mbc_track_get_cnt(track);

    % Extract the last point in the track
    p = track.points{cnt+1};

    % Ensure 'a' is positive and set starting parameters
    a = abs(a);
    s_01 = p.s1;  % Start position on s1-axis
    s_02 = p.s2;  % Start position on s2-axis
    psi_0 = p.psi; % Start orientation

    % Calculate endpoint x-coordinate based on the radius
    xe = sqrt(2 / a * abs(rad));
    
    % Adjust 'a' for negative radius (inverting curvature)
    if rad < 0 
        a = -a;
    end

    % Prepare integrand functions for calculating the endpoint positions
    integrandS1_SK = @(s) cos(a/2 * s.^2 + psi_0);
    integrandS2_SK = @(s) sin(a/2 * s.^2 + psi_0);
    integrandS1_OK = @(s, x) cos(-a/2 * s.^2 + a * x * s + psi_0);
    integrandS2_OK = @(s, x) sin(-a/2 * s.^2 + a * x * s + psi_0);

    % Calculate the endpoint positions based on whether the clothoid is opening or closing
    if opening == 1
        % For opening clothoid
        S1_xe = s_01 + integral(@(s) integrandS1_OK(s, xe), 0, xe);
        S2_xe = s_02 + integral(@(s) integrandS2_OK(s, xe), 0, xe);
    else 
        % For closing clothoid
        S1_xe = s_01 + integral(@(s) integrandS1_SK(s), 0, xe);
        S2_xe = s_02 + integral(@(s) integrandS2_SK(s), 0, xe);
    end

    % Update the track structure with the new endpoint
    track.points{cnt+2} = ...
        struct('s1', S1_xe, ...
               's2', S2_xe, ...
               'psi', psi_0 + rad, ...
               'x', p.x + abs(xe)); 

    % Update the track structure with the new clothoid segment details
    track.tracks{cnt+1} = struct('type', 'clothoid', 'a', a, 'xe', xe, 'w', w, 'opening', opening);
end
