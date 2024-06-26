%
% Mini-Auto-Drive
%
% System Parameters
%
% Copyright (c) 2020 Frank Traenkle
% http://www.modbas.de
%
% This file is part of Mini-Auto-Drive.
%
% Mini-Auto-Drive is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Mini-Auto-Drive is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Mini-Auto-Drive.  If not, see <http://www.gnu.org/licenses/>.

%% Clear workspace
clear all;

%% Global Variables
global P_dt P_display_dt;

%% Car and Belt Configuration
p_mad_car;

%% Sample Times
P_dt = 20e-3; % sample time of controller
P_sim_dt = 2e-3; % sample time of simulation
P_display_dt = 40e-3; % sample time of display

%% EKF e6
%e6_data;

%% Speed Controller
P_kr= 0.3440; % [s/m]
P_Ti = 0.2468; % [s] 
P_Ta = 20e-3; % [s] 

%% Lateral Controller
P_p_Tw = 300e-3; % [s]
P_k_p = 1; % [1/s] 

P_maneuvertype = true; % boolean
P_xManeuverEnd = 1; % [m]
P_ManeuverTypePathFollow = 0; 
P_ManeuverTypePark = 1; 

%% Safety Halt
P_safety_buffer = 0.05; % [m]
P_prevent_buffer = 0.1; % [m]
P_steering_extra = 0.1; % [m]

%% Create Race Track
% a1total = 2.7; % total surface width [ m ]
% a2total = 1.8; % total surface height [ m ]
% a1boundary = 0.05; % margin [ m ]
% a2boundary = 0.05; % margin [ m ]
% a1 = a1total - 2 * a1boundary; % total surface width [ m ]
% a2 = a2total - 2 * a2boundary; % total surface height [ m ]
% P_width = 0.25 * a2; % track P_width [ m ]
% 
% track = mbc_track_create(a1boundary + P_width, a2boundary + 0.5 * P_width, 0);
% track = mbc_straight_create(track, a1 - 2 * P_width, P_width);
% track = mbc_circle_create(track, 0.5 * P_width, pi, P_width);
% track = mbc_straight_create(track, a1 - 3 * P_width, P_width);
% track = mbc_circle_create(track, 0.5 * P_width, -pi, P_width);
% track = mbc_straight_create(track, a1 - 3 * P_width, P_width);
% track = mbc_circle_create(track, 0.5 * P_width, pi, P_width);
% track = mbc_straight_create(track, a1 - 2 * P_width, P_width);
% track = mbc_circle_create(track, 0.5 * P_width, 0.5 * pi, P_width);
% track = mbc_straight_create(track, a2 - 2 * P_width, P_width);
% track = mbc_circle_create(track, 0.5 * P_width, 0.5 * pi, P_width);
% track = mbc_track_display(track, 0.1, [ 0 a1total 0 a2total ]);
% path = track.center;
%spline = mbc_spline_create(track, 0.1, 0);

%% New Road Oval
% Road Surface
a1total = 2.70; % total surface width [ m ]
a2total = 1.8; % total surface height [ m ]
width = 0.2; % track width [ m ]

% Create Oval
track = mbc_track_create(0.150, 0.9, -pi/2);
track = mbc_straight_create(track, 0.22, width);
track = mbc_clothoid_create(track, 8, pi/4, width, 0);
track = mbc_clothoid_create(track, 8, pi/4, width, 1);
track = mbc_straight_create(track, 1.345, width);
track = mbc_clothoid_create(track, 8, pi/4, width, 0);
track = mbc_clothoid_create(track, 8, pi/4, width, 1);
track = mbc_straight_create(track, 0.443, width);
track = mbc_clothoid_create(track, 8, pi/4, width, 0);
track = mbc_clothoid_create(track, 8, pi/4, width, 1);
track = mbc_straight_create(track, 1.345, width);
track = mbc_clothoid_create(track, 8, pi/4, width, 0);
track = mbc_clothoid_create(track, 8, pi/4, width, 1);
track = mbc_straight_create(track, 0.223, width);
track = mbc_track_display(track, 0.1, [ 0 a1total 0 a2total ]);
path = track.center;

% a1total = 2.91; % total surface width [ m ]
% a2total = 2.31; % total surface height [ m ]
% width = 0.6; % track width [ m ]
% a1boundary = 0.05; % margin [ m ]
% a2boundary = 0.05; % margin [ m ]
% a1 = a1total - 2 * a1boundary; % total surface width [ m ]
% a2 = a2total - 2 * a2boundary; % total surface height [ m ]
% P_width = 0.45; % track P_width [ m ]
% % Create Oval
% track = mbc_track_create(0.330, 1.16, -pi/2);
% track = mbc_straight_create(track, 0.3, P_width);
% track = mbc_clothoid_create(track, 8, pi/4, P_width, 0);
% track = mbc_clothoid_create(track, 8, pi/4, P_width, 1);
% track = mbc_straight_create(track, 1.2, P_width);
% track = mbc_clothoid_create(track, 8, pi/4, P_width, 0);
% track = mbc_clothoid_create(track, 8, pi/4, P_width, 1);
% track = mbc_straight_create(track, 0.6, P_width);
% track = mbc_clothoid_create(track, 8, pi/4, P_width, 0);
% track = mbc_clothoid_create(track, 8, pi/4, P_width, 1);
% track = mbc_straight_create(track, 1.2, P_width);
% track = mbc_clothoid_create(track, 8, pi/4, P_width, 0);
% track = mbc_clothoid_create(track, 8, pi/4, P_width, 1);
% track = mbc_straight_create(track, 0.300, P_width);
% track = mbc_track_display(track, 0.1, [ 0 a1total 0 a2total ]);
% path = track.center;


%% New Road Straight
% Road Surface
% a1total = 6.0; % total surface width [ m ]
% a2total = 0.8; % total surface height [ m ]
% width = 0.6; % track width [ m ]
% 
% % Create Oval
% track = mbc_track_create(0.0, 0.4, 0);
% track = mbc_straight_create(track, 5.0, width);
% track = mbc_track_display(track, 0.1, [ 0 a1total 0 a2total ]);
% path = track.center;


%% Path for Lap Statistics
lappath = track.center;
P_lap_breakslen = uint32(length(lappath.points));
P_lap_points = zeros(SPLINE.Elements(2).Dimensions); 
P_lap_points(:,1:length(lappath.points)) = lappath.points;
P_lap_coefs = zeros(SPLINE.Elements(3).Dimensions);
P_lap_coefs(1:length(lappath.pp.coefs),:) = lappath.pp.coefs;
P_lap_segments = uint32(zeros(SPLINE.Elements(4).Dimensions)); 

%% Workspace variables for reference track generation in Simulink
P_w_breakslen = uint32(length(path.points));
P_w_points = zeros(SPLINE.Elements(2).Dimensions); 
P_w_points(:,1:length(path.points)) = path.points;
P_w_coefs = zeros(SPLINE.Elements(3).Dimensions);
P_w_coefs(1:length(path.pp.coefs),:) = path.pp.coefs;
P_w_segments = uint32(zeros(SPLINE.Elements(4).Dimensions)); 

% Init car display
mbc_car_display(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
