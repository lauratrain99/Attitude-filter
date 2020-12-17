function pos_n = pos_update(pos, vel, dt)
% pos_update: updates position in the navigation frame (lat, lon, h).
%
% INPUT:
%   pos,    3x1 position vector [lat lon h] (rad, rad, m).
%   vel,    3x1 NED velocities (m/s).
%   dt,     sampling interval.
%
% OUTPUT:
%   pos,    3x1 updated position vector [lat lon h] (rad, rad, m).
%
%   Copyright (C) 2014, Rodrigo Gonzalez, all rights reserved. 
%     
%   This file is part of NaveGo, an open-source MATLAB toolbox for 
%   simulation of integrated navigation systems.
%     
%   NaveGo is free software: you can redistribute it and/or modify
%   it under the terms of the GNU Lesser General Public License (LGPL) 
%   version 3 as published by the Free Software Foundation.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU Lesser General Public License for more details.
% 
%   You should have received a copy of the GNU Lesser General Public 
%   License along with this program. If not, see 
%   <http://www.gnu.org/licenses/>.
%
% References: 
%
%	Titterton, D.H. and Weston, J.L. (2004). Strapdown
% Inertial Navigation Technology (2nd Ed.). Institution
% of Engineering and Technology, USA. Eq. 3.79-3.81.
%
%	R. Gonzalez, J. Giribet, and H. Patiño. NaveGo: a 
% simulation framework for low-cost integrated navigation systems, 
% Journal of Control Engineering and Applied Informatics, vol. 17, 
% issue 2, pp. 110-120, 2015. Eq. 18.
%
% Version: 006
% Date:    2019/01/14
% Author:  Rodrigo Gonzalez <rodralez@frm.utn.edu.ar>
% URL:     https://github.com/rodralez/navego 

lat = pos(1); 
lon = pos(2); 
h   = pos(3);
vn  = vel(1);
ve  = vel(2);
vd  = vel(3);

%% Altitude

h_n  = h - (vd) * dt;

if h_n < 0
    h_n = abs(h);
    warning('pos_update: altitude is negative.')
end

%% Latitude

[RM,~] = radius(lat);

vn_c = vn / (RM + h_n);

lat_n = lat + (vn_c) * dt;

%% Longitude

[~, RN] = radius(lat_n);

ve_c  = ve / ((RN + h_n) * cos (lat_n));

lon_n = lon + (ve_c) * dt;

%% Position update

pos_n = [lat_n lon_n h_n];

end
