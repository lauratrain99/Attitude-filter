function ecef = ned2ecef(ned, llh_org)
% ecef2ned: converts from ECEF coordinates to NED coordinates.
%
% INPUTS
%   ned: Nx3 NED coordinates [X Y Z] (m, m, m)
%   llh_org: 1x3 system origin [lat, lon, h] (rad, rad, m).
%
% OUTPUTS
%   ecef: Nx3 ECEF coordinates [X Y Z] (m, m, m).
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
% Reference:
%
%
% Version: 002
% Date:    2019/01/16
% Author:  Rodrigo Gonzalez <rodralez@frm.utn.edu.ar>
% URL:     https://github.com/rodralez/navego

lat = llh_org(1);
lon = llh_org(2);

ecef_org = llh2ecef(llh_org)';

slat = sin(lat);
clat = cos(lat);
slon = sin(lon);
clon = cos(lon);

R = [  -slat*clon  -slat*slon   clat; ...
       -slon          clon         0; ... 
       -clat*clon  -clat*slon  -slat];

[MAX, N] = size(ned);
ecef = zeros(MAX, N);

for i=1:MAX
    
    ecef_t = ecef_org + R' * ned(i, :)';
    ecef (i, :) = ecef_t';
end

