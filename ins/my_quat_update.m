function [qua_n, DCMbn_n, euler] = my_quat_update(wb, qua, dt)
% att_update: updates attitude using quaternion or DCM.
%
% INPUT:
%   wb,         3x1 incremental turn-rates in body-frame (rad/s).
%   DCMbn,      3x3 body-to-nav DCM.
%   qua,        4x1 quaternion.
%   omega_ie_N, 3x1 Earth rate (rad/s).
%   omega_en_N, 3x1 Transport rate (rad/s).
%   dt,         1x1 IMU sampling interval (s).
%	att_mode,   attitude mode string.
%      'quaternion': attitude updated in quaternion format. Default value.
%             'dcm': attitude updated in Direct Cosine Matrix format.
%
% OUTPUT:
%   qua_n,      4x1 updated quaternion.
%   DCMbn_n,    3x3 updated body-to-nav DCM.
%   euler,      3x1 updated Euler angles (rad).
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
%	Crassidis, J.L. and Junkins, J.L. (2011). Optimal Esti-
% mation of Dynamic Systems, 2nd Ed. Chapman and Hall/CRC, USA.
% Eq. 7.39, p. 458.
%
% Version: 003
% Date:    2016/11/26
% Author:  Rodrigo Gonzalez <rodralez@frm.utn.edu.ar>
% URL:     https://github.com/rodralez/navego


%% Ignore other transport rate and Earth rate

wb_n = wb;

%% Quaternion update   

    qua_n   = qua_update(qua, wb_n, dt);    % Update quaternion
    qua_n   = qua_n / norm(qua_n);          % Brute-force normalization
    DCMbn_n = qua2dcm(qua_n);               % Update DCM
    DCMbn_n = reshape(DCMbn_n, 1, 9);       % Convert into a vector
    euler   = qua2euler(qua_n);    % Update Euler angles

end

