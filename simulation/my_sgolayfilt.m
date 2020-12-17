function samples_s = my_sgolayfilt (samples)
% my_sgolayfilt: Savitzky-Golay Filtering with variable frame length.
%
% INPUT
%		samples: Nx1 samples to be smoothed.
%
% OUTPUT
%		samples_s: Nx1 smoothed samples.
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
%   MATLAB help for sgolayfilt function.
%
% Version: 001
% Date:    2018/09/19
% Author:  Rodrigo Gonzalez <rodralez@frm.utn.edu.ar>
% URL:     https://github.com/rodralez/navego

samples_size = max(size(samples));

sgo_framelen = floor( samples_size / 5);

% sgo_framelen should not be greater than 21.
if (sgo_framelen > 21)
    sgo_framelen = 21;
end

% sgo_framelen must be odd
if ( mod(sgo_framelen, 2) == 0 )
    sgo_framelen = sgo_framelen +1;
end

sgo_order = ceil( sgo_framelen / 2);

samples_s = sgolayfilt(samples, sgo_order, sgo_framelen);

end
