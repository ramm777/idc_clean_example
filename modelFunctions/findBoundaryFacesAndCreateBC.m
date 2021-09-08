function bc = findBoundaryFacesAndCreateBC(G, varargin)


% SYNOPSIS:
%   bc = findBoundaryFacesAndCreateBC(G, varargin)
%
% DESCRIPTION: Find boundary faces and construct bc structure: bc type, value, sat 
%              (also see larger function in 3D)               
% 
% PARAMETERS:
%   G          - Grid structure
%   mode       - (optional) mode of collecting faces: 
%                               a) (default) 'pside', mrst standard way
%                               b) 'toleranceBased', find boundary faces based min/max dimensions of the grid
%   tolerance  - (optional) tolerance based on which the mode 'toleranceBased' will find the boundary faces
%   
% 
% RETURNS:
%   bc -  bc structure, where each item consists of: 
%                   a) faces
%                   b) value
%                   c) el_bc structure: disp_bc and force_bc (empty)
%
% EXAMPLE:
%
% SEE ALSO:
%
% AUTHOR:
%   Amanzhol Kubeyev, based on MRST existing routines

%{
Copyright 2009-2021 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}


    % -------------------------------------------------------------------------
    % Defaults and overwrite them
    
    opt     = struct( 'mode',          'pside',   ...
                      'tolerance',     0.000001 ); 
    opt     = merge_options(opt, varargin{:});

    
    % -------------------------------------------------------------------------
        
    switch opt.mode
        
        case('pside')
            
            oside = {'Left', 'Right', 'Back', 'Front'}; % 'Back' lower BC, 'Front' - upper BC
            bc = cell(4, 1);
            for i = 1 : numel(oside)
                bc{i} = pside([], G, oside{i}, 0);
                bc{i} = rmfield(bc{i}, 'type');
                bc{i} = rmfield(bc{i}, 'sat');
            end
               
        
        case('toleranceBased')
    
            for j = 1 : G.griddim
                Lmax = max(G.faces.centroids(:, j));
                Lmin = min(G.faces.centroids(:, j));
                x = [Lmin, Lmax];
                for i = 1 : 2

                    faces = find(abs(G.faces.centroids(:, j)-x(i)) < opt.tolerance);
                    assert(all(any(G.faces.neighbors(faces,:)==0,2)));

                    bc{i + (j - 1)*2} = addBC([], faces, 'pressure', 0);
                end
            end
    
      
    end

end