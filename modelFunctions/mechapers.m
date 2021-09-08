function [allmechapers, lensA] = mechapers(G, fnodes, s, allmechapers, lensA, varargin)
%
%
% SYNOPSIS:
%   [allmechapers, lensA] = mechapers(G, fnodes, s)
%
% DESCRIPTION: Calculate mech apertures at each simulation step. (Can be improved later by finding which nodes of intracells 
%              belong to RHSbody and LHSbody. But this will slow down the simulation).  
%                             
% 
% PARAMETERS:
%   G         - Grid structure
%   fnodes    - fracture nodes, must be sorted [LHSbody, RHSbody]
%   s         - simulation step
% 
% 
% RETURNS:
%   allmechapers  - all mechanical apertures at each simulation step
%   lensA         - lengths between the apertures 
%   
%
% EXAMPLE:
%
% SEE ALSO:
%
% AUTHOR:
%   Amanzhol Kubeyev

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
    
    
    opt     = struct( 'mode', 'standard', ...
                      'outputSameSize', 1); 
    opt     = merge_options(opt, varargin{:}); 
    
    
    % -------------------------------------------------------------------------
    
    Lcoords = G.nodes.coords(fnodes(:, 1), :); 
    Rcoords = G.nodes.coords(fnodes(:, 2), :);
    mecha = Rcoords(:, 1) -  Lcoords(:, 1);    
    
    
    switch opt.mode
        case('standard')
            
            if any(mecha < 0)
                warning('Negative apertures, forcing them to 0')
                mecha(mecha < 0) = 0;
            end
                   
        case('distance')
            
            distances = sqrt(sum(bsxfun(@minus, Rcoords, Lcoords).^2,2));
            
            if any(mecha < 0)
                warning('Negative apertures, forcing them to 0');
                
                neg = find(mecha < 0); 
                pos = find(mecha >= 0);
                mecha_n = nan(numel(fnodes(:, 1)), 1);
                mecha_n(neg, 1) = 0; 
                mecha_n(pos, 1) = distances(pos, 1);
                mecha = mecha_n; 
            else    
                mecha = distances; 
                
            end
 
    end
    
    if opt.outputSameSize == 1
        allmechapers(:, s) = mecha(2 : end);  % Delete the first raw to make size of matrix equal to size of lensA mat
    else
        allmechapers(:, s) = mecha;
    end
        
    % -------------------------------------------------------------------------
    
    lensL = abs(diff(Lcoords(:, 2)));
    lensR = abs(diff(Rcoords(:, 2)));        
    lensA(1 : numel(fnodes(:, 1)) - 1, s) = 0.5*(lensL + lensR);  
    
    
end    
    