function [z2 , jrc] = JRCtse(curve, res, stepsize, varargin)
%
% SYNOPSIS:
%   jrc = JRCtse(curve, res, stepsize)
%
% DESCRIPTION: 
%   - For a given curve (heights) calculate Joint Roughness Coefficient based on Tse and Cruden, (1979)
%   - By default, roughness is calculated from adjacent points, but with stepsize, resolution can be changed. 
% 
% PARAMETERS:
%   curve    - Vector of nodes' heights describing roughness 
%   res      - Resolution or single distance btwn nodes in metres
%   stepsize - Downsizing coefficient of the data. If data is original and not downsized stepsize=1.
%              If downsized data, i.e. every 10th datapoints are taken, then stepsize=10   
%
%   'mode'   - (optional), mode of calculation: 1)'original' (default)  
%                                               2)'reproduce' - will reproduce the geometry based on the 
%                                                               existing roughness of a fracture
%
% RETURNS:
%   jrc - joint roughness coefficient
%   
%
% EXAMPLE:
%
% SEE ALSO:
%   JRCmaerz()
% 
% AUTHOR:
%   Aman Kubeyev  

%{
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


    opt = struct('mode', 'original');      % Default calculation mode
    opt = merge_options(opt, varargin{:}); % Override default control options  
    
    
    assert(isvector(curve), 'Curve must be a vector'); % Assert that curve is a vector
    if ~iscolumn(curve)                                % Make vector column vector if it's not
        curve = curve'; 
    end
    
    
    c = [curve, [0; cumsum(repmat(stepsize*res, length(curve) - 1, 1))]];  % Reconstruct curve's y-coords
    
    
    checkCurveContinuity(c); % Check data
    assert(length(curve) > 1 && ~any(isnan(curve)), ...
                                'JRC error:  data is too short or NaNs')
                            
   

    %--------------------------------------------------------------------------
    % Calculate z2 and JRC

    
    if strcmp(opt.mode, 'original')
        
        %--------------------------------------------------------------------------
        % Calculate JRC 
        
        heightsum = sum(diff(curve).^2); % Sum the squares of the differences in adjacent y-coordinates over the entire length L                        
        dx2 = (stepsize*res).^2;         % suqare of the equal intervals distance Dx
        m = length(curve) - 1;           % number of intervals (one fewer than number of points)
        z2 = sqrt( 1./(m * dx2) * heightsum ); 
        jrc = 32.2 + 32.37 * log10(z2); 

    
    elseif strcmp(opt.mode, 'reproduce')
        
        
        %--------------------------------------------------------------------------
        % Reproduce curve by finding lot's of intersections with the original curve 
        
        
        avcx = mean(c(:, 1));      % Get average x-coord
        n    = 50*length(curve);   % Number of horizontal lines

        ys = linspace(min(c(:, 2)), max(c(:, 2)), n)';   % min y-coord to max ycoord

        % Length which defines min xcoord and max xcoords in horzlines
        horzlinesLen =  2*norm(c(1,:) - c(end,:));

              
        % Construct horizontal lines by vectorization [x1, x2; y1, y2]   
        a = repmat([avcx - horzlinesLen, avcx + horzlinesLen], n, 1); 
        b = [ys, ys]; 
        
        
        % For each horzline, get the intersection [x, y]
        for i = 1 : n
            horzlines = [a(i, :); b(i, :)]; 
            t = curvesIntersect(horzlines, c')'; % [x1;y1],[x2;y2]
            intnodes(i, :) = t(1, :);            % If there are two nodes, take the first 
        
        end

        
        curve = intnodes(:, 1);                 % Reproduced curve, here x-coords
        nres = intnodes(2, 2)- intnodes(1, 2);   % Reproduced resolution
        
        %--------------------------------------------------------------------------
        % Calculate JRC on a reproduced curve
        
        heightsum = sum(diff(curve).^2); % Sum the squares of the differences in adjacent y-coordinates over the entire length L                        
        dx2 = (nres).^2;         % suqare of the equal intervals distance Dx
        m = length(curve) - 1;           % number of intervals (one fewer than number of points)
        z2 = sqrt( 1./(m * dx2) * heightsum ); 
        jrc = 32.2 + 32.37 * log10(z2); 
        
       
        
    end



end


%% APPENDIX


%--------------------------------------------------------------------------
% Plot-superimpose horizontal lines
% figure(1); hold on
% for i = 1 : n
%     
%     coords = horzlines{i};
%     coords = coords'; 
%     
%     plot(coords(:, 1), coords(:, 2), '-o'); 
% 
% end


%--------------------------------------------------------------------------
% Construct horizontal lines by looping [x1, x2; y1, y2]   
        %for i = 1 : n
        %   horzlines{i, 1} = [avcx - horzlinesLen, avcx + horzlinesLen; ys(i), ys(i)];    
        %end
        

