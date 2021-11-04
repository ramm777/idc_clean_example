function checkCurveContinuity(curve)
%
% SYNOPSIS:
%   checkCurveContinuity(curve)
%
% DESCRIPTION: 
%   Assure that curve is contiuity or not self-intersecting.
%
%
% PARAMETERS:
%   curve       - curve coordinates as a matrix [x, y] where each item (x or y) is a vector of coords.
%   
% RETURNS:
%   error flag  - if the curve is self intersecting.  
%   
%
% EXAMPLE:
%
% SEE ALSO:
%
% AUTHOR:
%   Aidan Kubeyev  

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



    % Collect all segments in matrix segmentnumber-vs-[x1 y1, x2 y2]
    allsegments = nan(numel(curve(:, 1)) - 1, 4); 
    for j = 2 : numel(curve(:, 1)) 
        allsegments(j-1, :) = [curve(j, :), curve(j-1, :)]; 
    end



    for i = 1 : numel(allsegments(:, 1)) 

        % Compare each segment [x1 y1, x2 y2] towards all except current
        int = lineIntersect(allsegments(i, :), allsegments(1:end ~= i, :)); 
        assert(~any(int.isintersec), ...
            ['Curve is self-intersecting at segment=', num2str(i)]); 

    end




end