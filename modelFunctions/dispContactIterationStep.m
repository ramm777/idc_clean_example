function dispContactIterationStep(s, connodes)


% SYNOPSIS:
%   dispContactIterationStep(s, connodes)
%
% DESCRIPTION: Displays the contact iteration step during simulation with the relevant data
%                             
% 
% PARAMETERS:
%   s         - contact iteration step
%   connodes  - connection nodes in the contact interaction algorithm
% 
% 
% RETURNS:
%   displays the contact ireration step with contacting nodes data in the command window
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
            
        disp('------------------------------')
        disp(['Contact Iteration Step: ', num2str(s)])
        disp(['Connecting ', num2str(numel(connodes(:, 1))), ' pairs:'])
        if numel(connodes(:,1)) >= 2
            disp(num2str(connodes(1:2,:))) % display the first 3
            disp('........'); 
        else
            disp(num2str(connodes(1,:)))   % display the first one
        end
        
        
end