function h = plotMergeFigures(fig1, fig2, varargin)
%
% SYNOPSIS:
%   h = mergeFigures(fig1, fig2, varargin)
%
% DESCRIPTION: Function merges 2 or more figures into one. Figures which are to be 
%              merged must be in path of MATLAB for it to find them. 
% PARAMETERS:
%   fig1       - 1st figure's name 'your figure'
%   fig2       - 2nd figure's name to be merged to
%   ...        - more figures  
%   
% RETURNS:
%              singlee figure where figures were merged
%   
%
% EXAMPLE:
%
% SEE ALSO:
%
% AUTHOR:
%   Amanzhol Kubeyev  

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



    %----------------------------------------------------------------------
    % CONTROL
    
    % All figures must be closed before this function is executed
    close all
   
    % Screensize for future figures 
    ss = get(groot, 'Screensize'); 

    % Set default figure color
    set(groot,'defaultFigureColor','w')
    set(groot,'defaultAxesFontSize', 18)
    
    %----------------------------------------------------------------------
    % Open all specified figures
    
    openfig(fig1);    % 1
    openfig(fig2);    % 2 
    
    figures = {fig1; fig2}; 
  
    % if more than 2 figures are to be superimposed
    if nargin > 2
       
        % add figures fig3 fig4 etc.. 
        for i =  1 : length(varargin)
            
            openfig(varargin{i});
                        
            figures{2 + i, 1} = varargin{i}; 
            
        end
              
    end
    
    %----------------------------------------------------------------------
    % Merge all figures 
    
    for i = 2 : nargin
        copyobj(findobj(i - 1,'type','line'), findobj(i, 'type','axes'));
    end  
    
    %----------------------------------------------------------------------
    % Add legend: Warning always double check 
    
    figures = flip(figures); 
    legend(figures); 
            
    % close all figures except the last one
    close(1 : numel(figures(:, 1)) - 1); 
    
    axis auto
         

end