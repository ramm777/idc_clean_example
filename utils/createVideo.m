function createVideo(F, videoname)


% SYNOPSIS:
%   createVideo(F, videoname)
%
% DESCRIPTION: Create video based on the pre-recorded frames (F) at each simulation step
%              
%                             
% 
% PARAMETERS:
%   F         - pre-recorded frames of a figure at each simulation step, must be a cell structure  
%   videoname - name of the file of the video, usually mfilename
% 
% 
% RETURNS:
%             - saves the video in the current folder
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


    F2 = cell2mat(F);
   
    videoname = strcat(videoname, '.avi'); 
    video = VideoWriter(videoname); 
    video.FrameRate = 2;                    
    open(video); 
    writeVideo(video, F2); 
    close(video);
    
    disp(['Recorded video: ', videoname]); 

    
end