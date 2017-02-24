function color_vector=generate_color_map(numColors, control_pts)
% GENERATE_COLOR_MAP  Create color map
%
% Use as
%   generate_color_map(255, [340 225 90 10 10; 340 340 340 270 10; 340 340 340 340 340])
% where the first argument is the number of colors, and argument is a matrix containing
% the control points for the interpolation, one row per RGB, to generate the map with

    color_vector = [];

    
    for color=1:numColors
        cx = 10 + ((color-1)/(numColors-1))*(440);
        for j=1:size(control_pts,2)-1
            thisPt = control_pts(:,j);
            nextPt = control_pts(:,j+1);
            thisPtX=10 + (460-20) / (size(control_pts,2) - 1) * (j-1);
            nextPtX=10 + (460-20) / (size(control_pts,2) - 1) * (j);
                
            % Find the points to interpolate between
            if(cx >= thisPtX && cx <= nextPtX)
                cyInterp = thisPt + (cx - thisPtX)/(nextPtX - thisPtX)*(nextPt - thisPt);

                if(cyInterp > 350 - 10)
                    cyInterp = 350 - 10;
                elseif (cyInterp < 10)
                    cyInterp = 10;
                end
                color_vector(end+1,:) = round(255*(1-(cyInterp-10)/(350-20)));
                break

            else
                continue            
            end
        end
    end

