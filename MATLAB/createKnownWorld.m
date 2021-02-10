function [world, NumObstacles] = createKnownWorld(endcorner, origincorner, dim)
NumObstacles = 5;
  if dim == 2
  % check to make sure that the region is nonempty
    if (endcorner(1) <= origincorner(1)) || (endcorner(2) <= origincorner(2))
      disp('Not valid corner specifications!')
      world=[];
  % create world data structure
    else
    world.NumObstacles = NumObstacles;
    world.endcorner = endcorner;
    world.origincorner = origincorner;
    world.radius = [0 0 0 0 0];
    world.cx = [0 0 0 0 0];
    world.cy = [0 0 0 0 0];
    % create NumObstacles
    maxRadius = 3;

        world.radius(1) = maxRadius;
        cx = 6;
        cy = 6;
        world.cx(1) = cx;
        world.cy(1) = cy;

        world.radius(2) = maxRadius;
        cx = 15;
        cy = 14;
        world.cx(2) = cx;
        world.cy(2) = cy;

        world.radius(3) = maxRadius-1;
        cx = 10;
        cy = 11;
        world.cx(3) = cx;
        world.cy(3) = cy;

        world.radius(4) = maxRadius+0.5;
        cx = 5;
        cy = 15;
        world.cx(4) = cx;
        world.cy(4) = cy;

        world.radius(5) = maxRadius;
        cx = 15;
        cy = 5;
        world.cx(5) = cx;
        world.cy(5) = cy;
    end

  elseif dim == 3
        return;
   end
end