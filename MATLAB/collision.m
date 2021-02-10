function [collision_flag, Tf] = collision(parent, node, world, dim)

collision_flag = 0;

if collision_flag == 0 && dim == 2

        tf_temp = roots(opt_time(parent(1), parent(2), parent(3), parent(4), node(1), node(2), node(3), node(4)));
        tf_temp = tf_temp(imag(tf_temp) == 0);
        mask = tf_temp >=0;
        temp = tf_temp(mask);
        %NOTE
        Tf = min(temp);
        %Tf = sum(temp);

    
    states = @(t)state_traj(t, Tf, parent(1), parent(2), parent(3), parent(4), node(1), node(2), node(3), node(4));
    checkpoints = 10;
    state = zeros(2*dim, checkpoints);
    t = linspace(0, Tf, checkpoints + 2);
    for i = 1 : checkpoints
        % added real()
        state(:, i) = real(states(t(i + 1)));
    end
    traj = state(1 : 2, :);
    
    for j = 1 : checkpoints
    p = traj(:, j);
    
    for i=1:dim
        if (p(i)>world.endcorner(i))||(p(i)<world.origincorner(i))
            collision_flag = 1;
            return;    % break;
        end
    end
    
      % check each obstacle
      for i=1:world.NumObstacles
        if sum(([p(1);p(2)]-[world.cx(i); world.cy(i)]).*([p(1);p(2)]-[world.cx(i); world.cy(i)]))<(world.radius(i)+0.1)^2  % (norm([p(1);p(2)]-[world.cx(i); world.cy(i)])<=1*world.radius(i))
            collision_flag = 1;
            return;    % break;
        end
      end
    end

%%%%% dim=3 case has not been implemented yet %%%%%
elseif collision_flag == 0 && dim ==3
    collision_flag = 1;
    Tf = 0;
    
end

end