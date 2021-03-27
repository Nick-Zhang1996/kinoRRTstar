

function [i] =  RRTstar3D(dim, segmentLength, radius, show_output, samples)

coder.extrinsic('tic');
coder.extrinsic('toc');

firstSol=0;

% planning in state space
if dim == 2   
    state_dim = 6;
    start_cord = [2,2,0,0];
    goal_cord = [18,18,0,0];
elseif dim == 3
    state_dim = 9;
    start_cord = [2,2,2,0,0,0,0,0,0,];
    goal_cord = [18,8,8,0,0,0,0,0,0,];
end

% create random world
if dim == 2
    Size = 20;
    [world,~] = createKnownWorld(ones(1,dim)*Size,zeros(1,dim),dim);
elseif dim == 3
    [world,~] = createKnownWorld([20 10 10],zeros(1,dim),dim);
end
fprintf("world created... \n");

% Each Row Contains States, ConnectToEnd flag, Costz, ParentNodeIdx, and ChildNum
start_node = [start_cord,0,0,0,0];                
end_node = [goal_cord,0,0,0,0];

% establish tree starting with the start node
tree = [];
coder.varsize('tree');
tree = start_node;

%GChild  = [];
%coder.varsize('GChild')
GChild = zeros(samples, samples);

% tic

% check to see if start_node connects directly to end_node
if norm(start_node(1:dim)-end_node(1:dim))<segmentLength && collision(start_node,end_node,world,dim,state_dim)==0
  path = [start_node; end_node];
else
  if samples >0
  numPaths = 0;
  for i = 1:samples
      %fprintf("trying node... \n");
      [tree,GChild,flag] = extendTree(tree,GChild,end_node,segmentLength,radius,world,0,dim,state_dim,numPaths);
      %fprintf("nodes added");
      numPaths = numPaths + flag;

      if flag == 1 && numPaths==1
          min_cost = findMinCost(tree,end_node,dim,state_dim);
          toc
          fprintf("nodes: %.0f, min cost %.3f \n",i, min_cost);
      end
      
      if ismember(i, [400,600,800,1000])
         min_cost = findMinCost(tree,end_node,dim,state_dim);
         toc
         fprintf("nodes: %.0f, min cost %.3f \n",i, min_cost);
      end 

  end

  else
  i = 0;
  numPaths = 0;
  while numPaths < 1
      [tree,GChild,flag] = extendTree(tree,GChild,end_node,segmentLength,radius,world,0,dim,state_dim,numPaths);
      numPaths = numPaths + flag;
      i = i+1;
  end
  end

end


% if show_output == 1
% figure;
% [path_first,cost] = findMinimumPath(tree_small,end_node,dim,state_dim);
% Cost=[Cost, cost];
% plotExpandedTree(tree_small,dim,state_dim);
% plotWorld(world,dim);
% plotTraj(path_first,dim);
% figure;
% [path_400,cost] = findMinimumPath(tree_400,end_node,dim,state_dim);
% Cost=[Cost, cost];
% plotExpandedTree(tree_400,dim,state_dim);
% plotWorld(world,dim);
% plotTraj(path_400,dim);
% % figure;
% [path_600,cost] = findMinimumPath(tree_600,end_node,dim,state_dim);
% Cost=[Cost, cost];
% % plotExpandedTree(tree_600,dim,state_dim);
% % plotWorld(world,dim);
% % plotTraj(path_600,dim);
% % figure;
% [path_800,cost] = findMinimumPath(tree_800,end_node,dim,state_dim);
% Cost=[Cost, cost];
% % plotExpandedTree(tree_800,dim,state_dim);
% % plotWorld(world,dim);
% % plotTraj(path_800,dim);
% % figure;
% [path_1000,cost] = findMinimumPath(tree_1000,end_node,dim,state_dim);
% Cost=[Cost, cost];
% % % plotExpandedTree(tree_1000,dim,state_dim);
% % plotWorld(world,dim);
% % plotTraj(path_1000,dim);
% % figure;
% [path_1200,cost] = findMinimumPath(tree_1200,end_node,dim,state_dim);
% Cost=[Cost, cost];
% % plotExpandedTree(tree_1200,dim,state_dim);
% % plotWorld(world,dim);
% % plotTraj(path_1200,dim);
% % figure;
% [path_1400,cost] = findMinimumPath(tree_1400,end_node,dim,state_dim);
% Cost=[Cost, cost];
% % plotExpandedTree(tree_1400,dim,state_dim);
% % plotWorld(world,dim);
% % plotTraj(path_1400,dim);
% % figure;
% [path_1600,cost] = findMinimumPath(tree_1600,end_node,dim,state_dim);
% Cost=[Cost, cost];
% % plotExpandedTree(tree_1600,dim,state_dim);
% % plotWorld(world,dim);
% % plotTraj(path_1600,dim);
% % figure;
% [path_1800,cost] = findMinimumPath(tree_1800,end_node,dim,state_dim);
% Cost=[Cost, cost];
% % plotExpandedTree(tree_1800,dim,state_dim);
% % plotWorld(world,dim);
% % plotTraj(path_1800,dim);

% figure;
% [path,~] = findMinimumPath(tree,end_node,dim,state_dim);
% % Cost=[Cost, cost];
% % plotExpandedTree(tree,dim,state_dim);
% plotWorld(world,dim); hold on
% plotTraj(path,dim);
% end

end


function [new_tree, GChild, flag] = extendTree(tree, GChild, end_node, segmentLength, r, world, flag_chk, dim, state_dim, numPaths)   % segmentLength: maximum stepsize, r: neighbor radius
  new_node  = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
  new_tree = zeros(1,13);
  
  flag1 = 0;
  if numPaths == 0
      eta = 0;
  else
      eta = 0;
  end
  
  while flag1 == 0
    % select a random point
    if rand > eta
       randomPoint = zeros(1, state_dim);
       for i = 1 : dim
          randomPoint(1, i) = (world.endcorner(i) - world.origincorner(i)) * rand;
       end
    else
       randomPoint = end_node(1:state_dim);
    end

    % find node that is closest to randomPoint (Eucl. dist. between positions). 
    %fprintf("find node that is closest to randomPoint \n");
    % NOTE
%     tmp = tree(:, 1 : dim) - randomPoint(1 : dim);
    tmp_dist = zeros(size(tree,1), dim);
    for i = 1:size(tree,1)
      tmp_dist(i,1:dim) = tree(i, 1 : dim) - randomPoint(1 : dim);
    end
    
    sqr_dist = sqr_eucl_dist(tmp_dist, dim);
    [min_dist, idx] = min(sqr_dist);
    min_parent_idx = idx;
    
    Vect = randomPoint(1:dim)-tree(idx,1:dim);
    Vect=Vect/norm(Vect);
    
    % NOTE
    new_point = [0 0 0 0 0 0 0 0 0];
    % find new_point that is within the range of min_parent_idx in terms of segmentLength (Eucl. dist. between positions).
    %fprintf("find new_point that is within the range of min_parent_idx in terms of segmentLength (Eucl. dist. between positions)\n");
    if min_dist > segmentLength^2
        % generate a new point that is closest to randomPoint, segmentLength away from tree(idx,1:dim)
        new_point(1 : dim) = tree(idx, 1 : dim) + Vect * segmentLength;
    else
        new_point(1 : dim) = randomPoint(1 : dim);
    end

  % check if the new_point is in collision
  %fprintf("check if the new_point is in collision\n");
  if collision_point(new_point, world, dim) == 0
      [collision_flag, Tf, xf] = collisionFreeVel(tree(idx, :), new_point, world, dim, state_dim);    % this collision checking includes a steering function and forward simulation
                                                                                                                                         
    if collision_flag == 0  
      %fprintf("no collision \n");
      new_point(dim + 1 : state_dim) = real(xf(dim + 1 : state_dim));
      [min_cost, ~] = cost_npFreeVel(tree(idx,:), new_point, dim, state_dim, Tf);                 % total cost from root to new_point through its parent tree(idx,:)

%       tmp_dist = tree(:, 1 : dim) - new_point(1 : dim);
      tmp_dist = zeros(size(tree,1), dim);
      for i = 1:size(tree,1)
        tmp_dist(i,1:dim) = tree(i, 1 : dim) - new_point(1 : dim);
      end
      
      dist_sqr = sqr_eucl_dist(tmp_dist, dim);
         
      % find near neighbors   
      if dim == 2
        gamma = 60;
      elseif dim == 3
        gamma = 40; 
      end      
      nun = size(tree, 1);
      ner = gamma * ( log(nun + 1) / nun )^(1 / dim);
      r1 = min(ner, r);    
      near_idx = find(dist_sqr <= r1^2);
      
      size_near = size(near_idx, 1);
      %fprintf("choosing new parent\n")
      if size_near >= 1          
        for i = 1 : size_near                                                % choose parent node
            if near_idx(i) ~= idx
                [segcost, Tf] = segment_costFreeVel(tree(near_idx(i), :), new_point, dim);
                cost_near = tree(near_idx(i), state_dim + 2) + segcost;
                if cost_near + 0.001 < min_cost
                    [collision_flag, ~, xf] = collisionFreeVelKnowTf(tree(near_idx(i), :), new_point, world, dim, state_dim, Tf);
                    if collision_flag==0      
                        min_cost = cost_near;
                        min_parent_idx = near_idx(i);    
                        new_point(dim + 1 : state_dim) =   real(xf(dim + 1 : state_dim));
                    end
                end
            end
        end             
      end
      
      %fprintf("assembling new node \n");
      new_node = [new_point, 0 , min_cost, min_parent_idx, 0];
      new_tree = [tree; new_node];
      new_node_idx = size(new_tree, 1);      
      new_tree(min_parent_idx, state_dim + 4) = new_tree(min_parent_idx, state_dim + 4) + 1;     % ChildNum + 1
      %fprintf("updating GChild \n");
      GChild( min_parent_idx, new_tree(min_parent_idx, state_dim + 4) ) = new_node_idx;        % update GChild matrix
      
      %fprintf("rewiring \n");
      if size_near >= 1                                                  % rewire
        reduced_idx = near_idx;
        for j = 1 : size_near
            if reduced_idx(j) ~= min_parent_idx
                near_cost = new_tree(reduced_idx(j), state_dim + 2);
                [lcost, Tf] = segment_cost(new_point, new_tree(reduced_idx(j), :), dim);
                rnewcost = min_cost + lcost;
                if near_cost > rnewcost + 0.001
                    [collision_flag, ~] = collisionKnowTf(new_node, new_tree(reduced_idx(j), :), world, dim, state_dim, Tf);
                    if collision_flag == 0 
           
                    ecost = rnewcost - near_cost;                
                    GChild( new_tree(reduced_idx(j),state_dim+3), GChild( new_tree(reduced_idx(j),state_dim+3),: ) == reduced_idx(j) ) = -1;      % parent of reduced_idx(j) before rewire, change its child list.
                    new_tree(reduced_idx(j),state_dim+2) = rnewcost;           % update the cost and parent information of the node being rewired, reduced_idx(j)
                    new_tree(reduced_idx(j),state_dim+3) = new_node_idx;
                                
                    new_tree(new_node_idx, state_dim+4) = new_tree(new_node_idx, state_dim+4) + 1;         % add the node being rewired to the child list of the new added node, new_node_idx.
                    GChild( new_node_idx, new_tree(new_node_idx, state_dim+4) ) = reduced_idx(j) ;
                
                    brunchCost(new_tree, GChild, reduced_idx(j), ecost, dim, state_dim);       % update all cost of the descendant of the node being rewired
                
                    end
                end
            end
        end
      end
      %fprintf("rewiring done \n");
      flag1=1;
    end
    
  end
    
  end

  flag = 0;
  if flag_chk == 0
    % check to see if new node connects directly to end_node
    if  norm(new_node(1 : dim) - end_node(1 : dim)) < segmentLength 
        [collision_flag, ~] = collision(new_node, end_node, world, dim, state_dim);        
        if  collision_flag == 0
            flag = 1;
            new_tree(end, state_dim + 1) = 1;  % mark node as connecting to end.
        end
    end
  end
end

% update all cost of the descendant of the node being rewired, id_candidate
function brunchCost(new_tree, GChild, id_candidate, ecost, dim, state_dim)

   for i = 1 : new_tree(id_candidate,state_dim+4) 
        if GChild(id_candidate,i) ~= -1
            new_tree(GChild(id_candidate,i), state_dim+2) = new_tree(GChild(id_candidate,i),state_dim+2) + ecost;
            brunchCost(new_tree, GChild, GChild(id_candidate,i), ecost, dim, state_dim);
        end        
   end

end




