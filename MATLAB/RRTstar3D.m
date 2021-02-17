function [its, sizePath] =  RRTstar3D(dim, segmentLength, radius, random_world, show_output, samples) %#codegen

% firstSol=0;

% Planning in state space
if dim == 2   
    start_cord = [2,2,0,0];
    goal_cord = [18,18,0,0];
else
    start_cord = [2,18,2,0,0,0];
    goal_cord = [18,2,18,0,0,0];
end

% Create random world
Size = 20;
if random_world == 0
    [world,~] = createKnownWorld(ones(1,dim)*Size,zeros(1,dim),dim);
end

% Each Row Contains States, ConnectToEnd flag, Costz, ParentNodeIdx, and ChildNum
start_node = [start_cord,0,0,0,0];                
end_node = [goal_cord,0,0,0,0];

% Establish tree starting with the start node
tree = [];
coder.varsize('tree');
tree = start_node;


%coder.varsize('GChild')
GChild  = zeros(samples,samples);
%GChild = [0];
%tic

% check to see if start_node connects directly to end_node
if norm(start_node(1:dim)-end_node(1:dim))<segmentLength && collision(start_node,end_node,world,dim)==0
  path = [start_node; end_node];
else
  if samples >0
  its = 0;
  numPaths = 0;

  found_sol = 0;
  for i = 1:samples

      [tree,GChild,flag] = extendTree(tree,GChild,end_node,segmentLength,radius,world,0,dim);
  
      numPaths = numPaths + flag;
      its = its+1;
%       if mod(i,100) == 0
%           fprintf("%.0f nodes \n",i);
%       end
      
      % report first solution
      % its, time, cost
      if flag == 1 && numPaths==0
          min_cost = findMinCost(tree,end_node,dim);
          %toc
          fprintf("nodes: %.0f, min cost %.3f \n",i, min_cost);
      end
      
      if ismember(i, [400,1000,2000,3000,4000])
         min_cost = findMinCost(tree,end_node,dim);
         %toc
         fprintf("nodes: %.0f, min cost %.3f \n",i, min_cost);
      end 
      
%       if numPaths==1 && firstSol==0
%         toc  
%         firstSol=1;
%         tree_small = tree;
%         its
%       end
%       if its==500
%         toc
%         tree_500 = tree;
%       end
%       if its==1000
%         toc
%         tree_1000 = tree;
%       end
%       if its==1500
%         toc
%         tree_1500 = tree;
%       end
  end

  
  else
  its = 0;
  numPaths = 0;
  while numPaths < 1
      [tree,GChild,flag] = extendTree(tree,GChild,end_node,segmentLength,radius,world,0,dim);
      numPaths = numPaths + flag;
      its = its+1;
  end
%   its
  end

end

% toc
% run_time = toc

if show_output == 1
% figure;
% path_first = findMinimumPath(tree_small,end_node,dim);
% plotExpandedTree(world,tree_small,dim);
% plotWorld(world,path_first,dim);
% plotTraj(path_first,dim);
% figure;
% path_500 = findMinimumPath(tree_500,end_node,dim);
% plotExpandedTree(world,tree_500,dim);
% plotWorld(world,path_500,dim);
% plotTraj(path_500,dim);
% figure;
% path_1000 = findMinimumPath(tree_1000,end_node,dim);
% plotExpandedTree(world,tree_1000,dim);
% plotWorld(world,path_1000,dim);
% plotTraj(path_1000,dim);
% figure;
% path_1500 = findMinimumPath(tree_1500,end_node,dim);
% plotExpandedTree(world,tree_1500,dim);
% plotWorld(world,path_1500,dim);
% plotTraj(path_1500,dim);
% figure;

% figure;
path = findMinimumPath(tree,end_node,dim);
sizePath = size(path,1);
% plotExpandedTree(world,tree,dim);
% plotWorld(world,path,dim); hold on
% plotTraj(path,dim);
end

end


function [new_tree, GChild, flag] = extendTree(tree, GChild, end_node, segmentLength, r, world, flag_chk, dim)   % segmentLength: maximum stepsize, r: neighbor radius
  % NOTE
  new_node  = [0, 0, 0, 0, 0];
  new_tree = zeros(1,8);
  
  flag1 = 0;
  while flag1 == 0
    % select a random point
    if rand > 0
       randomPoint = zeros(1, 2 * dim);
       for i = 1 : dim
          randomPoint(1, i) = (world.endcorner(i) - world.origincorner(i)) * rand;
       end
    else
       randomPoint = end_node(1:2*dim);
    end

    % find node that is closest to randomPoint (Eucl. dist. between positions). 
    %tmp = tree(:, 1 : dim) - randomPoint(1 : dim);
    tmp = tree(:, 1 : dim);
    for ind = 1:size(tmp,1)
        tmp(ind,:) = tmp(ind,:) - randomPoint(1 : dim);
    end
    sqr_dist = sqr_eucl_dist(tmp, dim);
    [min_dist, idx] = min(sqr_dist);
    min_parent_idx = idx;
    
    Vect = randomPoint(1:dim)-tree(idx,1:dim);
    Vect=Vect/norm(Vect);
    new_point = [0 0 0 0];
    if dim == 2
        new_point(dim + 1 : 2 * dim) = 1.0 * [-1 + (1 - (-1)) * rand, -1 + (1 - (-1)) * rand];     % vmin = -1�� vmax = 1
%         if (size(new_point(dim + 1 : 2 * dim),1) ~= 1)
%             display("xx")
%         end
    elseif dim == 3
        new_point(dim + 1 : 2 * dim) = [-1 + (1 - (-1)) * rand, -1 + (1 - (-1)) * rand, -1 + (1 - (-1)) * rand];     % vmin = -1, vmax = 1
    end
    
    % find new_point that is within the range of min_parent_idx in terms of segmentLength (Eucl. dist. between positions).
    if min_dist > segmentLength^2
        % generate a new point that is closest to randomPoint, segmentLength away from tree(idx,1:dim)
        new_point(1 : dim) = tree(idx, 1 : dim) + Vect * segmentLength;
    else
        new_point(1 : dim) = randomPoint(1 : dim);
    end

  % check if the new_point is in collision
  if collision_point(new_point, world, dim) == 0
      [collision_flag, Tf] = collision(tree(idx, :), new_point, world, dim);    % this collision checking includes a steering function
                                                                                                                                         
    if collision_flag == 0  
        
      [min_cost, ~] = cost_np(tree(idx,:), new_point, dim, Tf);                 % total cost from root to new_point through its parent tree(idx,:)
      new_node  = [new_point, 0, min_cost, idx, 0];                             % new node candidate

      %tmp_dist = tree(:, 1 : dim) - new_point(1 : dim);
      tmp_dist = tree(:, 1 : dim);
      for ind = 1:size(tree,1)
          tmp_dist(ind,:) = tmp_dist(ind,:) - new_point(1 : dim);
      end
      dist_sqr = sqr_eucl_dist(tmp_dist, dim);
         
      % find near neighbors   
      if dim == 2
        gamma = 40;
      elseif dim == 3
        gamma = 24; 
      end      
      nun = size(tree, 1);
      ner = gamma * ( log(nun + 1) / nun )^(1 / dim);
      r1 = min(ner, r);    
      near_idx = find(dist_sqr <= r1^2);
      
      size_near = size(near_idx, 1);
      if size_near >= 1          
        for i = 1 : size_near                                                % choose parent node
            if near_idx(i) ~= idx
                [segcost, Tf] = segment_cost(tree(near_idx(i), :), new_point, dim);
                cost_near = tree(near_idx(i), 2 * dim + 2) + segcost;
                if cost_near + 0.01 < min_cost
                    [collision_flag, ~] = collisionKnowTf(tree(near_idx(i), :), new_node, world, dim, Tf);
                    if collision_flag==0      
                        min_cost = cost_near;
                        min_parent_idx = near_idx(i);              
                    end
                end
            end
        end             
      end
      
      new_node = [new_point, 0 , real(min_cost), min_parent_idx, 0];
      new_tree = [tree; new_node];
      % DEBUG check new node cost
      %fprintf("new node cost %.2f",real(min_cost));
      new_node_idx = size(new_tree, 1);      
      new_tree(min_parent_idx, 2 * dim + 4) = new_tree(min_parent_idx, 2 * dim + 4) + 1;     % ChildNum + 1

      GChild( min_parent_idx, new_tree(min_parent_idx, 2 * dim + 4) ) = new_node_idx;        % update GChild matrix

      if size_near >= 1                                                  % rewire
        reduced_idx = near_idx;
        for j = 1 : size_near
            if reduced_idx(j) ~= min_parent_idx
                near_cost = new_tree(reduced_idx(j), 2 * dim + 2);
                [lcost, Tf] = segment_cost(new_point, new_tree(reduced_idx(j), :), dim);
                rnewcost = min_cost + lcost;
                if near_cost > rnewcost + 0.01
                    [collision_flag, ~] = collisionKnowTf(new_node, new_tree(reduced_idx(j), :), world, dim, Tf);
                    if collision_flag == 0 
           
                    ecost = rnewcost - near_cost;                
                    GChild( new_tree(reduced_idx(j),2*dim+3), GChild( new_tree(reduced_idx(j),2*dim+3),: ) == reduced_idx(j) ) = -1;      % parent of reduced_idx(j) before rewire, change its child list.
                    new_tree(reduced_idx(j),2*dim+2) = real(rnewcost);           % update the cost and parent information of the node being rewired, reduced_idx(j)
                    new_tree(reduced_idx(j),2*dim+3) = new_node_idx;
                                
                    new_tree(new_node_idx, 2*dim+4) = new_tree(new_node_idx, 2*dim+4) + 1;         % add the node being rewired to the child list of the new added node, new_node_idx.
                    GChild( new_node_idx, new_tree(new_node_idx, 2*dim+4) ) = reduced_idx(j) ;
                
                    brunchCost(new_tree, GChild, reduced_idx(j), ecost, dim);       % update all cost of the descendant of the node being rewired
                
                    end
                end
            end
        end
      end
      flag1=1;
    end
    
  end
    
  end

  flag = 0;
  if flag_chk == 0
    % check to see if new node connects directly to end_node
    if  norm(new_node(1 : dim) - end_node(1 : dim)) < 0.2
        flag = 1;
        new_tree(end, 2 * dim + 1) = 1;  % mark node as connecting to end.
    elseif  norm(new_node(1 : dim) - end_node(1 : dim)) < segmentLength 
        [collision_flag, ~] = collision(new_node, end_node, world, dim);        
        if  collision_flag == 0
            flag = 1;
            new_tree(end, 2 * dim + 1) = 1;  % mark node as connecting to end.
        end
    end
  end
end

% update all cost of the descendant of the node being rewired, id_candidate
function brunchCost(new_tree, GChild, id_candidate, ecost, dim)

   for i = 1 : new_tree(id_candidate,2*dim+4) 
        if GChild(id_candidate,i) ~= -1
            new_tree(GChild(id_candidate,i), 2*dim+2) = new_tree(GChild(id_candidate,i),2*dim+2) + real(ecost);
            brunchCost(new_tree, GChild, GChild(id_candidate,i), ecost, dim);
        end        
   end

end




