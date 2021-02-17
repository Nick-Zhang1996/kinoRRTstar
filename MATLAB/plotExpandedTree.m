function plotExpandedTree(world,tree,dim)

global obj
    ind = size(tree,1);
    while ind>1
    node = tree(ind,:);
    parent_node = tree(node(2*dim+3),:);
    
        ind = ind - 1;
        if dim == 2
            
    in = num2cell([parent_node(1 : 4), node(1 : 4)]);
    Tf = roots(obj.eval_arrival_internal(in{:}));
    Tf = Tf(abs(imag(Tf)) < 0.0001);
    Tf = min(Tf(Tf >= 0));

    states = @(t)obj.eval_states_internal(t, Tf, in{:});
    dt = 0.2;
    Horizon = Tf / dt;
    t = linspace(0,Tf,Horizon);
    state = zeros(2*dim, length(t));
    for j = 1 : Horizon
        state(:, j) = states(t(j));
    end

    traj = state(1 : 2, :);
    

    p = plot(traj(1,1:size(t,2)),traj(2,1:size(t,2)));
        set(p,'Color','[0.5 0.5 0.5]','LineWidth',0.5);
        hold on;   
        plot(node(1),node(2),'Marker','.','MarkerEdgeColor','[0.5 0.5 0.5]');

        elseif dim == 3
            return;
        end
    end
end