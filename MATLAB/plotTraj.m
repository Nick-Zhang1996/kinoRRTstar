% plot trajectory from path and NN controller
function plotTraj(path,dim)

global obj
for i=1:size(path,1)-1
    
    if dim==2     
        
    in = num2cell([path(i, 1 : 4), path(i + 1, 1 : 4)]);
%     Tf = norm(path(i, 1 : 4) - path(i + 1, 1 : 4))/0.8;
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
    
    plot(path(i,1),path(i,2),'Marker','o','MarkerSize',8,'MarkerEdgeColor','[0.1 0.1 0.1]','LineWidth',1.5);
    plot(traj(1,1:size(t,2)),traj(2,1:size(t,2)),'color',[0.8500 0.3250 0.0980],'LineWidth',2);hold on
    xl = xlabel('$x (m)$','Interpreter','LaTeX');
    yl = ylabel('$y (m)$','Interpreter','LaTeX');
    set(xl,'FontSize',18);
    set(yl,'FontSize',18);
    set(gca,'FontSize',16,'FontName','Times');
%     figure(2)
%     plot(t,traj(1,1:size(t,2)),'LineWidth',2);hold on
%     plot(t,traj(2,1:size(t,2)),'LineWidth',2);
%     figure(3)
%     plot(t,traj(3,1:size(t,2)),'LineWidth',2);hold on
%     plot(t,traj(4,1:size(t,2)),'LineWidth',2);hold on
    
%     plot(traj(1,end),traj(2,end),'o','LineWidth',2,'MarkerSize',8);hold on
%     plot(path(i+1,1),path(i+1,2),'o','LineWidth',2,'MarkerSize',8);hold on

    elseif dim == 3
        return;
    end        
       
end
end