from kinoRRT import *

world = World(20,10,10)
obs1 = Box(6,0,0,6+2,5,10)
obs2 = Box(6,5,0,6+2,5+5,5)
obs3 = Box(12,5,0,12+2,5+5,10)
obs4 = Box(12,0,5,12+2,5,5+5)

world.addObstacle(obs1)
world.addObstacle(obs2)
world.addObstacle(obs3)
world.addObstacle(obs4)

start_node = Node(2,2,2)
goal_node = Node(18, 8, 8)

rrt = KinoRrtStar(world, start_node, goal_node, 600, 10)

try:
    rrt.run()
except (RuntimeError):
    print("error")
print("rrt run finished")

'''
waypoint_n = rrt.prepareSolution()

waypoints = []
for i in range(waypoint_n):
    p = rrt.getNextWaypoint()
    assert (p.valid)
    waypoints.append( (p.x,p.y,p.z) )

print(waypoints)
'''
