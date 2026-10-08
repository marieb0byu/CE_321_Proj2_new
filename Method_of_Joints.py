import sys
import Geometry_Operations as geom

# Bars at this node whose force has not been computed yet
def UnknownBars(node):
    return [bar for bar in node.bars if not bar.is_computed]

# A node is viable if it has one unknown bar, or two unknown bars that
# are not collinear (collinear bars can't both be solved from one node)
def NodeIsViable(node):
    unknown = UnknownBars(node)
    if len(unknown) == 1:
        return True
    if len(unknown) == 2:
        return abs(geom.SineBars(unknown[0], unknown[1])) > 1e-6
    return False

# Sum of all KNOWN forces at the node (loads, reactions, solved bars),
# resolved into the local frame whose x-axis points along local_x_bar
def KnownForcesInLocalFrame(node, local_x_bar):
    x_vec = geom.BarNodeToVector(node, local_x_bar)
    norm = geom.VectorTwoNorm(x_vec)
    load = [node.GetNetXForce(), node.GetNetYForce()]
    sum_x = geom.DotProduct(load, x_vec) / norm
    sum_y = geom.TwoDCrossProduct(x_vec, load) / norm
    for bar in node.bars:
        if bar.is_computed:
            bar_vec = geom.BarNodeToVector(node, bar)
            sum_x += bar.axial_load * geom.CosineVectors(x_vec, bar_vec)
            sum_y += bar.axial_load * geom.SineVectors(x_vec, bar_vec)
    return sum_x, sum_y

# Compute and store the force in local_x_bar from the sum of forces
# along its own axis, and mark it computed.
def SumOfForcesInLocalX(node, local_x_bar):
    sum_x, _ = KnownForcesInLocalFrame(node, local_x_bar)
    force = -sum_x
    local_x_bar.SetAxialLoad(force)
    local_x_bar.is_computed = True
    return force

# Compute and store the force in the second unknown bar from the sum of
# forces perpendicular to the first unknown bar (the local x bar), and
# mark it computed.
def SumOfForcesInLocalY(node, unknown_bars):
    x_bar, other_bar = unknown_bars[0], unknown_bars[1]
    _, sum_y = KnownForcesInLocalFrame(node, x_bar)
    force = -sum_y / geom.SineBars(x_bar, other_bar)
    other_bar.SetAxialLoad(force)
    other_bar.is_computed = True
    return force

def IterateUsingMethodOfJoints(nodes, bars):
    max_iterations = 10 * len(nodes) + 10
    iteration = 0
    while any(not bar.is_computed for bar in bars):
        iteration += 1
        if iteration > max_iterations:
            sys.exit("Method of joints did not converge; check the truss and supports")
        for node in nodes:
            if not NodeIsViable(node):
                continue
            unknown = UnknownBars(node)
            if len(unknown) == 2:
                SumOfForcesInLocalY(node, unknown)
            SumOfForcesInLocalX(node, unknown[0])
