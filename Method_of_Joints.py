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

# Force in local_x_bar from the sum of forces along its own axis.
# If another bar is still unknown, pass it and its (just solved) force.
def SumOfForcesInLocalX(node, local_x_bar, other_bar=None, other_force=0.0):
    sum_x, _ = KnownForcesInLocalFrame(node, local_x_bar)
    if other_bar is not None:
        sum_x += other_force * geom.CosineBars(local_x_bar, other_bar)
    return -sum_x

# Force in the second unknown bar from the sum of forces perpendicular to
# the first (the first bar has no component in local y)
def SumOfForcesInLocalY(node, unknown_bars):
    x_bar, other_bar = unknown_bars
    _, sum_y = KnownForcesInLocalFrame(node, x_bar)
    return -sum_y / geom.SineBars(x_bar, other_bar)

def IterateUsingMethodOfJoints(nodes, bars):
    while any(not bar.is_computed for bar in bars):
        progress = False
        for node in nodes:
            if not NodeIsViable(node):
                continue
            unknown = UnknownBars(node)
            x_bar = unknown[0]
            if len(unknown) == 2:
                other_bar = unknown[1]
                other_force = SumOfForcesInLocalY(node, unknown)
                other_bar.SetAxialLoad(other_force)
                other_bar.is_computed = True
                x_force = SumOfForcesInLocalX(node, x_bar, other_bar, other_force)
            else:
                x_force = SumOfForcesInLocalX(node, x_bar)
            x_bar.SetAxialLoad(x_force)
            x_bar.is_computed = True
            progress = True
        if not progress:
            sys.exit("No viable node found; method of joints cannot continue")
