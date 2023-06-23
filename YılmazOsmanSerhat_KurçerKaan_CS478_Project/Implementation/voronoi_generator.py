import math
import random
import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from IPython.display import clear_output


# RANDOMIZED INCREMENTAL ALGORITHM ----------------------------------------------------------

class DelaunayRI:

    def __init__(self, center=(0, 0), radius=9999):
        center = np.asarray(center)
        self.coords = [center+radius*np.array((-1, -1)),
                       center+radius*np.array((+1, -1)),
                       center+radius*np.array((+1, +1)),
                       center+radius*np.array((-1, +1))]
        self.triangles = {}
        self.circles = {}

        T1 = (0, 1, 3)
        T2 = (2, 3, 1)
        self.triangles[T1] = [T2, None, None]
        self.triangles[T2] = [T1, None, None]

        for t in self.triangles:
            self.circles[t] = self.circumcenter(t)

    def circumcenter(self, tri):
        pts = np.asarray([self.coords[v] for v in tri])
        pts2 = np.dot(pts, pts.T)
        A = np.bmat([[2 * pts2, [[1],
                                 [1],
                                 [1]]],
                     [[[1, 1, 1, 0]]]])

        b = np.hstack((np.sum(pts * pts, axis=1), [1]))
        x = np.linalg.solve(A, b)
        bary_coords = x[:-1]
        center = np.dot(bary_coords, pts)

        radius = np.sum(np.square(pts[0] - center))
        return (center, radius)

    def inCircleFast(self, tri, p):
        center, radius = self.circles[tri]
        return np.sum(np.square(center - p)) <= radius

    def addPoint(self, p):
        p = np.asarray(p)
        idx = len(self.coords)
        self.coords.append(p)

        # Search the triangle(s) whose circumcircle contains p
        bad_triangles = []
        for T in self.triangles:
            if self.inCircleFast(T, p):
                bad_triangles.append(T)

        boundary = []
        T = bad_triangles[0]
        edge = 0
        while True:
            tri_op = self.triangles[T][edge]
            if tri_op not in bad_triangles:
                boundary.append((T[(edge+1) % 3], T[(edge-1) % 3], tri_op))
                edge = (edge + 1) % 3

                if boundary[0][0] == boundary[-1][1]:
                    break
            else:
                edge = (self.triangles[tri_op].index(T) + 1) % 3
                T = tri_op

        for T in bad_triangles:
            del self.triangles[T]
            del self.circles[T]

        new_triangles = []
        for (e0, e1, tri_op) in boundary:
            T = (idx, e0, e1)

            self.circles[T] = self.circumcenter(T)
            self.triangles[T] = [tri_op, None, None]

            if tri_op:
                for i, neigh in enumerate(self.triangles[tri_op]):
                    if neigh:
                        if e1 in neigh and e0 in neigh:
                            self.triangles[tri_op][i] = T

            new_triangles.append(T)

        N = len(new_triangles)
        for i, T in enumerate(new_triangles):
            self.triangles[T][1] = new_triangles[(i+1) % N]
            self.triangles[T][2] = new_triangles[(i-1) % N]

        return self.exportDT()

    def exportTriangles(self):
        return [(a-4, b-4, c-4)
                for (a, b, c) in self.triangles if a > 3 and b > 3 and c > 3]

    def exportDT(self):
        # Filter out coordinates in the extended BBox
        coord = self.coords[4:]

        # Filter out triangles with any vertex in the extended BBox
        tris = [(a-4, b-4, c-4)
                for (a, b, c) in self.triangles if a > 3 and b > 3 and c > 3]
        return coord, tris

    def exportVoronoiRegions(self):
        useVertex = {i: [] for i in range(len(self.coords))}
        vor_coors = []
        index = {}
        for tidx, (a, b, c) in enumerate(sorted(self.triangles)):
            vor_coors.append(self.circles[(a, b, c)][0])

            useVertex[a] += [(b, c, a)]
            useVertex[b] += [(c, a, b)]
            useVertex[c] += [(a, b, c)]

            index[(a, b, c)] = tidx
            index[(c, a, b)] = tidx
            index[(b, c, a)] = tidx

        regions = {}
        for i in range(4, len(self.coords)):
            v = useVertex[i][0][0]
            r = []
            for _ in range(len(useVertex[i])):
                t = [t for t in useVertex[i] if t[0] == v][0]
                r.append(index[t])
                v = t[1]
            regions[i-4] = r

        return vor_coors, regions

def plot_triangles_RI( fig, ax, T):
    for triangle in T:
        if triangle[0] is None or triangle[1] is None or triangle[2] is None:
            continue
        polygon = Polygon(triangle, closed=True, fill=None, edgecolor='red')
        ax.add_patch(polygon)

    ax.set_aspect('equal')
    ax.autoscale_view()
    plt.show()

def plot_voronoi_RI( fig, ax, T):
    for region in T:
        polygon = Polygon(region, closed=True, fill=None, edgecolor='black')
        ax.add_patch(polygon)

    ax.set_aspect('equal')
    plt.show()


def get_triangles_RI(coordinates, triangles):
    all_triangles = []
    for triangle in triangles:
        points = [(coordinates[idx][0], coordinates[idx][1]) for idx in triangle]
        all_triangles.append(points)

    return all_triangles

def get_regions_RI(coordinates, regions):
    all_regions = []
    for region in regions.values():
        points = [(coordinates[idx][0], coordinates[idx][1]) for idx in region]
        all_regions.append(points)

    return all_regions


def calculate_voronoi_RI(points):
    fig, ax = plt.subplots()
    if len(points) >= 3:
        dt = DelaunayRI()
        for point in points:
            coordinates, triangles = dt.addPoint(point)
            triangles = get_triangles_RI(coordinates, triangles)

            plt.cla()
            x,y = zip(*points)
            plt.scatter(x,y)
            plot_triangles_RI(fig, ax, triangles)

            clear_output(wait = True)
            plt.pause(0.3)

        vor_coors, regions = dt.exportVoronoiRegions()
        voronoi_regions = get_regions_RI(vor_coors, regions)
        plot_voronoi_RI(fig, ax, voronoi_regions)

# FLIPPING ALGORITHM ----------------------------------------------------------

def findCircumCenter(P, Q, R):

    def lineFromPoints(P, Q, a, b, c):
        a = Q[1] - P[1]
        b = P[0] - Q[0]
        c = a * (P[0]) + b * (P[1])
        return a, b, c

    def perpendicularBisectorFromLine(P, Q, a, b, c):
        mid_point = [(P[0] + Q[0])//2, (P[1] + Q[1])//2]

        # c = -bx + ay
        c = -b * (mid_point[0]) + a * (mid_point[1])
        temp = a
        a = -b
        b = temp
        return a, b, c

    def lineLineIntersection(a1, b1, c1, a2, b2, c2):
        determinant = a1 * b2 - a2 * b1
        if (determinant == 0):
            return [(10.0)**19, (10.0)**19]
        else:
            x = (b2 * c1 - b1 * c2)//determinant
            y = (a1 * c2 - a2 * c1)//determinant
            return [x, y]

    if P is None or Q is None or R is None:
        return Q

    a, b, c = 0.0, 0.0, 0.0
    a, b, c = lineFromPoints(P, Q, a, b, c)

    e, f, g = 0.0, 0.0, 0.0
    e, f, g = lineFromPoints(Q, R, e, f, g)

    a, b, c = perpendicularBisectorFromLine(P, Q, a, b, c)
    e, f, g = perpendicularBisectorFromLine(Q, R, e, f, g)

    circumcenter = lineLineIntersection(a, b, c, e, f, g)

    if (circumcenter[0] == (10.0)**19 and circumcenter[1] == (10.0)**19):
        return Q
    else:
        return circumcenter

def is_illegal(pr, edge, pk):

    if pk is None:
        return False

    circumcenter = findCircumCenter(pr, edge[0], edge[1])
    if circumcenter is None:
        return False

    return math.dist(circumcenter, pr) >= math.dist(circumcenter, pk)

def find_adjacent_triangle(pr, pi, pj, T):
    pl = find_third_vertex(pr, pi, pj, T)
    if pl is not None:
        return pl, pi, pj
    else:
        return None, pi, pj

def flip_edge(edge1, edge2, T):
    for triangle in T:
        if edge1[0] in triangle and edge1[1] in triangle:
            # Remove the triangle from T
            T.remove(triangle)

    # Create the new triangles with the flipped edge
    triangle1 = [edge1[0], edge2[0], edge2[1]]
    triangle2 = [edge1[1], edge2[0], edge2[1]]

    # Add the new triangles to T
    T.append(triangle1)
    T.append(triangle2)

    return T

def find_initial_triangle(P):
    min_x = min(P, key=lambda p: p[0])[0]
    max_x = max(P, key=lambda p: p[0])[0]
    min_y = min(P, key=lambda p: p[1])[1]
    max_y = max(P, key=lambda p: p[1])[1]

    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2

    size = max(max_x - min_x, max_y - min_y)

    # Create the initial triangle vertices
    p_minus = [(int(center_x - 2 * size), int(center_y - size)),
               (int(center_x), int(center_y + 2 * size)),
               (int(center_x + 2 * size), int(center_y - size))]

    return p_minus

def find_triangle_containing_point(pr, T):
    for triangle in T:
        pi, pj, pk = triangle
        # Check if the point pr is inside the triangle
        if is_inside_triangle(pr, pi, pj, pk):
            return pi, pj, pk
    return None, None, None

def is_on_edge(point, edge):
    p1, p2 = edge
    if p1 is None or p2 is None:
        return False
    x1, y1 = p1
    x2, y2 = p2
    x, y = point

    slope = (y2 - y1) / (x2 - x1) if x2 != x1 else None

    if slope is not None:
        # Check if the point lies on the line defined by the edge
        expected_y = y1 + slope * (x - x1)
        if y == expected_y and min(x1, x2) <= x <= max(x1, x2):
            return True
    else:
        # Edge is vertical, check if the point has the same x-coordinate
        if x == x1 and min(y1, y2) <= y <= max(y1, y2):
            return True

    return False

def is_inside_triangle(pr, pi, pj, pk):
    x, y = pr
    if pi is None or pj is None or pk is None:
        return False
    x1, y1 = pi
    x2, y2 = pj
    x3, y3 = pk

    def calculate_triangle_area(x1, y1, x2, y2, x3, y3):
        return abs( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)) / 2

    # Calculate area of triangle ABC
    A = calculate_triangle_area(x1, y1, x2, y2, x3, y3)
    # Calculate area of triangle PBC
    A1 = calculate_triangle_area(x, y, x2, y2, x3, y3)
    # Calculate area of triangle PAC
    A2 = calculate_triangle_area(x1, y1, x, y, x3, y3)
    # Calculate area of triangle PAB
    A3 = calculate_triangle_area(x1, y1, x2, y2, x, y)

    return A == A1 + A2 + A3

def split_triangle(pr, pi, pj, pk, T):
    # Create the three new triangles by adding edges from pr to pi, pj, and pk
    triangle1 = [pi, pj, pr]
    triangle2 = [pj, pk, pr]
    triangle3 = [pk, pi, pr]

    # Replace the original triangle with the new triangles in T
    for triangle in T:
        if pi in triangle and pj in triangle and pk in triangle:
            T.remove(triangle)

    T.append(triangle1)
    T.append(triangle2)
    T.append(triangle3)

    return T

def find_third_vertex(pr, pi, pj, T):
    for triangle in T:
        if pi in triangle and pj in triangle and pr not in triangle:
            adjacent_vertices = [v for v in triangle if v not in (pi, pj)]
            if len(adjacent_vertices) > 0:
                return adjacent_vertices[0]

    return None

def split_triangles_with_edge(pr, pi, pj, pk, pl, T):
    triangle1 = [pr, pi, pl]
    triangle2 = [pr, pl, pj]
    triangle3 = [pr, pj, pk]
    triangle4 = [pr, pk, pi]

    # Replace the original triangle with the new triangles in T
    for triangle in T:
        if pi in triangle and pj in triangle and pk in triangle:
            T.remove(triangle) #T.remove([pi, pj, pk])
    for triangle in T:
        if pi in triangle and pj in triangle and pl in triangle:
            T.remove(triangle) #T.remove([pi, pj, pl])

    T.append(triangle1)
    T.append(triangle2)
    T.append(triangle3)
    T.append(triangle4)

    return T

def remove_initial_triangle(p_minus, T):
    filtered_T = []
    for triangle in T:
        if any(point in triangle for point in p_minus):
            continue
        filtered_T.append(triangle)
    return filtered_T

def generate_points(distribution, num_points):
    points = []

    if distribution == 'random':
        for _ in range(num_points):
            x = random.randint(0, 800)
            y = random.randint(0, 600)
            points.append((x, y))

    elif distribution == 'uniform':
        x_coords = np.random.uniform(0, 800, num_points)
        y_coords = np.random.uniform(0, 600, num_points)
        points = [(int(x), int(y)) for x, y in zip(x_coords, y_coords)]

    elif distribution == 'gaussian':
        x_mean, x_std = 400, 200
        y_mean, y_std = 300, 150
        x_coords = np.random.normal(x_mean, x_std, num_points)
        y_coords = np.random.normal(y_mean, y_std, num_points)
        points = [(int(x), int(y)) for x, y in zip(x_coords, y_coords)]

    return points


def LEGALIZEEDGE(pr, edge, T):
    pk, pi, pj = find_adjacent_triangle(pr, edge[0], edge[1], T)
    if pk is not None and is_illegal(pr, edge, pk):
        T = flip_edge((pi, pj), (pr, pk), T)
        LEGALIZEEDGE(pr, (pi, pk), T)
        LEGALIZEEDGE(pr, (pk, pj), T)

    return T

def DELAUNAYTRIANGULATION(P):
    fig, ax = plt.subplots()
    p_minus = find_initial_triangle(P)
    T = [p_minus]
    print(T)
    plot_triangles(fig, ax, T)
    P = random.sample(P, len(P))
    for r in range(len(P)):
        pr = P[r]
        pi, pj, pk = find_triangle_containing_point(pr, T)
        if is_on_edge(pr, (pi, pj)) or is_on_edge(pr, (pk, pj)) or is_on_edge(pr, (pi, pk)):
            pl = find_third_vertex(pr, pi, pj, T)
            T = split_triangles_with_edge(pr, pi, pj, pk, pl, T)
            T = LEGALIZEEDGE(pr, (pi, pl), T)
            T = LEGALIZEEDGE(pr, (pl, pj), T)
            T = LEGALIZEEDGE(pr, (pj, pk), T)
            T = LEGALIZEEDGE(pr, (pk, pi), T)
        else:
            T = split_triangle(pr, pi, pj, pk, T)
            T = LEGALIZEEDGE(pr, (pi, pj), T)
            T = LEGALIZEEDGE(pr, (pj, pk), T)
            T = LEGALIZEEDGE(pr, (pk, pi), T)


        plt.cla()
        x,y = zip(*P)
        plt.scatter(x,y)
        plot_triangles(fig, ax, T)
        clear_output(wait = True)
        plt.pause(0.3)

    T = remove_initial_triangle(p_minus, T)
    return T

def plot_triangles(fig, ax, T):
    for triangle in T:
        if triangle[0] is None or triangle[1] is None or triangle[2] is None:
            continue
        polygon = Polygon(triangle, closed=True, fill=None, edgecolor='red')
        ax.add_patch(polygon)

    ax.set_aspect('equal')
    ax.autoscale_view()
    plt.show()

def plot_voronoi(points, edges, T):
    """
    Plots the Voronoi diagram using the given set of points and edges.
    Args:
        points: List of points (x, y) representing the Voronoi sites.
        edges: List of tuples (start, end) representing the Voronoi edges.
    """
    x_coords, y_coords = zip(*points)
    fig, ax = plt.subplots()

    for edge in edges:
        start, end = edge
        x_values = [start[0], end[0]]
        y_values = [start[1], end[1]]
        ax.plot(x_values, y_values, 'k-', linewidth=1)

    ax.plot(x_coords, y_coords, 'ro')
    ax.set_aspect('equal')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('Voronoi Diagram')


    for triangle in T:
        if triangle[0] is None or triangle[1] is None or triangle[2] is None:
            continue
        polygon = Polygon(triangle, closed=True, fill=None, edgecolor='red')
        ax.add_patch(polygon)

    plt.show()

def VORONOI(T):
    edges = []
    points = []

    for triangle in T:
        a,b,c = triangle
        if a is None or b is None or c is None:
            continue

        c1 = findCircumCenter(a, b, c)
        if c1 is None:
            c1 = a
        points.append(c1)

        flag = True
        for triangle2 in T:
            c2 = findCircumCenter(triangle2[0], triangle2[1], triangle2[2])
            if c2 is None:
                c2 = triangle2[0]

            if a in triangle2 and b in triangle2 and c not in triangle2:
                edges.append((c1,c2))
                flag = False

        if flag:
            x_mid = (a[0] + b[0]) / 2.0
            y_mid = (a[1] + b[1]) / 2.0
            mid = [x_mid, y_mid]

            edges.append((c1,mid))


        flag = True
        for triangle2 in T:
            c2 = findCircumCenter(triangle2[0], triangle2[1], triangle2[2])
            if b in triangle2 and c in triangle2 and a not in triangle2:
                edges.append((c1,c2))
                flag = False

        if flag:
            x_mid = (b[0] + c[0]) / 2.0
            y_mid = (c[1] + c[1]) / 2.0
            mid = [x_mid, y_mid]

            edges.append((c1,mid))


        flag = True
        for triangle2 in T:
            c2 = findCircumCenter(triangle2[0], triangle2[1], triangle2[2])
            if a in triangle2 and c in triangle2 and b not in triangle2:
                edges.append((c1,c2))
                flag = False

        if flag:
            x_mid = (a[0] + c[0]) / 2.0
            y_mid = (a[1] + c[1]) / 2.0
            mid = [x_mid, y_mid]

            edges.append((c1,mid))


    return points, edges


print("Please enter the number of points")
point_count = int(input())

print("Please enter the distribution you want to use\n1 for uniform\n2 for random\n3 for gaussian")
dist = int(input())
if dist == 1:
    dist = 'uniform'
elif dist == 2:
    dist = 'random'
elif dist == 3:
    dist = 'gaussian'

P = generate_points(dist, point_count)


print("Please enter the algorithm you want to use\n1 for randomized incremental\n2 for flipping")
algo = int(input())
if algo == 2:
    tris = DELAUNAYTRIANGULATION(P)
    p, e = VORONOI(tris)
    plot_voronoi(P, e, tris)
elif algo == 1:
    calculate_voronoi_RI(P)
