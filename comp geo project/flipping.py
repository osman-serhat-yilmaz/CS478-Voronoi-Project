from scipy.spatial import Delaunay

def flip_edge(triangles, edges, tri_index, edge_index):
    # Get the indices of the triangles adjacent to the edge
    tri1_index, tri2_index = triangles[tri_index][edge_index], triangles[tri_index][(edge_index + 1) % 3]

    # Get the vertices of the edge
    vertex1, vertex2 = edges[(tri_index, edge_index)]

    # Find the third vertex of each triangle
    tri1_vertex = [v for v in triangles[tri1_index] if v != vertex1 and v != vertex2][0]
    tri2_vertex = [v for v in triangles[tri2_index] if v != vertex1 and v != vertex2][0]

    # Check if flipping the edge would result in a valid Delaunay triangulation
    if is_edge_delaunay(triangles, edges, tri1_index, tri1_vertex, tri2_index, tri2_vertex):
        # Remove the old edges
        del edges[(tri_index, edge_index)]
        del edges[(tri1_index, (edges[(tri1_index, 0)] == vertex2) + (edges[(tri1_index, 1)] == vertex2) * 2)]

        # Add the new edges
        edges[(tri_index, (edge_index + 2) % 3)] = tri2_vertex
        edges[(tri1_index, (edges[(tri1_index, 0)] == tri1_vertex) + (edges[(tri1_index, 1)] == tri1_vertex) * 2)] = vertex1

        # Update the triangle connectivity
        triangles[tri_index][edge_index] = tri2_index
        triangles[tri1_index][(edges[(tri1_index, 0)] == tri1_vertex) + (edges[(tri1_index, 1)] == tri1_vertex) * 2] = tri_index

        return True

    return False

def is_edge_delaunay(triangles, edges, tri1_index, tri1_vertex, tri2_index, tri2_vertex):
    # Check if the edge is shared by the two triangles
    if tri1_index not in triangles[tri2_index]:
        return False

    # Get the other vertices of the triangles
    tri1_other = [v for v in triangles[tri1_index] if v != tri1_vertex]
    tri2_other = [v for v in triangles[tri2_index] if v != tri2_vertex]

    # Get the coordinates of the vertices
    v1, v2, v3, v4 = tri1_vertex, tri1_other[0], tri2_vertex, tri2_other[0]

    # Check if the circumcircle is empty or contains both vertices of the other triangle
    return not point_in_circumcircle(v1, v2, v3, v4)

def point_in_circumcircle(v1, v2, v3, v4):
    # Calculate the determinant of the 4x4 matrix
    det = (v1[0] - v4[0]) * (v2[1] - v4[1]) * (v3[2] - v4[2]) + \
          (v2[0] - v4[0]) * (v3[1] - v4[1]) * (v1[2] - v4[2]) + \
          (v3[0] - v4[0]) * (v1[1] - v4[1]) * (v2[2] - v4[2]) - \
          (v1[2] - v4[2]) * (v2[1] - v4[1]) * (v3[0] - v4[0]) - \
          (v2[2] - v4[2]) * (v3[1] - v4[1]) * (v1[0] - v4[0]) - \
          (v3[2] - v4[2]) * (v1[1] - v4[1]) * (v2[0] - v4[0])

    return det > 0

def lexicographic_to_delaunay_triangulation(points, lex_triangles, lex_edges):
    # Create an empty list to store the triangles
    triangles = lex_triangles

    # Create a dictionary to store the edges
    it = iter(lex_edges)
    lex_edges_dct = dict(zip(it, it))
    edges = lex_edges_dct



    # Sort the points lexicographically
    sorted_points = sorted(points)

    # Process each additional point
    for i in range(4, len(sorted_points)):
        point = sorted_points[i]

        # Find the triangle containing the point
        for j, triangle in enumerate(triangles):
            vertex1, vertex2, vertex3 = triangle

            # Check if the point lies inside the triangle
            if point_in_triangle(point, sorted_points[vertex1], sorted_points[vertex2], sorted_points[vertex3]):
                # Create three new triangles by connecting the point to the triangle's vertices
                triangle1 = [vertex1, vertex2, i]
                triangle2 = [vertex2, vertex3, i]
                triangle3 = [vertex3, vertex1, i]

                triangles[j] = triangle1

                triangles.append(triangle2)
                triangles.append(triangle3)

                edges[(j, 0)] = len(triangles) - 2
                edges[(j, 1)] = len(triangles) - 1
                edges[(j, 2)] = len(triangles) - 3

                edges[(len(triangles) - 2, 0)] = len(triangles) - 1
                edges[(len(triangles) - 2, 1)] = len(triangles) - 3
                edges[(len(triangles) - 2, 2)] = j

                edges[(len(triangles) - 3, 0)] = len(triangles) - 2
                edges[(len(triangles) - 3, 1)] = j
                edges[(len(triangles) - 3, 2)] = len(triangles) - 1

                # Flip any edges that are not Delaunay
                flipped = True

                while flipped:
                    flipped = False

                    for tri_index, triangle in enumerate(triangles):
                        for edge_index in range(3):
                            if (tri_index, edge_index) in edges and flip_edge(triangles, edges, tri_index, edge_index):
                                flipped = True

                break

    # Return the resulting Delaunay triangulation
    return triangles

def point_in_triangle(point, vertex1, vertex2, vertex3):
    d1 = (point[0] - vertex3[0]) * (vertex1[1] - vertex3[1]) - (point[1] - vertex3[1]) * (vertex1[0] - vertex3[0])
    d2 = (point[0] - vertex1[0]) * (vertex2[1] - vertex1[1]) - (point[1] - vertex1[1]) * (vertex2[0] - vertex1[0])
    d3 = (point[0] - vertex2[0]) * (vertex3[1] - vertex2[1]) - (point[1] - vertex2[1]) * (vertex3[0] - vertex2[0])

    return (d1 >= 0 and d2 >= 0 and d3 >= 0) or (d1 <= 0 and d2 <= 0 and d3 <= 0)

import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

def plot_triangulation(points, triangles, title):
    # Create a triangulation object
    triangulation = Triangulation([point[0] for point in points], [point[1] for point in points], triangles)

    # Plot the triangulation
    plt.figure()
    plt.title(title)
    plt.triplot(triangulation, 'k-', linewidth=0.5)
    plt.plot([point[0] for point in points], [point[1] for point in points], 'ro', markersize=5)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.axis('equal')
    plt.show()

from scipy.spatial import ConvexHull

def lexicographic_triangulation(points):
    # Perform Delaunay triangulation
    triangulation = Delaunay(points)

    # Get the sorted triangles
    sorted_triangles = sorted(triangulation.simplices.tolist())

    # Get the edges from the triangulation
    edges = triangulation.neighbors.tolist()

    return sorted_triangles, edges


# Example usage
points = [(1, 0), (0, 1), (0.5, 0.5), (0.2, 0.4), (0.8, 0.9), (0.1,0.8), (2,1)]
triangles, edges = lexicographic_triangulation(points)

# Initial lexicographic triangulation
plot_triangulation(points, triangles, 'Initial Lexicographic Triangulation')

# Convert to Delaunay triangulation
delaunay_triangulation = lexicographic_to_delaunay_triangulation(points, triangles, edges)
plot_triangulation(points, delaunay_triangulation, 'Final Delaunay Triangulation')
