import tkinter as tk
from delaunayRandomIncremental import DelaunayRI

import tkinter as tk
from delaunayRandomIncremental import DelaunayRI
import random
import numpy as np

class VoronoiDiagramGUI:
    def __init__(self, root):
        self.voronoi_edges = []  # Initialize the voronoi_edges attribute
        self.voronoi_coors = {}  # Initialize the voronoi_coordinates attribute
        self.root = root
        self.canvas = tk.Canvas(root, width=800, height=600)
        self.canvas.pack()
        self.points = []
        self.triangles = []
        self.voronoi_edges = []

        self.canvas.bind('<Button-1>', self.add_point)
        self.canvas.bind('<Button-3>', self.calculate_voronoi)
        self.canvas.bind('<Button-2>', self.clear)

    def add_point(self, event):
        self.points.append((event.x, event.y))
        self.canvas.create_oval(event.x - 2, event.y - 2, event.x + 2, event.y + 2, fill='red')

    def generate_points(self, distribution, num_points):
        self.points = []

        if distribution == 'random':
            for _ in range(num_points):
                x = random.randint(0, 800)
                y = random.randint(0, 600)
                self.points.append((x, y))

        elif distribution == 'uniform':
            x_coords = np.random.uniform(0, 800, num_points)
            y_coords = np.random.uniform(0, 600, num_points)
            self.points = [(int(x), int(y)) for x, y in zip(x_coords, y_coords)]

        elif distribution == 'gaussian':
            x_mean, x_std = 400, 200
            y_mean, y_std = 300, 150
            x_coords = np.random.normal(x_mean, x_std, num_points)
            y_coords = np.random.normal(y_mean, y_std, num_points)
            self.points = [(int(x), int(y)) for x, y in zip(x_coords, y_coords)]

        self.draw_points()

    def calculate_voronoi(self, event):
        if len(self.points) >= 3:
            dt = DelaunayRI()
            for point in self.points:
                dt.addPoint(point)
            self.triangles = dt.exportDT()[1]
            self.voronoi_edges = dt.exportVoronoiRegions()[0]

            self.draw_voronoi()

    def draw_points(self):
        self.canvas.delete('points')
        for point in self.points:
            x, y = point
            self.canvas.create_oval(x - 2, y - 2, x + 2, y + 2, fill='red', tags='points')

    #def draw_voronoi(self):
    #    for edge in self.voronoi_edges:
    #        p1, p2 = edge[0], edge[1]
    #       self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill='blue')

    def draw_voronoi(self):
        min_x = min(self.vor_coors, key=lambda coor: coor[0])[0]
        max_x = max(self.vor_coors, key=lambda coor: coor[0])[0]
        min_y = min(self.vor_coors, key=lambda coor: coor[1])[1]
        max_y = max(self.vor_coors, key=lambda coor: coor[1])[1]

        for region in self.regions.values():
            points = [(self.vor_coors[idx][0], self.vor_coors[idx][1]) for idx in region]
            self.canvas.create_polygon(points, outline='black', fill='white')

        # Scale the canvas to fit the Voronoi diagram
        scale = 0.9 * min(800 / (max_x - min_x), 600 / (max_y - min_y))
        self.canvas.scale('all', 0, 0, scale, scale)

        self.root.mainloop()

    def clear(self, event):
        self.canvas.delete('all')
        self.points = []
        self.triangles = []
        self.voronoi_edges = []

if __name__ == '__main__':
    root = tk.Tk()
    app = VoronoiDiagramGUI(root)
    app.generate_points('random', 50)  # Example usage: generate 50 random points
    root.mainloop()
