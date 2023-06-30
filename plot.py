import json 
import numpy as np
import matplotlib.pyplot as plt 
import time 
class potential_flow_plotter:

    def __init__(self, json_file_conditions, json_file_data, i):
        self.json_file_conditions = json_file_conditions 
        self.json_file_data = json_file_data
        self.i = i  
        self.load_json_conditions()
        self.load_json_data()
        
    def load_json_conditions(self):
        with open(self.json_file_conditions, 'r') as json_handle:
            input_vals = json.load(json_handle)
            self.R = input_vals["geometry"]['cylinder_radius']
            self.x_start = input_vals["plot"]["x_start"]
            self.x_leading_edge = input_vals["plot"]["x_leading_edge"]
            self.x_trailing_edge = input_vals["plot"]["x_trailing_edge"]
            self.x_lower_limit = input_vals["plot"]["x_lower_limit"]
            self.x_upper_limit = input_vals["plot"]["x_upper_limit"]
            self.x_cent = input_vals["plot"]["x_cent"]
            self.delta_s = input_vals["plot"]["delta_s"]  # delta s is the time step value. 
            self.n_lines = input_vals["plot"]["n_lines"]
            self.delta_y = input_vals["plot"]["delta_y"]
            self.vel_inf = input_vals["Operating"]["vel_inf"]
            self.vortex_strength = input_vals["Operating"]["vortex_strength"]
            self.alpha = np.radians(input_vals["Operating"]["alpha[deg]"])

    def load_json_data_type(self):
        with open(self.json_file_data, 'r') as json_handle:
            input_vals = json.load(json_handle)
            self.plot_type = input_vals['plot_points'][self.i]

    def load_json_data(self):
        # plot = self.plot_type
        self.load_json_data_type()
        with open(self.plot_type, 'r') as csv_handle:
            points = [list(map(float, line.strip().split())) for line in csv_handle]
        self.plot_data = np.array(points)
    
    def plotter(self):
        """This function takes in the geometry characteristics specified in main and plots the geometry. It also plots the streamlines at delta_y increments above and below the stagnation lines."""
        plt.scatter(self.plot_data[:, 0], self.plot_data[:, 1], s = 2)  # Set the size of scatter points to 10
        plt.title("Potential Flow (Cylinder)")
        plt.xlim(self.x_lower_limit, self.x_upper_limit)
        # plt.xlim(self.x_lower_limit, self.x_upper_limit)
        plt.ylim(-self.R*2.5, self.R*2.5)        
        plt.gca().set_aspect('equal')

    def run(self):
        """This function takes all the above functions and runs them"""
        self.load_json_conditions()
        self.load_json_data_type()
        self.load_json_data()
        self.plotter()

if __name__ == "__main__":
    if __name__ == "__main__":
        time_1 = time.perf_counter()
        
        plot_object_geom = potential_flow_plotter('Streamline_A1.json','data.json', 0)
        plot_object_geom.run()

        plot_object_stags = potential_flow_plotter('Streamline_A1.json','data.json', 1)
        plot_object_stags.run()

        plot_object_upper_streams = potential_flow_plotter('Streamline_A1.json','data.json', 2)
        plot_object_upper_streams.run()

        plot_object_lower_streams = potential_flow_plotter('Streamline_A1.json','data.json', 3)
        plot_object_lower_streams.run()

        time_2 = time.perf_counter()
        print(f"Plot_generate_time:{time_2-time_1:0.4f} seconds")
        
        plt.show()

        







    
