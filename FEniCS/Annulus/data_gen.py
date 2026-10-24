# %%

# %%
import pyvista as pv
import pandas as pd
import numpy as np
import os
import meshio


# %%

# %%
e = np.e 
pi=np.pi
KA = 0.001
AK = 1000
alpha0_set = np.linspace(0.5, 1, 6) * np.pi
phi_set = np.array([0.52])

rho_set=np.array([1/e, 1/e**2, 1/np.sqrt(e)])
r_set = np.sqrt(rho_set)
r_set = np.append(r_set, rho_set**0.75)

# %%

# %%
def calculate_delta(mesh):   
    S = mesh.point_data['S_order_param']
    grad_qxx = mesh.compute_derivative(scalars='Qxx').point_data["gradient"]
    grad_qxy = mesh.compute_derivative(scalars='Qxy').point_data["gradient"]
    
    # Corrected formula to ensure the entire numerator is divided by S**2
    delta = ((grad_qxy[:, 1] * grad_qxx[:, 0]) - (grad_qxy[:, 0] * grad_qxx[:, 1])) / (S * S)     
    
    return delta

def calculate_free_energy(mesh):
    S = mesh.point_data['S_order_param']
    grad_qxx = mesh.compute_derivative(scalars='Qxx').point_data["gradient"]
    grad_qxy = mesh.compute_derivative(scalars='Qxy').point_data["gradient"]

    free_energy = 0.5*(grad_qxx[:, 0]*grad_qxx[:, 0] + grad_qxx[:, 1]* grad_qxx[:, 1] + grad_qxy[:, 0]* grad_qxy[:, 0]+ grad_qxy[:, 1]* grad_qxy[:, 1])

    return free_energy

def calculate_total_free_energy(mesh):
    # 1. Compute derivatives (Note: these return NEW mesh objects)
    # We can chain .point_data["gradient"] to get the array immediately
    grad_qxx = mesh.compute_derivative(scalars='Qxx').point_data["gradient"]
    grad_qxy = mesh.compute_derivative(scalars='Qxy').point_data["gradient"]

    # 2. Calculate the local free energy density (NumPy array)
    free_energy_density = 0.5 * (
        grad_qxx[:, 0]**2 + grad_qxx[:, 1]**2 + 
        grad_qxy[:, 0]**2 + grad_qxy[:, 1]**2
    )

    # 3. Create a lightweight copy of the mesh geometry to isolate this calculation
    # (Following the single-quantity optimization from earlier)
    energy_mesh = mesh.copy(deep=False)
    energy_mesh.clear_point_data()
    energy_mesh.clear_cell_data()
    
    # 4. Attach the array to the geometry
    energy_mesh.point_data['F_density'] = free_energy_density

    # 5. Integrate over the mesh to get the total energy
    integrated_result = energy_mesh.integrate_data()
    
    # 6. Extract the single scalar total
    total_free_energy = integrated_result.point_data['F_density'][0]
    
    return total_free_energy





# Initialize a list to accumulate dataframes from all files
all_data = []

for rho in rho_set[:]:
    for phi in phi_set:
        for r in r_set:
            for alpha0 in alpha0_set:
                xdmf_path = "./solutions/defects_solution_alpha0{0:.2f}_phi{1:.2f}_r{2:.2f}_KA{3:.4f}_rho{4:.2f}.xdmf".format(alpha0, phi, r, KA, rho)

                if not os.path.exists(xdmf_path):
                    print(f"File missing, skipping: {xdmf_path}")
                    continue

                print(f"Opening file: {xdmf_path}")
                
                t_list = []
                max_x_list, max_y_list = [], []
                min_x_list, min_y_list = [], []
                max_delta_list, min_delta_list = [], []  # Added to track max/min delta scalar values
                total_free_energy = []

                with meshio.xdmf.TimeSeriesReader(xdmf_path) as reader:
                    points, cells = reader.read_points_cells()
                    
                    # Ensure points are 3D for PyVista
                    if len(points.shape) == 2:
                        points = np.column_stack([points, np.zeros(len(points))])
                        
                    num_steps = reader.num_steps
                    print(f"Detected {num_steps} time steps.")

                    if isinstance(cells, list):
                        triangle_block = next((block for block in cells if block.type == "triangle"), cells[0])
                        cell_connectivity = triangle_block.data
                    else:
                        cell_connectivity = cells.data

                    cell_types = np.full(len(cell_connectivity), pv.CellType.TRIANGLE)
                    pv_cells = np.column_stack([np.full(len(cell_connectivity), 3), cell_connectivity])
                    base_mesh = pv.UnstructuredGrid(pv_cells.ravel(), cell_types, points)

                    for i in range(0, num_steps, 10):
                        t, point_data, _ = reader.read_data(i)
                        
                        mesh = base_mesh.copy(deep=False)
                        for name, data in point_data.items():
                            mesh.point_data[name] = data

                        delta = calculate_delta(mesh)
                        
                        
                        idx_max = np.argmax(delta)
                        idx_min = np.argmin(delta)
                        
                        val_max = delta[idx_max]
                        val_min = delta[idx_min]
                        
                        # Check threshold conditions using the calculated scalar limits
                        if val_max > 10 and val_min < -10:
                            t_list.append(t)
                            
                            max_x_list.append(mesh.points[idx_max, 0])
                            max_y_list.append(mesh.points[idx_max, 1])
                            
                            min_x_list.append(mesh.points[idx_min, 0])
                            min_y_list.append(mesh.points[idx_min, 1])
                            
                            max_delta_list.append(val_max)  # Store maximum delta value
                            min_delta_list.append(val_min)  # Store minimum delta value
                            total_F = calculate_total_free_energy(mesh)
                            total_free_energy.append(total_F)

                        if i % 100 == 0:
                            print(f"Frame {i} processed.")

                # Convert data to a pandas DataFrame if steps were captured
                if len(t_list) > 0:
                    t_arr = np.array(t_list)
                    max_x_arr = np.array(max_x_list)
                    max_y_arr = np.array(max_y_list)
                    min_x_arr = np.array(min_x_list)
                    min_y_arr = np.array(min_y_list)
                    max_delta_arr = np.array(max_delta_list)
                    min_delta_arr = np.array(min_delta_list)
                    total_F_arr = np.array(total_free_energy)
                    
                    # Compute polar coordinates
                    r_max = np.sqrt(max_x_arr**2 + max_y_arr**2)
                    r_min = np.sqrt(min_x_arr**2 + min_y_arr**2)
                    phi_max = np.arctan2(max_y_arr, max_x_arr)
                    phi_min = np.arctan2(min_y_arr, min_x_arr)
                    

                    #Adding rho array

                    rho_array =np.ones_like(r_max)*rho
                    # Create a DataFrame for the current file parameters
                    df_file = pd.DataFrame({
                        'alpha0': alpha0,
                        'phi_param': phi,
                        'r_param': r,
                        'KA': KA,
                        't': t_arr,
                        'max_x': max_x_arr,
                        'max_y': max_y_arr,
                        'min_x': min_x_arr,
                        'min_y': min_y_arr,
                        'r_max': r_max,
                        'r_min': r_min,
                        'phi_max': phi_max,
                        'phi_min': phi_min,
                        'max_delta': max_delta_arr,  # New column
                        'min_delta': min_delta_arr,
                        'rho': rho_array,
                        'total_free_energy':total_F_arr   # New column
                    })
                    
                    all_data.append(df_file)


if all_data:
    final_df = pd.concat(all_data, ignore_index=True)
    output_csv_path = "./fig7_defects_data.csv"
    final_df.to_csv(output_csv_path, index=False)
    print(f"Successfully saved all tracking data to {output_csv_path}")
else:
    print("No valid tracking data found across the specified conditions.")



