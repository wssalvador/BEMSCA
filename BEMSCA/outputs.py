import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
#    GLOBAL CONSTANTS
# -----------------------------------------------------------------------------
DEFAULT_BAR_WIDTH = 0.5
DPI = 600
THINING_FACTOR = 0.85
TOTAL_BAR_WIDTH = 0.15
Y_SHIFT_FACTOR = 1.1


# Function for organizing BEMSCA's outputs based on the obtained simulation results
def simulate_output(db_data, bio_params_name, simulation_results):
    # <>---------------- Terminal Output ----------------<>
    # Print headers
    print(f'\n||--------- PROPOSED WORKFLOW [{bio_params_name}] ---------||\n')
    print('<>--------- Planar Expansion Workflow ---------<>')

    # Organize data to make planar worflow easier to interpret for user
    platform = simulation_results.planar_expansion.planar_platform_data.name
    print(f'2D PLATFORM: {platform}\n')
    culture_volume = simulation_results.planar_expansion.planar_platform_data["Culture Volume"].item()
    planar_workflow_table = pd.DataFrame({'Surfaces': simulation_results.planar_expansion.planar_workflow})
    planar_workflow_table["Volume (L)"] = planar_workflow_table * culture_volume
    planar_workflow_table.index = [f'P{num}' for num in range(len(simulation_results.planar_expansion.planar_workflow))]
    print(planar_workflow_table)

    # Print important info
    print(f'\nINITIAL CELL NUMBER: {db_data["Bioprocess Parameters"].loc[bio_params_name, "Initial Cell Number"]:.2e}')
    print(f'TARGET CELL NUMBER: {simulation_results.planar_expansion.inoc_cells:.2e}')    
    print(f'DURATION: {simulation_results.planar_expansion.duration} days')
    print(f'COST: {simulation_results.planar_expansion.overall_cost:,.2f} €')
    
    # Print header
    print('\n<>------- Bioreactor Expansion Workflow -------<>')

    # Organize data to make bioreactor worflow easier to interpret for user
    print(f'BIOREACTORS: {simulation_results.bioreactor_expansion.bioreactors.index.tolist()}\n')
    bioreactor_workflow_table = simulation_results.bioreactor_expansion.bioreactor_workflow[1]
    bioreactor_workflow_table["Volume (L)"] = np.round(simulation_results.bioreactor_expansion.bioreactor_workflow[0], 3)
    bioreactor_workflow_table["Inoculated Cells"] = (bioreactor_workflow_table["Volume (L)"] * 1e3
                                                     * simulation_results.planar_expansion.seeding_density)
    bioreactor_workflow_table = bioreactor_workflow_table.astype({'Inoculated Cells': 'int64'})
    bioreactor_workflow_table.index = [f'C{num+1}' for num in range(len(bioreactor_workflow_table))]
    print(bioreactor_workflow_table)

    # Calculate the target cell number of bioreactor expansion (and overall bioprocess)
    fin_target_cell_number = db_data["Bioprocess Parameters"].loc[bio_params_name, "Target Cell Number"]

    # Print important info
    print(f'\nINITIAL CELL NUMBER: {simulation_results.planar_expansion.inoc_cells:.2e}')
    print(f'TARGET CELL NUMBER: {fin_target_cell_number:.2e}')
    print(f'''DURATION: {simulation_results.bioreactor_expansion.duration
                         - simulation_results.bioreactor_expansion.FIN_QUAL_DURATION} days''')
    print(f'FINAL QC DURATION: {simulation_results.bioreactor_expansion.FIN_QUAL_DURATION} days')
    print(f'COST: {simulation_results.bioreactor_expansion.overall_cost:,.2f} €')

    # Print header
    print('\n<>------------ Bioprocess Summary -----------<>\n')

    # Calculate relative category costs
    relative_consumables_cost = simulation_results.bioprocess_consumables_cost / simulation_results.bioprocess_overall_cost
    realtive_reagents_cost = simulation_results.bioprocess_reagents_cost / simulation_results.bioprocess_overall_cost
    relative_facility_cost = simulation_results.bioprocess_facility_cost / simulation_results.bioprocess_overall_cost
    relative_labor_cost = simulation_results.bioprocess_labor_cost / simulation_results.bioprocess_overall_cost

    # Round relative costs while ensuring they add up to 1
    relative_costs = [relative_consumables_cost, realtive_reagents_cost, relative_facility_cost, relative_labor_cost]
    decimal_parts = []
    for i, relative_cost in enumerate(relative_costs):
        decimal_parts.append(math.modf(relative_cost * 100)[0])
        relative_costs[i] = math.floor(relative_cost * 100) / 100

    while sum(relative_costs) < 1:
        maximum = max(decimal_parts)
        maximum_index = decimal_parts.index(maximum)
        relative_costs[maximum_index] += 0.01
        decimal_parts[maximum_index] = 0

    # Organize data to make bioprocess cost categories easier to interpret for user
    cost_categories = pd.DataFrame({'Consumables': [simulation_results.bioprocess_consumables_cost, relative_costs[0]],
                                   'Reagents': [simulation_results.bioprocess_reagents_cost, relative_costs[1]],
                                   'Facility': [simulation_results.bioprocess_facility_cost, relative_costs[2]],
                                   'Labor': [simulation_results.bioprocess_labor_cost, relative_costs[3]]},
                                   index=['Absolute Cost (€)', 'Relative Cost'])
    
    # Calculate relative medium costs (in %)    
    relative_medium_cost = simulation_results.bioprocess_medium_cost / simulation_results.bioprocess_overall_cost * 100    
    
    # Print cost categories table
    pd.options.display.float_format = "{:,.2f}".format
    print(cost_categories)

    # Calculate the average final cell number
    fin_cell_number_pd = (simulation_results.bioreactor_expansion.fold_increase_pds[-1]
                          * simulation_results.planar_expansion.inoc_cells)
    avg_fin_cell_number = np.average(fin_cell_number_pd)
    
    # Calculate the confidence level
    confidence_level = ((fin_cell_number_pd >= fin_target_cell_number).sum()
                        / simulation_results.bioreactor_expansion.SIMULATION_RUNS) * 100
    
    # Print important info
    print(f'\nAVERAGE FINAL CELL NUMBER: {avg_fin_cell_number:.2e}')
    print(f'CONFIDENCE LEVEL: {confidence_level:.1f}%')
    print(f'''OVERALL DURATION: {simulation_results.planar_expansion.duration
                                 + simulation_results.bioreactor_expansion.duration} days''')
    print(f'MEDIUM COST: {simulation_results.bioprocess_medium_cost:,.2f} € ({relative_medium_cost:.0f}%)')
    print(f'OVERALL COST: {simulation_results.bioprocess_overall_cost:,.2f} €\n')

    # <>---------------- Graphical Output ---------------<>
    # Create figure and axis
    figure, axes = plt.subplots(tight_layout=True)

    # Define graph title
    axes.set_title('Bioprocess costs by category', fontweight='bold')

    # Define position of bars
    y_pos = np.arange(len(cost_categories.columns))

    # Add bars to graph and define axes labels
    axes.barh(y_pos, cost_categories.loc["Absolute Cost (€)", :], DEFAULT_BAR_WIDTH, align='center', color='goldenrod')
    axes.set_yticks(y_pos)
    axes.set_yticklabels(cost_categories)
    axes.invert_yaxis()
    axes.set_xlabel('Cost (€)')
    
    # Add shading to reagents bar to highlight culture medium cost
    axes.barh(y_pos[1], simulation_results.bioprocess_medium_cost, DEFAULT_BAR_WIDTH,
              left=(simulation_results.bioprocess_reagents_cost - simulation_results.bioprocess_medium_cost),
              color='gold', label='Medium Cost')
    
    # Show graph legend
    axes.legend(loc='upper right')

    # Save figure to "results" folder
    plt.savefig(f'results/Cost_Categories_{bio_params_name}.png', dpi=DPI)

    return cost_categories, simulation_results.bioprocess_medium_cost


# Function for organizing BEMSCA's outputs when comparing two to five different conditions
def compare_output(db_data, bio_params_names, simulations, customization):
    # <>------------ User Graph Customization -----------<>            
    # Define base-case
    base_case_index = customization[0]
    
    # Define name
    name = customization[1]

    # Define colors
    primary_colors = db_data["Color Palettes"].loc[customization[2]]["Primary Colors"].split(";")

    secondary_colors = db_data["Color Palettes"].loc[customization[2]]["Secondary Colors"].split(";")
    
    # Define labels
    labels = customization[3]

    # Define medium labels
    medium_labels = []
    for label in labels:
        if label.split(" (")[0] not in medium_labels:
            medium_labels.append(label.split(" (")[0])
    
    # <>---------------- Graphical Output ---------------<>
    # Get organized outputs from each conditions simulation
    cost_categories = {}
    for bio_params_name, simulation in zip (bio_params_names, simulations):
        cost_categories[bio_params_name] = simulate_output(db_data, bio_params_name, simulation)

    # Create figure and axis (overwrite any existing figure)
    figure, axes = plt.subplots(tight_layout=True)

    # Define cost categories graph title
    axes.set_title(f'Impact of {name} on bioprocess category costs', fontweight='bold')
    
    # Adjust bar width according to number of conditions being compared
    bar_width = DEFAULT_BAR_WIDTH/2 * THINING_FACTOR**(len(simulations)-2)

    # Define position of bars
    y_pos = np.arange(len(cost_categories[bio_params_names[-1]][0].columns))

    # Determine appropriate initial y shift of bars according to number of conditions being compared
    if len(bio_params_names) % 2 == 0:
        y_shift = -bar_width * (Y_SHIFT_FACTOR/2) * (len(bio_params_names)/2)
    else:
        y_shift = -bar_width * Y_SHIFT_FACTOR * (len(bio_params_names)//2)

    # Construct cost categories graph
    for i, bio_params_name in enumerate(bio_params_names):
        # Add bars to graph
        axes.barh(y_pos + y_shift, cost_categories[bio_params_name][0].loc[['Absolute Cost (€)']].values[0],
                    bar_width, align='center', color=primary_colors[i], label=labels[i])
        
        # Add shading to reagents bar to highlight culture medium cost
        axes.barh(y_pos[1] + y_shift, cost_categories[bio_params_name][1], bar_width,
            left=(cost_categories[bio_params_name][0].loc['Absolute Cost (€)', 'Reagents']
                  - cost_categories[bio_params_name][1]), color=secondary_colors[i])
        
        # Update y shift
        y_shift += (bar_width * Y_SHIFT_FACTOR)

    # Define axes label
    axes.set_yticks(y_pos)
    axes.set_yticklabels(cost_categories[bio_params_names[-1]][0])
    axes.invert_yaxis()
    axes.set_xlabel('Cost (€)')

    # Show graph legend
    axes.legend(loc='upper right')

    # Define file name
    file_name = 'results/Cost_Categories_Comp'
    for bio_params_name in bio_params_names[1:]:
        file_name += f'_{bio_params_name}'
    file_name += '.png'

    # Save figure to "results" folder
    plt.savefig(file_name, dpi=DPI)

    # Clear axes in order to create new graph (with total costs)
    axes.clear()

    # Define total costs graph title
    axes.set_title(f'Impact of {name} on total bioprocess cost', fontweight='bold')

    # Adjust thick bar width according to number of conditions being compared
    bar_width = TOTAL_BAR_WIDTH * (THINING_FACTOR+0.05)**(len(simulations)-2)

    # Define y position of bars (for an aesthetic graph)
    y_pos = np.linspace(0.1+0.2**(len(simulations)-1), 0.9-0.2**(len(simulations)-1), len(simulations))

    # Create array with total costs
    total_costs = np.zeros(len(y_pos))
    for i, bio_params_name in enumerate(bio_params_names):
        total_costs[i] = cost_categories[bio_params_name][0].loc[['Absolute Cost (€)']].sum(axis=1).values[0]

    # Add bars to graph and define axes labels
    axes.set_ylim(0, 1)
    axes.barh(y_pos, total_costs, bar_width, align='center', color=primary_colors)
    axes.set_yticks(y_pos)
    axes.set_yticklabels(labels)
    axes.invert_yaxis()
    axes.set_xlabel('Total Cost (€)')

    # Add additional elements to graph
    for i, bio_params_name in enumerate(bio_params_names):
        # Add shading to reagents bar to highlight culture medium cost
        if i < len(medium_labels):
            bar = axes.barh(y_pos[i], cost_categories[bio_params_name][1], bar_width,
                            left=(total_costs[i] - cost_categories[bio_params_name][1]),
                            color=secondary_colors[i], label=f'{medium_labels[i]} medium cost')
        else:
            bar = axes.barh(y_pos[i], cost_categories[bio_params_name][1], bar_width,
                            left=(total_costs[i] - cost_categories[bio_params_name][1]),
                            color=secondary_colors[i])            
    
        # Calculate medium percentage and add to corresponding bar
        medium_percentage = cost_categories[bio_params_name][1] / total_costs[i] * 100
        axes.bar_label(bar, labels=['%.0f%%' % medium_percentage], label_type='center', fontweight='bold')

        # Calculate reduction percentage and add to corresponding bar
        if i != base_case_index:
            reduction_percentage = (1 - (total_costs[i] / total_costs[base_case_index])) * 100
            axes.bar_label(bar, labels=['(%.0f%%)' % reduction_percentage], padding=12)

    # Show graph legend
    axes.legend()

    # Define file name
    file_name = 'results/Total_Cost_Comp'
    for bio_params_name in bio_params_names[1:]:
        file_name += f'_{bio_params_name}'
    file_name += '.png'

    # Save figure to "results" folder
    plt.savefig(file_name, dpi=DPI)