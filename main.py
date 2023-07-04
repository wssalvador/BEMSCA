import numpy as np
import matplotlib.pyplot as plt
from database import *
from phases import *
from logistics import *
from utils import *


# Define Command Functions
def conditions_command(cc):
    print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
    print(f'Target Cell Number: {cc.target_cell_number:.2e}')
    print(f'Culture Medium: {cc.culture_medium}')
    print(f'Initial Inoculation Volume: {cc.initial_volume} mL')
    print(f'Simulation Source: {cc.simulation_source} et al.')
    print(f'Seeding Density: {cc.seeding_density:.0f} cells/mL')
    print(f'DS Supplementation: {cc.DS_supplementation}')
    print(f'--> See files in "simul" folder for respective distributions.')


def custom_command(cc):
    value_icn = input('Initial Cell Number: ')
    value_icn = int(float(value_icn))

    value_tcn = input('Target Cell Number: ')
    value_tcn = int(float(value_tcn))

    value_cm = input('Culture Medium [mTeSR1, E8, ihE8 (in-house E8) or B8]: ')
    value_cm = str(value_cm)

    value_iv = input('Initial Inoculation Volume (mL): ')
    value_iv = int(float(value_iv))

    custom_cc = CultureConditions(initnum=value_icn, targnum=value_tcn, medium=value_cm, initvol=value_iv,
                                  sim=cc.simulation_source, DS=cc.DS_supplementation)
    return custom_cc


def data_command(dt):
    print(f'---REAGENTS---\n{dt.reagents_costs}\n')
    print(f'---CONSUMABLES---')
    print(f'--[PEPS]--\n{dt.peps_costs}\n{dt.peps_parameters}\n')
    print(f'{dt.peps_cells.to_string(formatters=["{:.2e}".format, "{:.2e}".format])}\n')
    print(f'<6-well>\n{dt.well_parameters.to_string(index=False)}\n----------')
    print(f'--[VWBS]--\n{dt.vwbs_parameters}\n----------\n')
    print(f'---CULTURE MEDIA---\n{dt.culture_media_costs}')


def default_command():
    reset_cc = CultureConditions()
    return reset_cc


def help_command():
    print('"conditions" - shows the culture conditions currently set for the bioprocess')
    print('               model.')
    print('    "custom" - allows choosing custom culture conditions for the bioprocess')
    print('               model.')
    print('      "data" - shows the data table entries used by the bioprocess model.')
    print('   "default" - resets the default culture conditions for the bioprocess model.')
    print('      "exit" - terminates the program.')
    print('      "help" - shows a list of possible commands.')
    print('       "run" - runs the bioprocess model with the current conditions and')
    print('               simulation source.')
    print('    "source" - allows changing the simulation source for the bioprocess model.')
    print('     "study" - runs one of the preset comparisons of different experimental')
    print('               setups.')


def run_command(dt, ft, cc, fold_expansion_distribution, recovery_efficiency_distribution, init_sim_num):

    bioprocess = Bioprocess(dt, ft, cc, fold_expansion_distribution, recovery_efficiency_distribution,
                            init_sim_num)

    print('---PROPOSED WORKFLOW---')
    print('--[2D Expansion]--')
    print(f'{bioprocess.passaging_workflow}\n')
    print(f'DURATION: {bioprocess.cycle_info.at["2D Exp", "Duration (d)"]} days')
    print(f'EXECUTION TIME: {bioprocess.cycle_info.at["2D Exp", "Worktime (h)"]} hours')
    print(f'COST: {round(bioprocess.cycle_info.at["2D Exp", "Phase Cost (€)"], 2)} €')
    print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
    print(f'Average Final Cell Number: {bioprocess.cycle_info.at["2D Exp", "Cell Number"]:.2e}\n')

    print('--[Bioreactor Expansion]--')
    for expansion_cycle, current_workflow in enumerate(bioprocess.bioreactor_workflow):
        print(f'<Cycle C{expansion_cycle}>')
        print(f'{current_workflow}\n')
        initial_cell_number = cc.seeding_density * bioprocess.optimal_working_volumes[expansion_cycle]
        print(f'DURATION: {bioprocess.cycle_info.at[f"C{expansion_cycle}", "Duration (d)"]} days')
        print(f'EXECUTION TIME: {bioprocess.cycle_info.at[f"C{expansion_cycle}", "Worktime (h)"]} hours')
        print(f'COST: {round(bioprocess.cycle_info.at[f"C{expansion_cycle}", "Phase Cost (€)"], 2)} €')
        print(f'Initial Cell Number: {initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {bioprocess.cycle_info.at[f"C{expansion_cycle}", "Cell Number"]:.2e}\n')

    print('---BIOPROCESS RESULT---')
    print(f'CONFIDENCE LEVEL: {round(bioprocess.success_probability * 100, 1)} %')
    print(f'TOTAL DURATION: {bioprocess.total_duration} days')
    print(f'TOTAL EXECUTION TIME: {bioprocess.total_worktime} hours')
    print(f'TOTAL COST: {bioprocess.total_cost} €')
    print(f'Total Direct Cost: {bioprocess.total_direct_cost} €')
    print(f'Total Indirect Cost: {bioprocess.total_indirect_cost} €')
    print(f'Total Failure Cost: {bioprocess.total_failure_cost} €')
    print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
    print(f'Average Final Cell Number: {bioprocess.final_cell_number:.2e}\n')

    print('---BIOPROCESS INFO---')
    print(f'''{bioprocess.cycle_info.to_string(formatters=["{:.0f}".format, "{:.1f}".format, "{:.2e}".format,
          "{:.2f}".format, "{:.2f}".format])}\n''')
    print(f'{bioprocess.cycle_cost_categories}\n')
    print(f'{bioprocess.cycle_cost_stages}\n')
    print('---COST CATEGORIES---')
    print(f'Consumables: {round(bioprocess.total_consumables_cost, 2)} €')
    print(f'Reagents: {round(bioprocess.total_reagents_cost, 2)} €')
    print(f'Facility: {round(bioprocess.total_facility_cost, 2)} €')
    print(f'Labor: {round(bioprocess.total_labor_cost, 2)} €\n')

    # Bar Charts
    # create bar charts with cost distributions and save them in the respective folder
    bar_width = 0.5
    dots_per_inch = 600

    # Cost per Category
    # absolute cost distribution
    figure, axis = plt.subplots(tight_layout=True)

    axis.set_title('Bioprocess costs by category', fontweight='bold')

    cost_categories = ('Consumables', 'Reagents', 'Facility', 'Labor')
    y_pos = np.arange(len(cost_categories))
    costs = [bioprocess.total_consumables_cost, bioprocess.total_reagents_cost,
             bioprocess.total_facility_cost, bioprocess.total_labor_cost]

    axis.barh(y_pos, costs, bar_width, align='center', color='goldenrod')
    axis.set_yticks(y_pos)
    axis.set_yticklabels(cost_categories)
    axis.invert_yaxis()
    axis.set_xlabel('Cost (€)')

    plt.savefig('results/absolute_cost_category.png', dpi=dots_per_inch)

    # Cost per Stage
    # absolute cost distribution
    axis.clear()

    axis.set_title('Bioprocess costs by stage', fontweight='bold')

    cost_categories = ('Pre-Inoculation', 'Expansion', 'Quality Control', 'Harvesting')
    y_pos = np.arange(len(cost_categories))
    costs = [bioprocess.total_pre_inoculation_cost, bioprocess.total_expansion_cost,
             bioprocess.total_quality_control_cost, bioprocess.total_harvesting_cost]
    axis.barh(y_pos, costs, bar_width, align='center', color='goldenrod')
    axis.set_yticks(y_pos)
    axis.set_yticklabels(cost_categories)
    axis.invert_yaxis()
    axis.set_xlabel('Cost (€)')

    plt.savefig('results/absolute_cost_stage.png', dpi=dots_per_inch)

    # Total Costs
    bar_width = 0.6
    axis.clear()

    axis.set_title('Total bioprocess cost', fontweight='bold')

    labels = ['mTeSR1']
    colors = ['red', 'red', 'goldenrod', 'red', 'red']
    y_pos = [1, 2, 3, 4, 5]
    costs = [0, 0, bioprocess.total_cost - bioprocess.total_medium_cost, 0, 0]

    axis.barh(y_pos, costs, bar_width, align='center', color=colors)
    axis.set_yticks([3])
    axis.set_yticklabels(labels)
    axis.invert_yaxis()
    axis.set_xlabel('Total Cost (€)')

    # add medium cost seperately from other reagent costs
    medium_percentages = [bioprocess.total_medium_cost / bioprocess.total_cost * 100]

    bar = axis.barh(y_pos[2], bioprocess.total_medium_cost, bar_width, left=costs[2],
                           color='gold', label='Medium Cost')
    axis.bar_label(bar, labels=['%.0f%%' % medium_percentages[0]], label_type='center',
                   fontweight='bold')
    axis.legend()

    plt.savefig('results/total_cost.png', dpi=dots_per_inch)

    # final cell number distribution
    figure, axis = plt.subplots()
    bin_number = 200
    simulation_number = int(1e5)
    axis.clear()
    figure.suptitle('Final Cell Number Distribution')
    axis.hist(bioprocess.final_cycle_cell_number_distribution, bin_number, color='firebrick')
    axis.set_xlabel('Cell Number')
    axis.set_ylabel('Frequency')
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=simulation_number))
    plt.savefig('simul/final_cell_number.png', dpi=dots_per_inch)


def source_command(cc):
    value_ss = input('Simulation Source [Nogueira, Borys]: ')
    value_ss = str(value_ss)
    if value_ss == 'Nogueira':
        value_DS = input('DS Supplementation [True, False]: ')
        if value_DS == 'True':
            value_DS = True
        else:
            value_DS = False
    else:
        value_DS = False

    custom_cc = CultureConditions(initnum=cc.initial_cell_number, targnum=cc.target_cell_number,
                                  medium=cc.culture_medium, initvol=cc.initial_volume, sim=value_ss, DS=value_DS)
    return custom_cc


def study_command(dt, ft, cc):
    value_st = input('Select study to run [DS, density, B8, qual, general]: ')
    if value_st == 'DS' or value_st == 'general':
        # simulate bioprocesses relevant to the comparison of medium without and with DS
        # without DS
        without_DS_cc = CultureConditions(initnum=cc.initial_cell_number, targnum=cc.target_cell_number,
                                  medium='mTeSR1', initvol=cc.initial_volume, sim='Nogueira', DS=False)

        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(dt, without_DS_cc)

        mTeSR1_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                recovery_efficiency_distribution, init_sim_num)

        print('---PROPOSED WORKFLOW (-DS)---')
        print('--[2D Expansion]--')
        print(f'{mTeSR1_bioprocess.passaging_workflow}\n')
        print(f'DURATION: {mTeSR1_bioprocess.cycle_info.at["2D Exp", "Duration (d)"]} days')
        print(f'EXECUTION TIME: {mTeSR1_bioprocess.cycle_info.at["2D Exp", "Worktime (h)"]} hours')
        print(f'COST: {round(mTeSR1_bioprocess.cycle_info.at["2D Exp", "Phase Cost (€)"], 2)} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {mTeSR1_bioprocess.cycle_info.at["2D Exp", "Cell Number"]:.2e}\n')

        print('--[Bioreactor Expansion]--')
        for expansion_cycle, current_workflow in enumerate(mTeSR1_bioprocess.bioreactor_workflow):
            print(f'<Cycle C{expansion_cycle}>')
            print(f'{current_workflow}\n')
            initial_cell_number = cc.seeding_density * mTeSR1_bioprocess.optimal_working_volumes[expansion_cycle]
            print(f'DURATION: {mTeSR1_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Duration (d)"]} days')
            print(f'EXECUTION TIME: {mTeSR1_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Worktime (h)"]} hours')
            print(f'COST: {round(mTeSR1_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Phase Cost (€)"], 2)} €')
            print(f'Initial Cell Number: {initial_cell_number:.2e}')
            print(f'Average Final Cell Number: {mTeSR1_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Cell Number"]:.2e}\n')

        print('---BIOPROCESS RESULT---')
        print(f'CONFIDENCE LEVEL: {round(mTeSR1_bioprocess.success_probability * 100, 1)} %')
        print(f'TOTAL DURATION: {mTeSR1_bioprocess.total_duration} days')
        print(f'TOTAL EXECUTION TIME: {mTeSR1_bioprocess.total_worktime} hours')
        print(f'TOTAL COST: {mTeSR1_bioprocess.total_cost} €')
        print(f'Total Direct Cost: {mTeSR1_bioprocess.total_direct_cost} €')
        print(f'Total Indirect Cost: {mTeSR1_bioprocess.total_indirect_cost} €')
        print(f'Total Failure Cost: {mTeSR1_bioprocess.total_failure_cost} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {mTeSR1_bioprocess.final_cell_number:.2e}\n')

        print('---BIOPROCESS INFO---')
        print(f'''{mTeSR1_bioprocess.cycle_info.to_string(formatters=["{:.0f}".format, "{:.1f}".format, "{:.2e}".format,
              "{:.2f}".format, "{:.2f}".format])}\n''')
        print(f'{mTeSR1_bioprocess.cycle_cost_categories}\n')
        print(f'{mTeSR1_bioprocess.cycle_cost_stages}\n')
        print('---COST CATEGORIES---')
        print(f'Consumables: {round(mTeSR1_bioprocess.total_consumables_cost, 2)} €')
        print(f'Reagents: {round(mTeSR1_bioprocess.total_reagents_cost, 2)} €')
        print(f'Facility: {round(mTeSR1_bioprocess.total_facility_cost, 2)} €')
        print(f'Labor: {round(mTeSR1_bioprocess.total_labor_cost, 2)} €\n')

        # with DS
        with_DS_cc = CultureConditions(initnum=cc.initial_cell_number, targnum=cc.target_cell_number,
                                  medium='mTeSR1', initvol=cc.initial_volume, sim='Nogueira', DS=True)

        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(dt, with_DS_cc)

        mTeSR1_DS_bioprocess = Bioprocess(dt, ft, with_DS_cc, fold_expansion_distribution,
                                recovery_efficiency_distribution, init_sim_num)

        print('---PROPOSED WORKFLOW (+DS)---')
        print('--[2D Expansion]--')
        print(f'{mTeSR1_DS_bioprocess.passaging_workflow}\n')
        print(f'DURATION: {mTeSR1_DS_bioprocess.cycle_info.at["2D Exp", "Duration (d)"]} days')
        print(f'EXECUTION TIME: {mTeSR1_DS_bioprocess.cycle_info.at["2D Exp", "Worktime (h)"]} hours')
        print(f'COST: {round(mTeSR1_DS_bioprocess.cycle_info.at["2D Exp", "Phase Cost (€)"], 2)} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {mTeSR1_DS_bioprocess.cycle_info.at["2D Exp", "Cell Number"]:.2e}\n')

        print('--[Bioreactor Expansion]--')
        for expansion_cycle, current_workflow in enumerate(mTeSR1_DS_bioprocess.bioreactor_workflow):
            print(f'<Cycle C{expansion_cycle}>')
            print(f'{current_workflow}\n')
            initial_cell_number = cc.seeding_density * mTeSR1_DS_bioprocess.optimal_working_volumes[expansion_cycle]
            print(f'DURATION: {mTeSR1_DS_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Duration (d)"]} days')
            print(f'EXECUTION TIME: {mTeSR1_DS_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Worktime (h)"]} hours')
            print(f'COST: {round(mTeSR1_DS_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Phase Cost (€)"], 2)} €')
            print(f'Initial Cell Number: {initial_cell_number:.2e}')
            print(f'Average Final Cell Number: {mTeSR1_DS_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Cell Number"]:.2e}\n')

        print('---BIOPROCESS RESULT---')
        print(f'CONFIDENCE LEVEL: {round(mTeSR1_DS_bioprocess.success_probability * 100, 1)} %')
        print(f'TOTAL DURATION: {mTeSR1_DS_bioprocess.total_duration} days')
        print(f'TOTAL EXECUTION TIME: {mTeSR1_DS_bioprocess.total_worktime} hours')
        print(f'TOTAL COST: {mTeSR1_DS_bioprocess.total_cost} €')
        print(f'Total Direct Cost: {mTeSR1_DS_bioprocess.total_direct_cost} €')
        print(f'Total Indirect Cost: {mTeSR1_DS_bioprocess.total_indirect_cost} €')
        print(f'Total Failure Cost: {mTeSR1_DS_bioprocess.total_failure_cost} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {mTeSR1_DS_bioprocess.final_cell_number:.2e}\n')

        print('---BIOPROCESS INFO---')
        print(f'''{mTeSR1_DS_bioprocess.cycle_info.to_string(formatters=["{:.0f}".format, "{:.1f}".format, "{:.2e}".format,
              "{:.2f}".format, "{:.2f}".format])}\n''')
        print(f'{mTeSR1_DS_bioprocess.cycle_cost_categories}\n')
        print(f'{mTeSR1_DS_bioprocess.cycle_cost_stages}\n')
        print('---COST CATEGORIES---')
        print(f'Consumables: {round(mTeSR1_DS_bioprocess.total_consumables_cost, 2)} €')
        print(f'Reagents: {round(mTeSR1_DS_bioprocess.total_reagents_cost, 2)} €')
        print(f'Facility: {round(mTeSR1_DS_bioprocess.total_facility_cost, 2)} €')
        print(f'Labor: {round(mTeSR1_DS_bioprocess.total_labor_cost, 2)} €\n')

        # Bar Charts
        # create bar charts with cost distributions and save them in the respective folder
        bar_width = 0.25
        dots_per_inch = 600

        # Cost per Category
        # absolute cost distribution
        figure, axis = plt.subplots(tight_layout=True)

        axis.set_title('Impact of DS on bioprocess costs by category', fontweight='bold')

        cost_categories = ('Consumables', 'Reagents', 'Facility', 'Labor')
        y_pos = np.arange(len(cost_categories))
        mTeSR1_costs = [mTeSR1_bioprocess.total_consumables_cost,
                 mTeSR1_bioprocess.total_reagents_cost, mTeSR1_bioprocess.total_facility_cost, mTeSR1_bioprocess.total_labor_cost]
        mTeSR1_DS_costs = [mTeSR1_DS_bioprocess.total_consumables_cost,
                 mTeSR1_DS_bioprocess.total_reagents_cost, mTeSR1_DS_bioprocess.total_facility_cost, mTeSR1_DS_bioprocess.total_labor_cost]

        axis.barh(y_pos - bar_width * 0.55, mTeSR1_costs, bar_width, align='center', color='goldenrod',
                  label='mTeSR1')
        axis.barh(y_pos + bar_width * 0.55, mTeSR1_DS_costs, bar_width, align='center', color='firebrick',
                  label='mTeSR1+DS')
        axis.legend()
        axis.set_yticks(y_pos)
        axis.set_yticklabels(cost_categories)
        axis.invert_yaxis()
        axis.set_xlabel('Cost (€)')

        plt.savefig('studies/DS_comparison_absolute_cost_category.png', dpi=dots_per_inch)

        # Cost per Stage
        # absolute cost distribution
        axis.clear()

        axis.set_title('Impact of DS on bioprocess costs by stage', fontweight='bold')

        cost_categories = ('Pre-Inoculation', 'Expansion', 'Quality Control', 'Harvesting')
        y_pos = np.arange(len(cost_categories))
        mTeSR1_costs = [mTeSR1_bioprocess.total_pre_inoculation_cost, mTeSR1_bioprocess.total_expansion_cost,
                        mTeSR1_bioprocess.total_quality_control_cost, mTeSR1_bioprocess.total_harvesting_cost]
        mTeSR1_DS_costs = [mTeSR1_DS_bioprocess.total_pre_inoculation_cost,
                        mTeSR1_DS_bioprocess.total_expansion_cost, mTeSR1_DS_bioprocess.total_quality_control_cost, mTeSR1_DS_bioprocess.total_harvesting_cost]
        axis.barh(y_pos - bar_width * 0.55, mTeSR1_costs, bar_width, align='center', color='goldenrod',
                  label='mTeSR1')
        axis.barh(y_pos + bar_width * 0.55, mTeSR1_DS_costs, bar_width, align='center', color='firebrick',
                  label='mTeSR1+DS')
        axis.legend()
        axis.set_yticks(y_pos)
        axis.set_yticklabels(cost_categories)
        axis.invert_yaxis()
        axis.set_xlabel('Cost (€)')

        plt.savefig('studies/DS_comparison_absolute_cost_stage.png', dpi=dots_per_inch)

        # Total Costs
        bar_width = 0.6
        axis.clear()

        axis.set_title('Impact of DS on total bioprocess cost', fontweight='bold')

        labels = ['mTeSR1', 'mTeSR1+DS']
        colors = ['red', 'goldenrod', 'red', 'firebrick', 'red']
        y_pos = [1, 2, 3, 4, 5]
        costs = [0, mTeSR1_bioprocess.total_cost - mTeSR1_bioprocess.total_medium_cost, 0,
                 mTeSR1_DS_bioprocess.total_cost - mTeSR1_DS_bioprocess.total_medium_cost, 0]

        axis.barh(y_pos, costs, bar_width, align='center', color=colors)
        axis.set_yticks([2, 4])
        axis.set_yticklabels(labels)
        axis.invert_yaxis()
        axis.set_xlabel('Total Cost (€)')

        # add medium cost seperately from other reagent costs
        medium_percentages = [mTeSR1_bioprocess.total_medium_cost / mTeSR1_bioprocess.total_cost * 100,
                              mTeSR1_DS_bioprocess.total_medium_cost / mTeSR1_DS_bioprocess.total_cost * 100]
        reduction_percentages = [(1 - (mTeSR1_DS_bioprocess.total_cost / mTeSR1_bioprocess.total_cost)) * 100]

        mTeSR1_bar = axis.barh(y_pos[1], mTeSR1_bioprocess.total_medium_cost, bar_width, left=costs[1],
                               color='gold', label='mTeSR1 medium cost')
        mTeSR1_DS_bar = axis.barh(y_pos[3], mTeSR1_DS_bioprocess.total_medium_cost, bar_width,
                                  left=costs[3], color='salmon', label='mTeSR1+DS medium cost')
        axis.bar_label(mTeSR1_DS_bar, labels=['(%.0f%%)' % reduction_percentages[0]], padding=12)
        axis.bar_label(mTeSR1_bar, labels=['%.0f%%' % medium_percentages[0]], label_type='center',
                       fontweight='bold')
        axis.bar_label(mTeSR1_DS_bar, labels=['%.0f%%' % medium_percentages[1]], label_type='center',
                       fontweight='bold')
        axis.legend()

        plt.savefig('studies/DS_comparison_total_cost.png', dpi=dots_per_inch)

    if value_st == 'B8' or value_st == 'general':
        # simulate bioprocesses relevant to the comparison of B8 with mTeSR1+DS considering different efficacies
        # mTeSR1+DS
        without_DS_cc = CultureConditions(initnum=cc.initial_cell_number, targnum=cc.target_cell_number,
                                  medium='mTeSR1', initvol=cc.initial_volume, sim='Nogueira', DS=False)

        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(dt, without_DS_cc)

        mTeSR1_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                recovery_efficiency_distribution, init_sim_num)

        print('---PROPOSED WORKFLOW (mTeSR1)---')
        print('--[2D Expansion]--')
        print(f'{mTeSR1_bioprocess.passaging_workflow}\n')
        print(f'DURATION: {mTeSR1_bioprocess.cycle_info.at["2D Exp", "Duration (d)"]} days')
        print(f'EXECUTION TIME: {mTeSR1_bioprocess.cycle_info.at["2D Exp", "Worktime (h)"]} hours')
        print(f'COST: {round(mTeSR1_bioprocess.cycle_info.at["2D Exp", "Phase Cost (€)"], 2)} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {mTeSR1_bioprocess.cycle_info.at["2D Exp", "Cell Number"]:.2e}\n')

        print('--[Bioreactor Expansion]--')
        for expansion_cycle, current_workflow in enumerate(mTeSR1_bioprocess.bioreactor_workflow):
            print(f'<Cycle C{expansion_cycle}>')
            print(f'{current_workflow}\n')
            initial_cell_number = cc.seeding_density * mTeSR1_bioprocess.optimal_working_volumes[expansion_cycle]
            print(f'DURATION: {mTeSR1_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Duration (d)"]} days')
            print(f'EXECUTION TIME: {mTeSR1_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Worktime (h)"]} hours')
            print(f'COST: {round(mTeSR1_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Phase Cost (€)"], 2)} €')
            print(f'Initial Cell Number: {initial_cell_number:.2e}')
            print(f'Average Final Cell Number: {mTeSR1_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Cell Number"]:.2e}\n')

        print('---BIOPROCESS RESULT---')
        print(f'CONFIDENCE LEVEL: {round(mTeSR1_bioprocess.success_probability * 100, 1)} %')
        print(f'TOTAL DURATION: {mTeSR1_bioprocess.total_duration} days')
        print(f'TOTAL EXECUTION TIME: {mTeSR1_bioprocess.total_worktime} hours')
        print(f'TOTAL COST: {mTeSR1_bioprocess.total_cost} €')
        print(f'Total Direct Cost: {mTeSR1_bioprocess.total_direct_cost} €')
        print(f'Total Indirect Cost: {mTeSR1_bioprocess.total_indirect_cost} €')
        print(f'Total Failure Cost: {mTeSR1_bioprocess.total_failure_cost} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {mTeSR1_bioprocess.final_cell_number:.2e}\n')

        print('---BIOPROCESS INFO---')
        print(f'''{mTeSR1_bioprocess.cycle_info.to_string(formatters=["{:.0f}".format, "{:.1f}".format, "{:.2e}".format,
              "{:.2f}".format, "{:.2f}".format])}\n''')
        print(f'{mTeSR1_bioprocess.cycle_cost_categories}\n')
        print(f'{mTeSR1_bioprocess.cycle_cost_stages}\n')
        print('---COST CATEGORIES---')
        print(f'Consumables: {round(mTeSR1_bioprocess.total_consumables_cost, 2)} €')
        print(f'Reagents: {round(mTeSR1_bioprocess.total_reagents_cost, 2)} €')
        print(f'Facility: {round(mTeSR1_bioprocess.total_facility_cost, 2)} €')
        print(f'Labor: {round(mTeSR1_bioprocess.total_labor_cost, 2)} €\n')

        B8_cc = CultureConditions(initnum=cc.initial_cell_number, targnum=cc.target_cell_number,
                                  medium='B8', initvol=cc.initial_volume, sim='Nogueira', DS=False)

        # B8 0.75x efficacy (and 8 days of expansion instead of 7)
        B8_cc.expansion_duration = 8
        B8_cc.fold_increase_average = without_DS_cc.fold_increase_average * 0.75
        B8_cc.fold_increase_std = without_DS_cc.fold_increase_std * 0.75

        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(dt, B8_cc)

        for index, value in enumerate(fold_expansion_distribution):
            fold_expansion_distribution[index] *= 0.75

        print(f"BEFORE: {ft.overall_daily_facility_costs}")

        ft.overall_daily_facility_costs = round(ft.overall_daily_facility_costs + 6
                     * dt.vwbs_base_units.at['PBS 3MAG', 'Cost (€)'] / dt.equipment_depreciation_period, 2)

        print(f"AFTER: {ft.overall_daily_facility_costs}")

        B8_075x_bioprocess = Bioprocess(dt, ft, B8_cc, fold_expansion_distribution,
                                recovery_efficiency_distribution, init_sim_num)

        ft.overall_daily_facility_costs = round(ft.overall_daily_facility_costs - 6
                     * dt.vwbs_base_units.at['PBS 3MAG', 'Cost (€)'] / dt.equipment_depreciation_period, 2)

        print(f"BEFORE: {ft.overall_daily_facility_costs}")

        print('---PROPOSED WORKFLOW (B8 0.75x)---')
        print('--[2D Expansion]--')
        print(f'{B8_075x_bioprocess.passaging_workflow}\n')
        print(f'DURATION: {B8_075x_bioprocess.cycle_info.at["2D Exp", "Duration (d)"]} days')
        print(f'EXECUTION TIME: {B8_075x_bioprocess.cycle_info.at["2D Exp", "Worktime (h)"]} hours')
        print(f'COST: {round(B8_075x_bioprocess.cycle_info.at["2D Exp", "Phase Cost (€)"], 2)} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {B8_075x_bioprocess.cycle_info.at["2D Exp", "Cell Number"]:.2e}\n')

        print('--[Bioreactor Expansion]--')
        for expansion_cycle, current_workflow in enumerate(B8_075x_bioprocess.bioreactor_workflow):
            print(f'<Cycle C{expansion_cycle}>')
            print(f'{current_workflow}\n')
            initial_cell_number = cc.seeding_density * B8_075x_bioprocess.optimal_working_volumes[expansion_cycle]
            print(f'DURATION: {B8_075x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Duration (d)"]} days')
            print(f'EXECUTION TIME: {B8_075x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Worktime (h)"]} hours')
            print(f'COST: {round(B8_075x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Phase Cost (€)"], 2)} €')
            print(f'Initial Cell Number: {initial_cell_number:.2e}')
            print(f'Average Final Cell Number: {B8_075x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Cell Number"]:.2e}\n')

        print('---BIOPROCESS RESULT---')
        print(f'CONFIDENCE LEVEL: {round(B8_075x_bioprocess.success_probability * 100, 1)} %')
        print(f'TOTAL DURATION: {B8_075x_bioprocess.total_duration} days')
        print(f'TOTAL EXECUTION TIME: {B8_075x_bioprocess.total_worktime} hours')
        print(f'TOTAL COST: {B8_075x_bioprocess.total_cost} €')
        print(f'Total Direct Cost: {B8_075x_bioprocess.total_direct_cost} €')
        print(f'Total Indirect Cost: {B8_075x_bioprocess.total_indirect_cost} €')
        print(f'Total Failure Cost: {B8_075x_bioprocess.total_failure_cost} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {B8_075x_bioprocess.final_cell_number:.2e}\n')

        print('---BIOPROCESS INFO---')
        print(f'''{B8_075x_bioprocess.cycle_info.to_string(formatters=["{:.0f}".format, "{:.1f}".format, "{:.2e}".format,
              "{:.2f}".format, "{:.2f}".format])}\n''')
        print(f'{B8_075x_bioprocess.cycle_cost_categories}\n')
        print(f'{B8_075x_bioprocess.cycle_cost_stages}\n')
        print('---COST CATEGORIES---')
        print(f'Consumables: {round(B8_075x_bioprocess.total_consumables_cost, 2)} €')
        print(f'Reagents: {round(B8_075x_bioprocess.total_reagents_cost, 2)} €')
        print(f'Facility: {round(B8_075x_bioprocess.total_facility_cost, 2)} €')
        print(f'Labor: {round(B8_075x_bioprocess.total_labor_cost, 2)} €\n')

        # B8 1x efficacy (and standard 7 days of expansion)
        B8_cc.expansion_duration = 7
        B8_cc.fold_increase_average = without_DS_cc.fold_increase_average * 1
        B8_cc.fold_increase_std = without_DS_cc.fold_increase_std * 1

        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(dt, B8_cc)

        B8_1x_bioprocess = Bioprocess(dt, ft, B8_cc, fold_expansion_distribution,
                                recovery_efficiency_distribution, init_sim_num)

        print('---PROPOSED WORKFLOW (B8 1x)---')
        print('--[2D Expansion]--')
        print(f'{B8_1x_bioprocess.passaging_workflow}\n')
        print(f'DURATION: {B8_1x_bioprocess.cycle_info.at["2D Exp", "Duration (d)"]} days')
        print(f'EXECUTION TIME: {B8_1x_bioprocess.cycle_info.at["2D Exp", "Worktime (h)"]} hours')
        print(f'COST: {round(B8_1x_bioprocess.cycle_info.at["2D Exp", "Phase Cost (€)"], 2)} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {B8_1x_bioprocess.cycle_info.at["2D Exp", "Cell Number"]:.2e}\n')

        print('--[Bioreactor Expansion]--')
        for expansion_cycle, current_workflow in enumerate(B8_1x_bioprocess.bioreactor_workflow):
            print(f'<Cycle C{expansion_cycle}>')
            print(f'{current_workflow}\n')
            initial_cell_number = cc.seeding_density * B8_1x_bioprocess.optimal_working_volumes[expansion_cycle]
            print(f'DURATION: {B8_1x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Duration (d)"]} days')
            print(f'EXECUTION TIME: {B8_1x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Worktime (h)"]} hours')
            print(f'COST: {round(B8_1x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Phase Cost (€)"], 2)} €')
            print(f'Initial Cell Number: {initial_cell_number:.2e}')
            print(f'Average Final Cell Number: {B8_1x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Cell Number"]:.2e}\n')

        print('---BIOPROCESS RESULT---')
        print(f'CONFIDENCE LEVEL: {round(B8_1x_bioprocess.success_probability * 100, 1)} %')
        print(f'TOTAL DURATION: {B8_1x_bioprocess.total_duration} days')
        print(f'TOTAL EXECUTION TIME: {B8_1x_bioprocess.total_worktime} hours')
        print(f'TOTAL COST: {B8_1x_bioprocess.total_cost} €')
        print(f'Total Direct Cost: {B8_1x_bioprocess.total_direct_cost} €')
        print(f'Total Indirect Cost: {B8_1x_bioprocess.total_indirect_cost} €')
        print(f'Total Failure Cost: {B8_1x_bioprocess.total_failure_cost} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {B8_1x_bioprocess.final_cell_number:.2e}\n')

        print('---BIOPROCESS INFO---')
        print(f'''{B8_1x_bioprocess.cycle_info.to_string(formatters=["{:.0f}".format, "{:.1f}".format, "{:.2e}".format,
              "{:.2f}".format, "{:.2f}".format])}\n''')
        print(f'{B8_1x_bioprocess.cycle_cost_categories}\n')
        print(f'{B8_1x_bioprocess.cycle_cost_stages}\n')
        print('---COST CATEGORIES---')
        print(f'Consumables: {round(B8_1x_bioprocess.total_consumables_cost, 2)} €')
        print(f'Reagents: {round(B8_1x_bioprocess.total_reagents_cost, 2)} €')
        print(f'Facility: {round(B8_1x_bioprocess.total_facility_cost, 2)} €')
        print(f'Labor: {round(B8_1x_bioprocess.total_labor_cost, 2)} €\n')

        # B8 1.5x efficacy (and 6 days of expansion instead of 7)
        B8_cc.expansion_duration = 6
        B8_cc.fold_increase_average = without_DS_cc.fold_increase_average * 1.5
        B8_cc.fold_increase_std = without_DS_cc.fold_increase_std * 1.5

        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(dt, B8_cc)

        for index, value in enumerate(fold_expansion_distribution):
            fold_expansion_distribution[index] *= 1.5

        B8_15x_bioprocess = Bioprocess(dt, ft, B8_cc, fold_expansion_distribution,
                                recovery_efficiency_distribution, init_sim_num)

        print('---PROPOSED WORKFLOW (B8 1.5x)---')
        print('--[2D Expansion]--')
        print(f'{B8_15x_bioprocess.passaging_workflow}\n')
        print(f'DURATION: {B8_15x_bioprocess.cycle_info.at["2D Exp", "Duration (d)"]} days')
        print(f'EXECUTION TIME: {B8_15x_bioprocess.cycle_info.at["2D Exp", "Worktime (h)"]} hours')
        print(f'COST: {round(B8_15x_bioprocess.cycle_info.at["2D Exp", "Phase Cost (€)"], 2)} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {B8_15x_bioprocess.cycle_info.at["2D Exp", "Cell Number"]:.2e}\n')

        print('--[Bioreactor Expansion]--')
        for expansion_cycle, current_workflow in enumerate(B8_15x_bioprocess.bioreactor_workflow):
            print(f'<Cycle C{expansion_cycle}>')
            print(f'{current_workflow}\n')
            initial_cell_number = cc.seeding_density * B8_15x_bioprocess.optimal_working_volumes[expansion_cycle]
            print(f'DURATION: {B8_15x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Duration (d)"]} days')
            print(f'EXECUTION TIME: {B8_15x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Worktime (h)"]} hours')
            print(f'COST: {round(B8_15x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Phase Cost (€)"], 2)} €')
            print(f'Initial Cell Number: {initial_cell_number:.2e}')
            print(f'Average Final Cell Number: {B8_15x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Cell Number"]:.2e}\n')

        print('---BIOPROCESS RESULT---')
        print(f'CONFIDENCE LEVEL: {round(B8_15x_bioprocess.success_probability * 100, 1)} %')
        print(f'TOTAL DURATION: {B8_15x_bioprocess.total_duration} days')
        print(f'TOTAL EXECUTION TIME: {B8_15x_bioprocess.total_worktime} hours')
        print(f'TOTAL COST: {B8_15x_bioprocess.total_cost} €')
        print(f'Total Direct Cost: {B8_15x_bioprocess.total_direct_cost} €')
        print(f'Total Indirect Cost: {B8_15x_bioprocess.total_indirect_cost} €')
        print(f'Total Failure Cost: {B8_15x_bioprocess.total_failure_cost} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {B8_15x_bioprocess.final_cell_number:.2e}\n')

        print('---BIOPROCESS INFO---')
        print(f'''{B8_15x_bioprocess.cycle_info.to_string(formatters=["{:.0f}".format, "{:.1f}".format, "{:.2e}".format,
              "{:.2f}".format, "{:.2f}".format])}\n''')
        print(f'{B8_15x_bioprocess.cycle_cost_categories}\n')
        print(f'{B8_15x_bioprocess.cycle_cost_stages}\n')
        print('---COST CATEGORIES---')
        print(f'Consumables: {round(B8_15x_bioprocess.total_consumables_cost, 2)} €')
        print(f'Reagents: {round(B8_15x_bioprocess.total_reagents_cost, 2)} €')
        print(f'Facility: {round(B8_15x_bioprocess.total_facility_cost, 2)} €')
        print(f'Labor: {round(B8_15x_bioprocess.total_labor_cost, 2)} €\n')

        # B8 2.0x efficacy (and 5 days of expansion instead of 7)
        B8_cc.expansion_duration = 5
        B8_cc.fold_increase_average = without_DS_cc.fold_increase_average * 2
        B8_cc.fold_increase_std = without_DS_cc.fold_increase_std * 2

        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(dt, B8_cc)

        for index, value in enumerate(fold_expansion_distribution):
            fold_expansion_distribution[index] *= 2.0

        B8_20x_bioprocess = Bioprocess(dt, ft, B8_cc, fold_expansion_distribution,
                                recovery_efficiency_distribution, init_sim_num)

        print('---PROPOSED WORKFLOW (B8 2x)---')
        print('--[2D Expansion]--')
        print(f'{B8_20x_bioprocess.passaging_workflow}\n')
        print(f'DURATION: {B8_20x_bioprocess.cycle_info.at["2D Exp", "Duration (d)"]} days')
        print(f'EXECUTION TIME: {B8_20x_bioprocess.cycle_info.at["2D Exp", "Worktime (h)"]} hours')
        print(f'COST: {round(B8_20x_bioprocess.cycle_info.at["2D Exp", "Phase Cost (€)"], 2)} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {B8_20x_bioprocess.cycle_info.at["2D Exp", "Cell Number"]:.2e}\n')

        print('--[Bioreactor Expansion]--')
        for expansion_cycle, current_workflow in enumerate(B8_20x_bioprocess.bioreactor_workflow):
            print(f'<Cycle C{expansion_cycle}>')
            print(f'{current_workflow}\n')
            initial_cell_number = cc.seeding_density * B8_20x_bioprocess.optimal_working_volumes[expansion_cycle]
            print(f'DURATION: {B8_20x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Duration (d)"]} days')
            print(f'EXECUTION TIME: {B8_20x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Worktime (h)"]} hours')
            print(f'COST: {round(B8_20x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Phase Cost (€)"], 2)} €')
            print(f'Initial Cell Number: {initial_cell_number:.2e}')
            print(f'Average Final Cell Number: {B8_20x_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Cell Number"]:.2e}\n')

        print('---BIOPROCESS RESULT---')
        print(f'CONFIDENCE LEVEL: {round(B8_20x_bioprocess.success_probability * 100, 1)} %')
        print(f'TOTAL DURATION: {B8_20x_bioprocess.total_duration} days')
        print(f'TOTAL EXECUTION TIME: {B8_20x_bioprocess.total_worktime} hours')
        print(f'TOTAL COST: {B8_20x_bioprocess.total_cost} €')
        print(f'Total Direct Cost: {B8_20x_bioprocess.total_direct_cost} €')
        print(f'Total Indirect Cost: {B8_20x_bioprocess.total_indirect_cost} €')
        print(f'Total Failure Cost: {B8_20x_bioprocess.total_failure_cost} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {B8_20x_bioprocess.final_cell_number:.2e}\n')

        print('---BIOPROCESS INFO---')
        print(f'''{B8_20x_bioprocess.cycle_info.to_string(formatters=["{:.0f}".format, "{:.1f}".format, "{:.2e}".format,
              "{:.2f}".format, "{:.2f}".format])}\n''')
        print(f'{B8_20x_bioprocess.cycle_cost_categories}\n')
        print(f'{B8_20x_bioprocess.cycle_cost_stages}\n')
        print('---COST CATEGORIES---')
        print(f'Consumables: {round(B8_20x_bioprocess.total_consumables_cost, 2)} €')
        print(f'Reagents: {round(B8_20x_bioprocess.total_reagents_cost, 2)} €')
        print(f'Facility: {round(B8_20x_bioprocess.total_facility_cost, 2)} €')
        print(f'Labor: {round(B8_20x_bioprocess.total_labor_cost, 2)} €\n')

        # Bar Charts
        # create bar charts with cost distributions and save them in the respective folder
        bar_width = 0.15
        dots_per_inch = 600

        # Cost per Category
        # absolute cost distribution
        figure, axis = plt.subplots(tight_layout=True)

        axis.set_title('Impact of B8 efficacy on bioprocess costs by category', fontweight='bold')

        cost_categories = ('Consumables', 'Reagents', 'Facility', 'Labor')
        y_pos = np.arange(len(cost_categories))
        B8_075x_costs = [B8_075x_bioprocess.total_consumables_cost, B8_075x_bioprocess.total_reagents_cost,
                         B8_075x_bioprocess.total_facility_cost, B8_075x_bioprocess.total_labor_cost]
        mTeSR1_costs = [mTeSR1_bioprocess.total_consumables_cost, mTeSR1_bioprocess.total_reagents_cost,
                           mTeSR1_bioprocess.total_facility_cost, mTeSR1_bioprocess.total_labor_cost]
        B8_1x_costs = [B8_1x_bioprocess.total_consumables_cost, B8_1x_bioprocess.total_reagents_cost,
                       B8_1x_bioprocess.total_facility_cost, B8_1x_bioprocess.total_labor_cost]
        B8_15x_costs = [B8_15x_bioprocess.total_consumables_cost, B8_15x_bioprocess.total_reagents_cost,
                        B8_15x_bioprocess.total_facility_cost, B8_15x_bioprocess.total_labor_cost]
        B8_20x_costs = [B8_20x_bioprocess.total_consumables_cost, B8_20x_bioprocess.total_reagents_cost,
                        B8_20x_bioprocess.total_facility_cost, B8_20x_bioprocess.total_labor_cost]

        axis.barh(y_pos - bar_width * 2.2, B8_075x_costs, bar_width, align='center', color='rosybrown',
                  label='B8 (0.75x)')
        axis.barh(y_pos - bar_width * 1.1, mTeSR1_costs, bar_width, align='center', color='goldenrod',
                  label='mTeSR1')
        axis.barh(y_pos, B8_1x_costs, bar_width, align='center', color='indianred',
                  label='B8 (1x)')
        axis.barh(y_pos + bar_width * 1.1, B8_15x_costs, bar_width, align='center', color='lightcoral',
                  label='B8 (1.5x)')
        axis.barh(y_pos + bar_width * 2.2, B8_20x_costs, bar_width, align='center', color='lightsalmon',
                  label='B8 (2x)')
        axis.set_yticks(y_pos)
        axis.set_yticklabels(cost_categories)
        axis.invert_yaxis()
        axis.set_xlabel('Cost (€)')
        axis.legend()

        plt.savefig('studies/B8_comparison_absolute_cost_category.png', dpi=dots_per_inch)

        # Cost per Stage
        # absolute cost distribution
        axis.clear()

        axis.set_title('Impact of B8 efficacy on bioprocess costs by stage', fontweight='bold')

        cost_categories = ('Pre-Inoculation', 'Expansion', 'Quality Control', 'Harvesting')
        y_pos = np.arange(len(cost_categories))
        B8_075x_costs = [B8_075x_bioprocess.total_pre_inoculation_cost,
                       B8_075x_bioprocess.total_expansion_cost, B8_075x_bioprocess.total_quality_control_cost,
                       B8_075x_bioprocess.total_harvesting_cost]
        mTeSR1_costs = [mTeSR1_bioprocess.total_pre_inoculation_cost,
                        mTeSR1_bioprocess.total_expansion_cost, mTeSR1_bioprocess.total_quality_control_cost, mTeSR1_bioprocess.total_harvesting_cost]
        B8_1x_costs = [B8_1x_bioprocess.total_pre_inoculation_cost,
                       B8_1x_bioprocess.total_expansion_cost, B8_1x_bioprocess.total_quality_control_cost,
                       B8_1x_bioprocess.total_harvesting_cost]
        B8_15x_costs = [B8_15x_bioprocess.total_pre_inoculation_cost,
                       B8_15x_bioprocess.total_expansion_cost, B8_15x_bioprocess.total_quality_control_cost,
                       B8_15x_bioprocess.total_harvesting_cost]
        B8_20x_costs = [B8_20x_bioprocess.total_pre_inoculation_cost,
                       B8_20x_bioprocess.total_expansion_cost, B8_20x_bioprocess.total_quality_control_cost,
                       B8_20x_bioprocess.total_harvesting_cost]

        axis.barh(y_pos - bar_width * 2.2, B8_075x_costs, bar_width, align='center', color='rosybrown',
                  label='B8 (0.75x)')
        axis.barh(y_pos - bar_width * 1.1, mTeSR1_costs, bar_width, align='center', color='goldenrod',
                  label='mTeSR1 costs')
        axis.barh(y_pos, B8_1x_costs, bar_width, align='center', color='indianred',
                  label='B8 (1x)')
        axis.barh(y_pos + bar_width * 1.1, B8_15x_costs, bar_width, align='center', color='lightcoral',
                  label='B8 (1.5x)')
        axis.barh(y_pos + bar_width * 2.2, B8_20x_costs, bar_width, align='center', color='lightsalmon',
                  label='B8 (2x)')
        axis.legend()
        axis.set_yticks(y_pos)
        axis.set_yticklabels(cost_categories)
        axis.invert_yaxis()
        axis.set_xlabel('Cost (€)')

        plt.savefig('studies/B8_comparison_absolute_cost_stage.png', dpi=dots_per_inch)

        # Total Costs
        bar_width = 0.5
        axis.clear()

        axis.set_title('Impact of B8 efficacy on total bioprocess cost', fontweight='bold')

        labels = ['B8 (0.75x)', 'mTeSR1', 'B8 (1x)', 'B8 (1.5x)', 'B8 (2x)']
        colors = ['rosybrown', 'goldenrod', 'indianred', 'lightcoral', 'lightsalmon']
        y_pos = np.arange(len(labels))
        costs = [B8_075x_bioprocess.total_cost - B8_075x_bioprocess.total_medium_cost,
                 mTeSR1_bioprocess.total_cost - mTeSR1_bioprocess.total_medium_cost,
                 B8_1x_bioprocess.total_cost - B8_1x_bioprocess.total_medium_cost,
                 B8_15x_bioprocess.total_cost - B8_15x_bioprocess.total_medium_cost,
                 B8_20x_bioprocess.total_cost - B8_20x_bioprocess.total_medium_cost]

        axis.barh(y_pos, costs, bar_width, align='center', color=colors)
        axis.set_yticks(y_pos)
        axis.set_yticklabels(labels)
        axis.invert_yaxis()
        axis.set_xlabel('Total Cost (€)')

        # add medium cost seperately from other reagent costs
        medium_percentages = [B8_075x_bioprocess.total_medium_cost / B8_075x_bioprocess.total_cost * 100,
                              mTeSR1_bioprocess.total_medium_cost / mTeSR1_bioprocess.total_cost * 100,
                              B8_1x_bioprocess.total_medium_cost / B8_1x_bioprocess.total_cost * 100,
                              B8_15x_bioprocess.total_medium_cost / B8_15x_bioprocess.total_cost * 100,
                              B8_20x_bioprocess.total_medium_cost / B8_20x_bioprocess.total_cost * 100]
        reduction_percentages = [(1 - (B8_075x_bioprocess.total_cost / mTeSR1_bioprocess.total_cost)) * 100,
                                 (1 - (B8_1x_bioprocess.total_cost / mTeSR1_bioprocess.total_cost)) * 100,
                                 (1 - (B8_15x_bioprocess.total_cost / mTeSR1_bioprocess.total_cost)) * 100,
                                 (1 - (B8_20x_bioprocess.total_cost / mTeSR1_bioprocess.total_cost)) * 100]

        B8_075x_bar = axis.barh(y_pos[0], B8_075x_bioprocess.total_medium_cost, bar_width, left=costs[0],
                               color='peachpuff', label='B8 medium cost')
        mTeSR1_bar = axis.barh(y_pos[1], mTeSR1_bioprocess.total_medium_cost, bar_width,
                                  left=costs[1], color='gold', label='mTeSR1 medium cost')
        B8_1x_bar = axis.barh(y_pos[2], B8_1x_bioprocess.total_medium_cost, bar_width, left=costs[2],
                               color='peachpuff')
        B8_15x_bar = axis.barh(y_pos[3], B8_15x_bioprocess.total_medium_cost, bar_width, left=costs[3],
                               color='peachpuff')
        B8_20x_bar = axis.barh(y_pos[4], B8_20x_bioprocess.total_medium_cost, bar_width, left=costs[4],
                               color='peachpuff')
        axis.bar_label(B8_075x_bar, labels=['(%.0f%%)' % reduction_percentages[0]], padding=12)
        axis.bar_label(B8_1x_bar, labels=['(%.0f%%)' % reduction_percentages[1]], padding=12)
        axis.bar_label(B8_15x_bar, labels=['(%.0f%%)' % reduction_percentages[2]], padding=12)
        axis.bar_label(B8_20x_bar, labels=['(%.0f%%)' % reduction_percentages[3]], padding=12)

        axis.bar_label(B8_075x_bar, labels=['%.0f%%' % medium_percentages[0]], label_type='center',
                       fontweight='bold')
        axis.bar_label(mTeSR1_bar, labels=['%.0f%%' % medium_percentages[1]], label_type='center',
                       fontweight='bold')
        axis.bar_label(B8_1x_bar, labels=['%.0f%%' % medium_percentages[2]], label_type='center',
                       fontweight='bold')
        axis.bar_label(B8_15x_bar, labels=['%.0f%%' % medium_percentages[3]], label_type='center',
                       fontweight='bold')
        axis.bar_label(B8_20x_bar, labels=['%.0f%%' % medium_percentages[4]], label_type='center',
                       fontweight='bold')
        axis.legend()

        plt.savefig('studies/B8_comparison_total_cost.png', dpi=dots_per_inch)

    if value_st == 'qual' or value_st == 'general':
        # simulate bioprocesses relevant to the comparison of different quality control failure detection
        # probabilities
        without_DS_cc = CultureConditions(initnum=cc.initial_cell_number, targnum=cc.target_cell_number,
                                  medium='mTeSR1', initvol=cc.initial_volume, sim='Nogueira', DS=False)

        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(dt, without_DS_cc)

        # 20% failure chance / no intermediate quality control (0% failure detection probability)
        without_DS_cc.bioprocess_failure_prob = 0.2
        without_DS_cc.qual_ctrl_detection_prob = 0

        mTeSR1_20_0_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        intermediate_qual_ctrl_cost = (mTeSR1_20_0_bioprocess.cycle_cost_stages['Quality Control'].sum()
                                       - mTeSR1_20_0_bioprocess.cycle_cost_stages['Quality Control'].max())

        mTeSR1_20_0_bioprocess_cost = round(mTeSR1_20_0_bioprocess.total_cost - intermediate_qual_ctrl_cost, 2)

        # 30% failure chance / no intermediate quality control (0% failure detection probability)
        without_DS_cc.bioprocess_failure_prob = 0.3
        without_DS_cc.qual_ctrl_detection_prob = 0

        mTeSR1_30_0_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        intermediate_qual_ctrl_cost = (mTeSR1_30_0_bioprocess.cycle_cost_stages['Quality Control'].sum()
                                       - mTeSR1_30_0_bioprocess.cycle_cost_stages['Quality Control'].max())

        mTeSR1_30_0_bioprocess_cost = round(mTeSR1_30_0_bioprocess.total_cost - intermediate_qual_ctrl_cost, 2)

        # 40% failure chance / no intermediate quality control (0% failure detection probability)
        without_DS_cc.bioprocess_failure_prob = 0.4
        without_DS_cc.qual_ctrl_detection_prob = 0

        mTeSR1_40_0_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        intermediate_qual_ctrl_cost = (mTeSR1_40_0_bioprocess.cycle_cost_stages['Quality Control'].sum()
                                       - mTeSR1_40_0_bioprocess.cycle_cost_stages['Quality Control'].max())

        mTeSR1_40_0_bioprocess_cost = round(mTeSR1_40_0_bioprocess.total_cost - intermediate_qual_ctrl_cost, 2)

        # 20% failure chance / 10% failure detection probability
        without_DS_cc.bioprocess_failure_prob = 0.2
        without_DS_cc.qual_ctrl_detection_prob = 0.1

        mTeSR1_20_10_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # 20% failure chance / 20% failure detection probability
        without_DS_cc.qual_ctrl_detection_prob = 0.2

        mTeSR1_20_20_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # 20% failure chance / 30% failure detection probability
        without_DS_cc.qual_ctrl_detection_prob = 0.3

        mTeSR1_20_30_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # 20% failure chance / 40% failure detection probability
        without_DS_cc.qual_ctrl_detection_prob = 0.4

        mTeSR1_20_40_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # 30% failure chance / 10% failure detection probability
        without_DS_cc.bioprocess_failure_prob = 0.3
        without_DS_cc.qual_ctrl_detection_prob = 0.1

        mTeSR1_30_10_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # 30% failure chance / 20% failure detection probability
        without_DS_cc.qual_ctrl_detection_prob = 0.2

        mTeSR1_30_20_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # 30% failure chance / 30% failure detection probability
        without_DS_cc.qual_ctrl_detection_prob = 0.3

        mTeSR1_30_30_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # 30% failure chance / 40% failure detection probability
        without_DS_cc.qual_ctrl_detection_prob = 0.4

        mTeSR1_30_40_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # 40% failure chance / 10% failure detection probability
        without_DS_cc.bioprocess_failure_prob = 0.4
        without_DS_cc.qual_ctrl_detection_prob = 0.1

        mTeSR1_40_10_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # 40% failure chance / 20% failure detection probability
        without_DS_cc.qual_ctrl_detection_prob = 0.2

        mTeSR1_40_20_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # 40% failure chance / 30% failure detection probability
        without_DS_cc.qual_ctrl_detection_prob = 0.3

        mTeSR1_40_30_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # 40% failure chance / 40% failure detection probability
        without_DS_cc.qual_ctrl_detection_prob = 0.4

        mTeSR1_40_40_bioprocess = Bioprocess(dt, ft, without_DS_cc, fold_expansion_distribution,
                                  recovery_efficiency_distribution, init_sim_num)

        # Bar Charts
        # create bar charts with saved cost distributions and save them in the respective folder
        bar_width = 0.2
        dots_per_inch = 600

        figure, axis = plt.subplots(tight_layout=True)

        axis.set_title('Impact of quality control on total bioprocess cost', fontweight='bold')

        labels = ['20% batch failure rate', '30% batch failure rate', '40% batch failure rate']
        colors = ['rebeccapurple', 'slateblue', 'mediumpurple', 'plum']
        y_pos = np.arange(len(labels))
        fdp_10_costs = [-(mTeSR1_20_0_bioprocess_cost - mTeSR1_20_10_bioprocess.total_cost)
                        / mTeSR1_20_0_bioprocess_cost * 100,
                        -(mTeSR1_30_0_bioprocess_cost - mTeSR1_30_10_bioprocess.total_cost)
                        / mTeSR1_30_0_bioprocess_cost * 100,
                        -(mTeSR1_40_0_bioprocess_cost - mTeSR1_40_10_bioprocess.total_cost)
                        / mTeSR1_40_0_bioprocess_cost * 100]
        fdp_20_costs = [-(mTeSR1_20_0_bioprocess_cost - mTeSR1_20_20_bioprocess.total_cost)
                        / mTeSR1_20_0_bioprocess_cost * 100,
                        -(mTeSR1_30_0_bioprocess_cost - mTeSR1_30_20_bioprocess.total_cost)
                        / mTeSR1_30_0_bioprocess_cost * 100,
                        -(mTeSR1_40_0_bioprocess_cost - mTeSR1_40_20_bioprocess.total_cost)
                        / mTeSR1_40_0_bioprocess_cost * 100]
        fdp_30_costs = [-(mTeSR1_20_0_bioprocess_cost - mTeSR1_20_30_bioprocess.total_cost)
                        / mTeSR1_20_0_bioprocess_cost * 100,
                        -(mTeSR1_30_0_bioprocess_cost - mTeSR1_30_30_bioprocess.total_cost)
                        / mTeSR1_30_0_bioprocess_cost * 100,
                        -(mTeSR1_40_0_bioprocess_cost - mTeSR1_40_30_bioprocess.total_cost)
                        / mTeSR1_40_0_bioprocess_cost * 100]
        fdp_40_costs = [-(mTeSR1_20_0_bioprocess_cost - mTeSR1_20_40_bioprocess.total_cost)
                        / mTeSR1_20_0_bioprocess_cost * 100,
                        -(mTeSR1_30_0_bioprocess_cost - mTeSR1_30_40_bioprocess.total_cost)
                        / mTeSR1_30_0_bioprocess_cost * 100,
                        -(mTeSR1_40_0_bioprocess_cost - mTeSR1_40_40_bioprocess.total_cost)
                        / mTeSR1_40_0_bioprocess_cost * 100]

        print(f'!!!{-(mTeSR1_30_0_bioprocess_cost - mTeSR1_30_20_bioprocess.total_cost)/ mTeSR1_30_0_bioprocess_cost * 100}!!!')

        axis.barh(y_pos - bar_width * 1.65, fdp_10_costs, bar_width, align='center', color='rebeccapurple',
                  label='10% FDR')
        axis.barh(y_pos - bar_width * 0.55, fdp_20_costs, bar_width, align='center', color='slateblue',
                  label='20% FDR')
        axis.barh(y_pos + bar_width * 0.55, fdp_30_costs, bar_width, align='center', color='mediumpurple',
                  label='30% FDR')
        axis.barh(y_pos + bar_width * 1.65, fdp_40_costs, bar_width, align='center', color='plum',
                  label='40% FDR')
        axis.legend()
        axis.set_yticks(y_pos)
        axis.set_yticklabels(labels)
        axis.invert_yaxis()
        axis.set_xlim(-12, 12)
        axis.set_xlabel('Relative Change from No Quality Control Cost (%)')

        plt.savefig('studies/qual_crtl_comparison_total_cost.png', dpi=dots_per_inch)

    if value_st == 'density' or value_st == 'general':
        # simulate bioprocesses relevant to the comparison of high and low seeding densities
        # high seeding density
        high_density_cc = CultureConditions(initnum=cc.initial_cell_number, targnum=cc.target_cell_number,
                                  medium='mTeSR1', initvol=cc.initial_volume, sim='Nogueira', DS=False)

        #print(ft.overall_daily_facility_costs)

        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(dt, high_density_cc)

        high_density_bioprocess = Bioprocess(dt, ft, high_density_cc, fold_expansion_distribution,
                                     recovery_efficiency_distribution, init_sim_num)

        print('---PROPOSED WORKFLOW (High Density)---')
        print('--[2D Expansion]--')
        print(f'{high_density_bioprocess.passaging_workflow}\n')
        print(f'DURATION: {high_density_bioprocess.cycle_info.at["2D Exp", "Duration (d)"]} days')
        print(f'EXECUTION TIME: {high_density_bioprocess.cycle_info.at["2D Exp", "Worktime (h)"]} hours')
        print(f'COST: {round(high_density_bioprocess.cycle_info.at["2D Exp", "Phase Cost (€)"], 2)} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {high_density_bioprocess.cycle_info.at["2D Exp", "Cell Number"]:.2e}\n')

        print('--[Bioreactor Expansion]--')
        for expansion_cycle, current_workflow in enumerate(high_density_bioprocess.bioreactor_workflow):
            print(f'<Cycle C{expansion_cycle}>')
            print(f'{current_workflow}\n')
            initial_cell_number = (high_density_cc.seeding_density
                                       * high_density_bioprocess.optimal_working_volumes[expansion_cycle])
            print(f'DURATION: {high_density_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Duration (d)"]} days')
            print(f'WORKTIME: {high_density_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Worktime (h)"]} hours')
            print(f'COST: {round(high_density_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Phase Cost (€)"], 2)} €')
            print(f'Initial Cell Number: {initial_cell_number:.2e}')
            print(f'Final Cell Number: {high_density_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Cell Number"]:.2e}\n')

        print('---BIOPROCESS RESULT---')
        print(f'CONFIDENCE LEVEL: {round(high_density_bioprocess.success_probability * 100, 1)} %')
        print(f'TOTAL DURATION: {high_density_bioprocess.total_duration} days')
        print(f'TOTAL EXECUTION TIME: {high_density_bioprocess.total_worktime} hours')
        print(f'TOTAL COST: {high_density_bioprocess.total_cost} €')
        print(f'Total Direct Cost: {high_density_bioprocess.total_direct_cost} €')
        print(f'Total Indirect Cost: {high_density_bioprocess.total_indirect_cost} €')
        print(f'Total Failure Cost: {high_density_bioprocess.total_failure_cost} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {high_density_bioprocess.final_cell_number:.2e}\n')

        print('---BIOPROCESS INFO---')
        print(f'''{high_density_bioprocess.cycle_info.to_string(formatters=["{:.0f}".format, "{:.1f}".format, "{:.2e}".format,
              "{:.2f}".format, "{:.2f}".format])}\n''')
        print(f'{high_density_bioprocess.cycle_cost_categories}\n')
        print(f'{high_density_bioprocess.cycle_cost_stages}\n')
        print('---COST CATEGORIES---')
        print(f'Consumables: {round(high_density_bioprocess.total_consumables_cost, 2)} €')
        print(f'Reagents: {round(high_density_bioprocess.total_reagents_cost, 2)} €')
        print(f'Facility: {round(high_density_bioprocess.total_facility_cost, 2)} €')
        print(f'Labor: {round(high_density_bioprocess.total_labor_cost, 2)} €\n')

        # low seeding density
        low_density_cc = CultureConditions(initnum=cc.initial_cell_number, targnum=cc.target_cell_number,
                                 medium='mTeSR1', initvol=cc.initial_volume, sim='Borys', DS=False)

        ft.overall_daily_facility_costs = round(ft.overall_daily_facility_costs + 6
                     * dt.vwbs_base_units.at['PBS 3MAG', 'Cost (€)'] / dt.equipment_depreciation_period, 2)

        #print(ft.overall_daily_facility_costs)

        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(dt, low_density_cc)

        low_density_bioprocess = Bioprocess(dt, ft, low_density_cc, fold_expansion_distribution,
                                    recovery_efficiency_distribution, init_sim_num)

        print('---PROPOSED WORKFLOW (Low Density)---')
        print('--[2D Expansion]--')
        print(f'{low_density_bioprocess.passaging_workflow}\n')
        print(f'DURATION: {low_density_bioprocess.cycle_info.at["2D Exp", "Duration (d)"]} days')
        print(f'EXECUTION TIME: {low_density_bioprocess.cycle_info.at["2D Exp", "Worktime (h)"]} hours')
        print(f'COST: {round(low_density_bioprocess.cycle_info.at["2D Exp", "Phase Cost (€)"], 2)} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {low_density_bioprocess.cycle_info.at["2D Exp", "Cell Number"]:.2e}\n')

        print('--[Bioreactor Expansion]--')
        for expansion_cycle, current_workflow in enumerate(low_density_bioprocess.bioreactor_workflow):
            print(f'<Cycle C{expansion_cycle}>')
            print(f'{current_workflow}\n')
            initial_cell_number = (low_density_cc.seeding_density
                                       * low_density_bioprocess.optimal_working_volumes[expansion_cycle])
            print(f'DURATION: {low_density_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Duration (d)"]} days')
            print(f'WORKTIME: {low_density_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Worktime (h)"]} hours')
            print(f'COST: {round(low_density_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Phase Cost (€)"], 2)} €')
            print(f'Initial Cell Number: {initial_cell_number:.2e}')
            print(f'Final Cell Number: {low_density_bioprocess.cycle_info.at[f"C{expansion_cycle}", "Cell Number"]:.2e}\n')

        print('---BIOPROCESS RESULT---')
        print(f'CONFIDENCE LEVEL: {round(low_density_bioprocess.success_probability * 100, 1)} %')
        print(f'TOTAL DURATION: {low_density_bioprocess.total_duration} days')
        print(f'TOTAL EXECUTION TIME: {low_density_bioprocess.total_worktime} hours')
        print(f'TOTAL COST: {low_density_bioprocess.total_cost} €')
        print(f'Total Direct Cost: {low_density_bioprocess.total_direct_cost} €')
        print(f'Total Indirect Cost: {low_density_bioprocess.total_indirect_cost} €')
        print(f'Total Failure Cost: {low_density_bioprocess.total_failure_cost} €')
        print(f'Initial Cell Number: {cc.initial_cell_number:.2e}')
        print(f'Average Final Cell Number: {low_density_bioprocess.final_cell_number:.2e}\n')

        print('---BIOPROCESS INFO---')
        print(f'''{low_density_bioprocess.cycle_info.to_string(formatters=["{:.0f}".format, "{:.1f}".format, "{:.2e}".format,
              "{:.2f}".format, "{:.2f}".format])}\n''')
        print(f'{low_density_bioprocess.cycle_cost_categories}\n')
        print(f'{low_density_bioprocess.cycle_cost_stages}\n')
        print('---COST CATEGORIES---')
        print(f'Consumables: {round(low_density_bioprocess.total_consumables_cost, 2)} €')
        print(f'Reagents: {round(low_density_bioprocess.total_reagents_cost, 2)} €')
        print(f'Facility: {round(low_density_bioprocess.total_facility_cost, 2)} €')
        print(f'Labor: {round(low_density_bioprocess.total_labor_cost, 2)} €\n')

        # Bar Charts
        # create bar charts with cost distributions and save them in the respective folder
        bar_width = 0.25
        dots_per_inch = 600

        # Cost per Category
        # absolute cost distribution
        figure, axis = plt.subplots(tight_layout=True)

        axis.set_title('Impact of inoculation density on bioprocess costs by category', fontweight='bold')

        cost_categories = ('Consumables', 'Reagents', 'Facility', 'Labor')
        y_pos = np.arange(len(cost_categories))
        high_density_costs = [high_density_bioprocess.total_consumables_cost,
                 high_density_bioprocess.total_reagents_cost, high_density_bioprocess.total_facility_cost, high_density_bioprocess.total_labor_cost]
        low_density_costs = [low_density_bioprocess.total_consumables_cost,
                 low_density_bioprocess.total_reagents_cost, low_density_bioprocess.total_facility_cost, low_density_bioprocess.total_labor_cost]

        axis.barh(y_pos - bar_width * 0.55, high_density_costs, bar_width, align='center', color='goldenrod',
                  label='high density')
        axis.barh(y_pos + bar_width * 0.55, low_density_costs, bar_width, align='center', color='peru',
                  label='low density')
        axis.legend()
        axis.set_yticks(y_pos)
        axis.set_yticklabels(cost_categories)
        axis.invert_yaxis()
        axis.set_xlabel('Cost (€)')

        plt.savefig('studies/seeding_density_comparison_absolute_cost_category.png', dpi=dots_per_inch)

        # Cost per Stage
        # absolute cost distribution
        axis.clear()

        axis.set_title('Impact of inoculation density on bioprocess costs by stage', fontweight='bold')

        cost_categories = ('Pre-Inoculation', 'Expansion', 'Quality Control', 'Harvesting')
        y_pos = np.arange(len(cost_categories))
        high_density_costs = [high_density_bioprocess.total_pre_inoculation_cost,
                              high_density_bioprocess.total_expansion_cost,
                              high_density_bioprocess.total_quality_control_cost,
                              high_density_bioprocess.total_harvesting_cost]
        low_density_costs = [low_density_bioprocess.total_pre_inoculation_cost,
                             low_density_bioprocess.total_expansion_cost,
                             low_density_bioprocess.total_quality_control_cost,
                             low_density_bioprocess.total_harvesting_cost]
        axis.barh(y_pos - bar_width * 0.55, high_density_costs, bar_width, align='center', color='goldenrod',
                  label='high density')
        axis.barh(y_pos + bar_width * 0.55, low_density_costs, bar_width, align='center', color='peru',
                  label='low density')
        axis.legend()
        axis.set_yticks(y_pos)
        axis.set_yticklabels(cost_categories)
        axis.invert_yaxis()
        axis.set_xlabel('Cost (€)')

        plt.savefig('studies/seeding_density_comparison_absolute_cost_stage.png', dpi=dots_per_inch)

        # Total Costs
        bar_width = 0.6
        axis.clear()

        axis.set_title('Impact of inoculation density on total bioprocess cost', fontweight='bold')

        labels = ['High Density', 'Low Density']
        colors = ['red', 'goldenrod', 'red', 'peru', 'red']
        y_pos = [1, 2, 3, 4, 5]
        costs = [0, high_density_bioprocess.total_cost - high_density_bioprocess.total_medium_cost, 0,
                 low_density_bioprocess.total_cost - low_density_bioprocess.total_medium_cost, 0]

        axis.barh(y_pos, costs, bar_width, align='center', color=colors)
        axis.set_yticks([2, 4])
        axis.set_yticklabels(labels)
        axis.invert_yaxis()
        axis.set_xlabel('Total Cost (€)')

        # add medium cost seperately from other reagent costs
        medium_percentages = [high_density_bioprocess.total_medium_cost / high_density_bioprocess.total_cost * 100,
                              low_density_bioprocess.total_medium_cost / low_density_bioprocess.total_cost * 100]
        reduction_percentages = [(1 - (low_density_bioprocess.total_cost / high_density_bioprocess.total_cost))
                                 * 100]

        high_density_bar = axis.barh(y_pos[1], high_density_bioprocess.total_medium_cost, bar_width, left=costs[1],
                               color='gold', label='high density medium cost')
        low_density_bar = axis.barh(y_pos[3], low_density_bioprocess.total_medium_cost, bar_width,
                                  left=costs[3], color='peachpuff', label='low density medium cost')
        axis.bar_label(low_density_bar, labels=['(%.0f%%)' % reduction_percentages[0]], padding=12)
        axis.bar_label(high_density_bar, labels=['%.0f%%' % medium_percentages[0]], label_type='center',
                       fontweight='bold')
        axis.bar_label(low_density_bar, labels=['%.0f%%' % medium_percentages[1]], label_type='center',
                       fontweight='bold')
        axis.legend()

        plt.savefig('studies/seeding_density_comparison_total_cost.png', dpi=dots_per_inch)


# Print Program Header
print('\n>-------------------------------------<=>-------------------------------------<')
print('                              Welcome to BEMSCA.                              ')
print(' Type "run" to simulate an hiPSC expansion bioprocess with the default culture ')
print('    conditions or "study" to use one of the preset comparisons of different    ')
print('           experimental setups. For a list of commands type "help".           ')
print('>-------------------------------------<=>-------------------------------------<')

# Clear Previous Output
clear_folders()

# Initialize Data Table and Default Culture Conditions
data_table = DataTable()
facility = Facility(data_table)
culture_conditions = CultureConditions()

# Initialize Stochastic Distributions Relevant to Simulation Source
(fold_expansion_distribution, recovery_efficiency_distribution, init_sim_num) = simul_distributions(data_table,
                                                                                      culture_conditions)

# Start Input Loop
command = ''

while command != 'exit':

    command = input('\n>>> ')

    if command == 'conditions':
        conditions_command(culture_conditions)
    elif command == 'custom':
        culture_conditions = custom_command(culture_conditions)
        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(data_table, culture_conditions)
    elif command == 'data':
        data_command(data_table)
    elif command == 'default':
        culture_conditions = default_command()
        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(data_tabel, culture_conditions)
    elif command == 'help':
        help_command()
    elif command == 'run':
        run_command(data_table, facility, culture_conditions, fold_expansion_distribution,
                    recovery_efficiency_distribution, init_sim_num)
    elif command == 'source':
        culture_conditions = source_command(culture_conditions)
        (fold_expansion_distribution, recovery_efficiency_distribution,
         init_sim_num) = simul_distributions(data_table, culture_conditions)
    elif command == 'study':
        study_command(data_table, facility, culture_conditions)
    elif command != 'exit':
        print('Invalid command.')
