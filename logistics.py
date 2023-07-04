import math
import collections
import numpy as np
import pandas as pd
from random import shuffle
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


# Build Stochastic Distributions Relevant to Simulation Source
def simul_distributions(dt, cc):
    # create distributions for fold expansion and recovery efficiency according to simulation source, as well as
    # cell number at the end of an expansion phase for the chosen planar platform
    simulation_number = int(1e5)
    fold_expansion_distribution = []
    recovery_efficiency_distribution = []
    cell_number_2D_distribution = []

    for index in range(simulation_number):
        fold_expansion_distribution.append(cc.fold_increase())
        recovery_efficiency_distribution.append(cc.recovery_efficiency())
        cell_number_2D_distribution.append(cc.well_final_cell_number(dt))   
    
    # create histograms and save them in their respective folders
    bin_number = 200
    dots_per_inch = 600
    figure, axis = plt.subplots()
    
    # fold expansion
    figure.suptitle('Fold Expansion Distribution')
    axis.hist(fold_expansion_distribution, bin_number, color='firebrick')
    axis.set_xlabel('Fold Expansion')
    axis.set_ylabel('Frequency')
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=simulation_number))
    plt.savefig('simul/fold_expansion.png', dpi=dots_per_inch)
    
    # recovery efficiency
    axis.clear()
    figure.suptitle('Recovery Efficiency Distribution')
    axis.hist(recovery_efficiency_distribution, bin_number, color='firebrick')
    axis.set_xlabel('Recovery Efficiency')    
    axis.set_ylabel('Frequency')    
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=simulation_number))
    plt.savefig('simul/recovery_efficiency.png', dpi=dots_per_inch)
    
    # 2D platform cell number
    axis.clear()
    figure.suptitle('Cell Number per Planar Expansion Platform Distribution')
    axis.hist(cell_number_2D_distribution, bin_number, color='firebrick')
    axis.set_xlabel('Cell Number per Planar Expansion Platform')    
    axis.set_ylabel('Frequency')    
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=simulation_number))
    plt.savefig('simul/cell_number_2Dplatform_.png', dpi=dots_per_inch)     
    
    return fold_expansion_distribution, recovery_efficiency_distribution, simulation_number


def well_amount(dt, cc, initial_well_amount, inoculation_cell_number):
    # determine the necessary amount of wells to obtain enough cells for inoculation taking stochasticity of
    # cell growth into account
    necessary_well_amount = 0
    current_cell_number = 0
    while current_cell_number < inoculation_cell_number:
        necessary_well_amount += 1
        current_cell_number += math.ceil(cc.well_final_cell_number(dt))

    # estimate the necessary expansion factor to reach the necessary well amount from the initial well amount
    expansion_factor = math.ceil(necessary_well_amount / initial_well_amount)

    # calculate how many passages are necessary to achieve the required expansion factor
    necessary_passage_number = 0
    while expansion_factor > dt.maximum_passaging_ratio ** necessary_passage_number:
        necessary_passage_number += 1  

    # determine the ideal combination of passaging ratios to obtain the required expansion factor
    aux1 = expansion_factor
    necessary_passages = 0
    while aux1 > 1:
        for passaging_ratio in reversed(range(2, dt.maximum_passaging_ratio + 1)):
            aux2 = aux1
            if aux1 % passaging_ratio == 0:
                aux1 /= passaging_ratio
                necessary_passages += 1
                break

        # if aux1 is not divisible by a possible passaging ratio, or if the amount of passages excedes
        # the necessary passage number, increment the expansion factor and restart the process
        if aux1 == aux2 or necessary_passages > necessary_passage_number:
            expansion_factor += 1
            aux1 = expansion_factor
            necessary_passages = 0

    # take into account the surplus wells that are required to satisfy the combination of passaging ratios
    while necessary_well_amount < expansion_factor * initial_well_amount:            
        necessary_well_amount += 1
        current_cell_number += math.ceil(cc.well_final_cell_number(dt))
            
    return necessary_well_amount, current_cell_number


def workflow_calculations(dt, cc):
    # Well Amount
    # determine maximum necessary well amount to achieve cell number for the initial inoculation
    simulation_number = int(1e3)
    well_amount_distribution = []
    cell_number_distribution = []
    
    for simul in range(simulation_number):
        (necessary_well_amount, cell_number) = well_amount(dt, cc, cc.initial_well_amount,
                                                           cc.inoculation_cell_number)
        well_amount_distribution.append(necessary_well_amount)
        cell_number_distribution.append(cell_number)
        
    occurrences = collections.Counter(well_amount_distribution)
    keys = occurrences.keys()
    necessary_well_amount = max(keys)
    
    average_wells_cell_number = math.ceil(sum(cell_number_distribution) / len(cell_number_distribution))
    
    # create bar chart for the necessary well amount and save it in its respective folder
    bin_number = 200
    bar_width = 0.5
    dots_per_inch = 600
    figure, axis = plt.subplots()
    
    # well amount
    possible_well_amounts = {'well amount': [], 'frequency': []}
    for key in keys:
        possible_well_amounts['well amount'].append(str(key))
        
    possible_well_amounts['well amount'].sort()
    for possible_well_amount in possible_well_amounts['well amount']:
        possible_well_amounts['frequency'].append(occurrences[int(possible_well_amount)])
    
    figure.suptitle('Necessary Well Amount to Inoculate Initial Volume')
    axis.bar(possible_well_amounts['well amount'], possible_well_amounts['frequency'], bar_width,
             color='firebrick')
    axis.set_xlabel('Necessary Well Amount')    
    axis.set_ylabel('Frequency')
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=simulation_number))
    plt.savefig('workflows/well_amount.png', dpi=dots_per_inch)    
    
    # Bioreactor Volumes
    simulation_number = int(1e5)
    desired_total_fold_increase = cc.target_cell_number / cc.inoculation_cell_number
    desired_fold_increase_3MAG = dt.vwbs_parameters.at['PBS 3MAG', 'Min Vol (mL)'] / cc.initial_volume
    lower_limit_fold_increase_3MAG = dt.vwbs_parameters.at['PBS 3MAG', 'Max Vol (mL)'] / cc.initial_volume
    upper_limit_fold_increase_3MAG = 2 * dt.vwbs_parameters.at['PBS 3MAG', 'Min Vol (mL)'] / cc.initial_volume
    
    initial_cycles_weights = [0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.02, 1.04, 1.06, 1.08, 1.1]
    for cycle_number in range(1, 11):
        for initial_cycles_weight in initial_cycles_weights:
            desired_cycle_fold_increase = desired_total_fold_increase**(1 / cycle_number)
            desired_cycle_fold_increases = []
            if cycle_number == 1:
                desired_cycle_fold_increases.append(desired_cycle_fold_increase)
            else:
                for cycle in range(cycle_number - 1):
                    desired_cycle_fold_increases.append(desired_cycle_fold_increase * initial_cycles_weight)
                desired_cycle_fold_increases.append(desired_cycle_fold_increase
                                                    * initial_cycles_weight**(-(cycle_number - 1)))

                # make sure the bioreactor inoculated for the last cycle is a PBS 3MAG
                penultimate_cum_fold_increase = 1
                for desired_fold_increase in desired_cycle_fold_increases[:-1]:
                    penultimate_cum_fold_increase *= desired_fold_increase

                if penultimate_cum_fold_increase < desired_fold_increase_3MAG:
                    desired_cycle_fold_increase = desired_fold_increase_3MAG**(1 / (cycle_number - 1))
                    for index in range(cycle_number - 1):
                        desired_cycle_fold_increases[index] = desired_cycle_fold_increase
                    desired_cycle_fold_increases[-1] = desired_total_fold_increase / desired_fold_increase_3MAG
                
                # consider the volume gap for which it is impossible to inoculate 2 PBS 3MAG
                if lower_limit_fold_increase_3MAG < penultimate_cum_fold_increase < upper_limit_fold_increase_3MAG:
                    desired_cycle_fold_increase = upper_limit_fold_increase_3MAG**(1 / (cycle_number - 1))
                    for index in range(cycle_number - 1):
                        desired_cycle_fold_increases[index] = desired_cycle_fold_increase
                    desired_cycle_fold_increases[-1] = desired_total_fold_increase / upper_limit_fold_increase_3MAG
            
            optimal_working_volumes = [cc.initial_volume]
            rounded_optimal_working_volumes = [cc.initial_volume] 
            for index, cycle_fold_increase in enumerate(desired_cycle_fold_increases[:-1]):
                optimal_working_volumes.append(optimal_working_volumes[index] * cycle_fold_increase)
                rounded_optimal_working_volumes.append(math.ceil(round(optimal_working_volumes[index]
                                                                       * cycle_fold_increase / 10, 2)) * 10)
                
            for index, volume in enumerate(rounded_optimal_working_volumes[1:]):
                desired_cycle_fold_increases[index] = volume / rounded_optimal_working_volumes[index]

            cycle_cell_numbers = []
            for cycle in range(cycle_number):
                cycle_cell_numbers.append([])               
                
            failures = 0                
            for simul in range(simulation_number):
                fold_increases = []
                cell_number = cc.seeding_density * cc.initial_volume 
                for cycle in range(cycle_number):
                    fold_increases.append(cc.fold_increase() * cc.recovery_efficiency())

                lag = 0
                for index, fold_increase in enumerate(fold_increases):
                    cell_number *= fold_increases[index]
                    cycle_cell_numbers[index].append(cell_number)
                    cell_number /= fold_increases[index]
                    
                    if lag > 0:
                        cum_fold_increase = 1
                        desired_fold_increase = 1
                        for i in range(lag + 1):
                            cum_fold_increase *= fold_increases[index - i]
                            desired_fold_increase *= desired_cycle_fold_increases[index - i]
                        if cum_fold_increase >= desired_fold_increase:
                            fold_increases[index] = desired_fold_increase / (cum_fold_increase / fold_increase)
                            lag = 0
                        else:
                            lag += 1

                    elif fold_increase >= desired_cycle_fold_increases[index]:
                        if index != (cycle_number - 1):
                            fold_increases[index] = desired_cycle_fold_increases[index]
                    else:
                        lag += 1
                    
                    cell_number *= fold_increases[index]

                if lag > 0:
                    failures += 1
                    
#                print(fold_increases, initial_cycles_weight, cell_number, failures)
            
            cycle_average_cell_numbers = [sum(cycle_cell_number) / simulation_number 
                                          for cycle_cell_number in cycle_cell_numbers]

            success_probability = (simulation_number - failures) / simulation_number
            
#            print('--')
#            print(initial_cycles_weight)
#            print(desired_cycle_fold_increases)
#            print(success_probability)
            
            if success_probability >= cc.threshold:
                success_probability = round(round(((simulation_number - failures) / simulation_number) / 0.001, 0)
                                        * 0.001, 3)
                break
                
        if success_probability >= cc.threshold:
            success_probability = round(round(((simulation_number - failures) / simulation_number) / 0.001, 0)
                                    * 0.001, 3)            
#            print(cycle_number, success_probability)
#            print(initial_cycles_weight)
#            print(desired_cycle_fold_increases)
#            print(rounded_optimal_working_volumes)            
            break
        
    return (necessary_well_amount, average_wells_cell_number, rounded_optimal_working_volumes, success_probability,
           cycle_average_cell_numbers, cycle_cell_numbers[-1])


def create_2D_workflow(dt, initial_well_amount, necessary_well_amount):
    passaging_workflow = {'Passaging Ratios': [], 'Passage Well Amount': []}
    
    # determine the ideal combination of passaging ratios to obtain the required expansion factor
    expansion_factor = math.ceil(necessary_well_amount / initial_well_amount)
    while expansion_factor > 1:
        for passaging_ratio in reversed(range(2, dt.maximum_passaging_ratio + 1)):
            if expansion_factor % passaging_ratio == 0:
                expansion_factor /= passaging_ratio
                passaging_workflow['Passaging Ratios'].append(passaging_ratio)
                break

    # begin the process with the lower passaging ratios
    passaging_workflow['Passaging Ratios'].sort()

    # determine the well amount used for each passage of hiPSCs according to previously determined passaging
    # ratios
    passage_index = []
    passage_number = 0
    current_well_amount = initial_well_amount
    for passage in passaging_workflow['Passaging Ratios']:
        current_well_amount *= passage
        passaging_workflow['Passage Well Amount'].append(current_well_amount)
        passage_number += 1
        passage_index.append(f'P{passage_number}')

    # add initial well amount to the passaging workflow
    passage_index.insert(0, 'P0')
    passaging_workflow['Passaging Ratios'].insert(0, np.nan)
    passaging_workflow['Passage Well Amount'].insert(0, initial_well_amount)

    # convert to pandas dataframe
    passaging_workflow = pd.DataFrame(passaging_workflow, index=passage_index)
    
    return passage_number, passaging_workflow


def create_3D_workflow(dt, optimal_working_volumes):
    expansion_cycles = 0
    bioreactor_workflow = []
    for index, optimal_working_volume in enumerate(optimal_working_volumes):
        volume = optimal_working_volume
        
        # determine the combinations of vwbs for the optimal working volume defined for each expansion cycle
        vwbs_types = pd.DataFrame({'Amount': [0, 0, 0], 'Working Volume (mL)': [0, 0, 0],
                                   'Associated Costs (€)': [.0, .0, .0]}, index=dt.vwbs)

        possible_3MAG = math.floor(volume / dt.vwbs_parameters.at['PBS 3MAG', 'Min Vol (mL)'])
        if possible_3MAG > 0:
            fraction_3MAG = volume / dt.vwbs_parameters.at['PBS 3MAG', 'Max Vol (mL)']
            minimum_3MAG = math.ceil(fraction_3MAG)

            # take into account that 3 PBS 0.5MAG are less costly than an extra PBS 3MAG (except for last cycle)
            if index != (len(optimal_working_volumes) - 1):
                surplus_fraction = fraction_3MAG % 1 * dt.vwbs_parameters.at['PBS 3MAG', 'Max Vol (mL)']
                if (surplus_fraction != 0 and surplus_fraction <= 3
                 * dt.vwbs_parameters.at['PBS 0.5MAG', 'Max Vol (mL)']):
                    minimum_3MAG -= 1
                    fraction_3MAG = minimum_3MAG

            vwbs_types.at['PBS 3MAG', 'Amount'] += minimum_3MAG
            vwbs_types.at['PBS 3MAG', 'Associated Costs (€)'] += (minimum_3MAG
                                                             * dt.vwbs_parameters.at['PBS 3MAG', 'Cost (€)'])
            volume_per_vwb = math.ceil(dt.vwbs_parameters.at['PBS 3MAG', 'Max Vol (mL)']
                                         * (fraction_3MAG / minimum_3MAG)) 
            vwbs_types.at['PBS 3MAG', 'Working Volume (mL)'] = volume_per_vwb
            volume -= volume_per_vwb * minimum_3MAG

        possible_05MAG = math.floor(volume / dt.vwbs_parameters.at['PBS 0.5MAG', 'Min Vol (mL)'])
        if possible_05MAG > 0:
            fraction_05MAG = volume / dt.vwbs_parameters.at['PBS 0.5MAG', 'Max Vol (mL)']
            minimum_05MAG = math.ceil(fraction_05MAG)

            # take into account that 1 PBS 0.1MAG is less costly than an extra PBS 0.5MAG
            surplus_fraction = fraction_05MAG % 1 * dt.vwbs_parameters.at['PBS 0.5MAG', 'Max Vol (mL)']
            if surplus_fraction != 0 and surplus_fraction <= dt.vwbs_parameters.at['PBS 0.1MAG', 'Max Vol (mL)']:
                minimum_05MAG -= 1
                fraction_05MAG = minimum_05MAG

            vwbs_types.at['PBS 0.5MAG', 'Amount'] = minimum_05MAG
            vwbs_types.at['PBS 0.5MAG', 'Associated Costs (€)'] = (minimum_05MAG
                                                                   * dt.vwbs_parameters.at['PBS 0.5MAG', 'Cost (€)'])
            volume_per_vwb = math.ceil(dt.vwbs_parameters.at['PBS 0.5MAG', 'Max Vol (mL)']
                                         * (fraction_05MAG / minimum_05MAG)) 
            vwbs_types.at['PBS 0.5MAG', 'Working Volume (mL)'] = volume_per_vwb
            volume -= volume_per_vwb * minimum_05MAG

            # take into account that there must be at least 60 mL left for the inoculation of a PBS 0.1MAG
            if volume > 0 and volume < 60:
                volume_deficit = 60 - volume
                volume_per_vwb -= math.ceil(volume_deficit / vwbs_types.at['PBS 0.5MAG', 'Amount'])
                vwbs_types.at['PBS 0.5MAG', 'Working Volume (mL)'] = volume_per_vwb
                volume = 60 + math.floor(vwbs_types.at['PBS 0.5MAG', 'Amount'] / 2)

        # take into account that the inoculation of a PBS 0.1MAG requires a minimum of 60 mL
        if volume > 0 and volume < 60:                
            volume = 60

        possible_01MAG = math.floor(volume / dt.vwbs_parameters.at['PBS 0.1MAG', 'Min Vol (mL)'])
        if possible_01MAG > 0:
            fraction_01MAG = volume / dt.vwbs_parameters.at['PBS 0.1MAG', 'Max Vol (mL)']
            minimum_01MAG = math.ceil(fraction_01MAG)
            vwbs_types.at['PBS 0.1MAG', 'Amount'] = minimum_01MAG
            vwbs_types.at['PBS 0.1MAG', 'Associated Costs (€)'] = (minimum_01MAG
                                                              * dt.vwbs_parameters.at['PBS 0.1MAG', 'Cost (€)'])
            volume_per_vwb = math.ceil(dt.vwbs_parameters.at['PBS 0.1MAG', 'Max Vol (mL)']
                                         * (fraction_01MAG / minimum_01MAG))

            # take into account that the inoculation of a second PBS 0.1MAG requires a minimum of 60 mL
            if volume_per_vwb < 60:
                volume_per_vwb = 60

            vwbs_types.at['PBS 0.1MAG', 'Working Volume (mL)'] = volume_per_vwb
            volume -= volume_per_vwb * minimum_01MAG

        expansion_cycles += 1
        bioreactor_workflow.append(vwbs_types)
    
    return expansion_cycles, bioreactor_workflow