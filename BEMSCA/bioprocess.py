import math
import numpy as np
import numpy.random as nprand
import pandas as pd
import utils


# -----------------------------------------------------------------------------
#    GLOBAL CONSTANTS
# -----------------------------------------------------------------------------
INI_QC_CELLS = 4e6 # number of cells used for initial quality control
INT_FC_ANTIBODIES = ["OCT4", "SOX2"]
INT_IMMUNO_ANTIBODIES = ["OCT4", "SOX2"]
ROCKI = "Y_27632"
SUR_FC_ANTIBODIES = ["TRA-1-60-PE", "SSEA-4-PE"]
SUR_IMMUNO_ANTIBODIES = ["TRA-1-60", "SSEA-4"]
YEAR_TO_DAYS = 365.25 # days


# -----------------------------------------------------------------------------
#    CLASSES
# -----------------------------------------------------------------------------
# Define class for simulation of planar expansion
class PlanarExpansion():
    # Initializer of class object
    def __init__(self, db_data, bio_params, d_facility_cost, d_labor_cost):
        # <>---------------- Class Constants ----------------<>        
        self.MAX_PASS_RATIO = 6
        self.MIN_PASS_RATIO = 3 
        self.N_WASHES = 3   
        self.PASS_DURATION = 4 # days
        self.RECOVERY_MEDIUM = 0.01 # L

        # <>---------- Important Object Attributes ----------<>
        self.culture_medium = db_data["Expansion Simulations"].loc[bio_params["Expansion Simulation"],
                                                                     "Culture Medium"].item()

        self.seeding_density = db_data["Expansion Simulations"].loc[bio_params["Expansion Simulation"],
                                                            "Seeding Density"].item()
        
        self.min_bioreactor_volume = db_data["Bioreactors"].loc[bio_params["Bioreactors"].item().split(";"),
                                                                "Min Volume"].min().item() * 1e3        

        self.planar_platform_data = db_data["2D Platforms"].loc[bio_params["2D Platform"].item()]

        self.coating_substrate = bio_params["Coating Substrate"].item()

        self.dissociation_enzyme = bio_params["Dissociation Enzyme"].item()

        # <>------------------- Main Body -------------------<>
        # Calculate number of cells necessary for bioreactor inoculation (seeding density * bioreactor volume)
        self.inoc_cells = self.seeding_density * self.min_bioreactor_volume

        # Determine optimal planar expansion workflow
        self.planar_workflow = self.determine_planar_workflow()

        # Determine total surfaces and duration of optimal workflow
        self.total_surfaces = sum(self.planar_workflow)
        self.duration = len(self.planar_workflow) * self.PASS_DURATION

        # Determine planar expansion cost
        self.determine_planar_expansion_cost(db_data, d_facility_cost, d_labor_cost)


    # Function for determining the optimal planar expansion workflow
    def determine_planar_workflow(self):
        # Add flow cytometry cells to target cell number for planar expansion (3 million)
        target_cells = self.inoc_cells + INI_QC_CELLS

        # Assuming that planar platform confluency follows a normal distribution:
        # Let desired_cells = X, surface_confluency = u, confluency_std = std and surface_number = s
        # z = (X - s*u) / s*std, where z = -3 implies that P(cells >= X) = 99.85%
        # So s = X/(u - 3*std) gives a high certainty that enough cells are obtained for bioreactor inoculation
        
        # Determine minimum number of surfaces of selected 2D platforms required for bioreactor inoculation
        surface_confluency = self.planar_platform_data["Surface Confluency"].item()
        confluency_std = self.planar_platform_data["Confluency STD"].item()

        final_surfaces = math.ceil(target_cells / (surface_confluency - 3*confluency_std))

        # Determine minimum number of passages required to obtain target cell number
        passages = math.ceil(math.log(final_surfaces, self.MAX_PASS_RATIO))

        # Determine optimal planar expansion workflow (number of surfaces for each passage)
        # It is more economical to leave large passage ratios towards the end of planar expansion
        planar_workflow = [final_surfaces]
        surfaces = final_surfaces
        for passage in range(1, passages):
            surfaces = math.ceil(surfaces/self.MAX_PASS_RATIO)
            if surfaces < self.MIN_PASS_RATIO:
                surfaces = self.MIN_PASS_RATIO

            planar_workflow.insert(0, surfaces)

        # Cells are initially thawed onto a single surface
        planar_workflow.insert(0, 1)

        return planar_workflow
    

    # Function for determining planar expansion cost
    def determine_planar_expansion_cost(self, db_data, d_facility_cost, d_labor_cost):
        # <>---------------- Consumables Costs ---------------<>
        # Calculate platform cost based on number of necessary surfaces
        N_platforms = 0
        for surfaces in self.planar_workflow:
            N_platforms += math.ceil(surfaces / self.planar_platform_data["Surfaces"].item())
        T_platform_cost = N_platforms * self.planar_platform_data["Cost"].item()
               
        # <>----------------- Reagents Costs -----------------<>
        # Calculate total coating volume, along with corresponding costs
        T_coating_volume = self.planar_platform_data["Coating Volume"].item() * self.total_surfaces
        T_coating_cost = T_coating_volume * db_data["Reagents"].loc[self.coating_substrate].item()        
        
        # Calculate total EDTA volume along with corresponding costs
        T_EDTA_volume = (self.planar_platform_data["Washing Volume"].item()
                            * (self.total_surfaces - self.planar_workflow[-1]) * self.N_WASHES)
        T_EDTA_cost = T_EDTA_volume * db_data["Reagents"].loc["EDTA"].item()

        # Calculate total medium volume, including medium used for thawing cells, along with corresponding cost
        T_medium_volume = self.planar_platform_data["Culture Volume"].item() * self.total_surfaces * self.PASS_DURATION
        T_medium_volume += self.RECOVERY_MEDIUM
        self.T_medium_cost = T_medium_volume * db_data["Reagents"].loc[self.culture_medium].item()

        # Calculate total ROCK Inhibitor cost (added to medium for 1h before cell dissociation)
        T_ROCKi_volume = self.planar_platform_data["Culture Volume"].item() * self.planar_workflow[-1]
        T_ROCKi_cost = T_ROCKi_volume * db_data["Reagents"].loc[ROCKI].item()

        # Calculate DPBS volume along with corresponding costs (used for washing before dissociation enzyme)
        T_DPBS_volume = self.planar_platform_data["Washing Volume"].item() * self.planar_workflow[-1]
        T_DPBS_cost = T_DPBS_volume * db_data["Reagents"].loc["DPBS"].item()

        # Calculate total dissociation enzyme volume along with corresponding cost
        T_diss_enz_volume = self.planar_platform_data["Washing Volume"].item() * self.planar_workflow[-1]
        T_diss_enz_cost = T_diss_enz_volume * db_data["Reagents"].loc[self.dissociation_enzyme].item()

        # <>--------------- Costs by Category ---------------<>
        # Calculate total category costs of planar expansion
        self.T_consumables_cost = T_platform_cost
        self.T_reagents_cost = (T_coating_cost + T_EDTA_cost + self.T_medium_cost
                                + T_ROCKi_cost + T_DPBS_cost + T_diss_enz_cost)
        self.T_facility_cost = d_facility_cost * self.duration
        self.T_labor_cost = d_labor_cost * self.duration
        
        # Calculate overall planar expansion cost
        self.overall_cost = self.T_consumables_cost + self.T_reagents_cost + self.T_facility_cost + self.T_labor_cost


# Define class for simulation of bioreactor expansion
class BioreactorExpansion():
    # Initializer of class object
    def __init__(self, db_data, bio_params, d_facility_cost, d_labor_cost):
        # <>---------------- Class Constants ----------------<>        
        self.DISS_ENZ_VOL_RATIO = 0.2
        self.DECREASE_RATIO = 0.99
        self.FIN_QUAL_DURATION = 3 # days
        self.MIN_FINAL_VOLUME = 1.8 # L
        self.SIMULATION_RUNS = int(1e5)

        # <>---------- Important Object Attributes ----------<>
        self.expansion_simulation = db_data["Expansion Simulations"].loc[bio_params["Expansion Simulation"]]

        self.recovery_simulation = db_data["Recovery Simulations"].loc[bio_params["Recovery Simulation"]]

        self.bioreactors = db_data["Bioreactors"].loc[bio_params["Bioreactors"].item().split(";")]

        self.culture_medium = self.expansion_simulation["Culture Medium"].item()

        self.coating_substrate = bio_params["Coating Substrate"].item()

        self.dissociation_enzyme = bio_params["Dissociation Enzyme"].item()        

        # <>------------------- Main Body -------------------<>
        # Calculate number of cells at first bioreactor inoculation (seeding density * bioreactor volume)
        self.initial_cells = self.expansion_simulation["Seeding Density"] * self.bioreactors["Min Volume"].min() * 1e3

        # Calculate desired total fold increase
        self.tfi = bio_params["Target Cell Number"].item() / self.initial_cells.item()

        # Determine optimal bioreactor expansion workflow ([0] -> cycle medium volumes, [1] -> required bioreactors)
        self.bioreactor_workflow = self.determine_bioreactor_workflow(bio_params)

        # Determine bioreactor expansion duration
        self.duration = len(self.bioreactor_workflow[0]) * self.expansion_simulation["Bioreactor Culture Time"].item()

        # Add final quality control duration
        self.duration += self.FIN_QUAL_DURATION

        # Determine bioreactor expansion cost
        self.determine_bioreactor_expansion_cost(db_data, d_facility_cost, d_labor_cost)


    # Function for determining the optimal bioreactor expansion workflow
    def determine_bioreactor_workflow(self, bio_params):
        # Calculate std of fold expansion and recovery efficiency distributions
        self.fold_exp_std = utils.sem_to_std(self.expansion_simulation["Fold Expansion SEM"].item(),
                                                    self.expansion_simulation["Experiment Sample Size"].item())
        
        self.recovery_eff_std = utils.sem_to_std(self.recovery_simulation["Recovery Efficiency SEM"].item(),
                                                        self.recovery_simulation["Experiment Sample Size"].item())
        
        # !!! Recovery efficiency std was divided by 3 in order to create a beta distribution that makes sense
        # (not overly weighted towards 1). This was done after observing that the original paper reports a sem
        # which is too high when considering the reported average and that the possible values range from 0 to 1. If
        # other data is used the following line of code should not be replicated:
        self.recovery_eff_std /= 3    
        
        # Calculate alpha and beta parameters of recovery efficiency beta distribution
        (alpha, beta) = utils.alpha_beta(self.recovery_simulation["Recovery Efficiency AVG"].item(),
                                                self.recovery_eff_std)

        # Check how many bioreactor expansion cycles are required to obtain the target cell number
        # while respecting the minimum threshold
        required_cycles = 1
        while True:
            fold_increase_pd = []
            # Simulate total fold increase for the stipulated number of simulation runs
            for s in range(0, self.SIMULATION_RUNS):
                fold_increase = 1
                for cycle in range(1, required_cycles+1):
                    # For each cycle draw random samples for fold expansion and recovery efficiency and calculate
                    # cycle fold increase (fold expansion x recovery efficiency)
                    fold_exp = nprand.normal(self.expansion_simulation["Fold Expansion AVG"].item(), self.fold_exp_std)
                    recovery_eff = nprand.beta(alpha, beta)
                    fold_increase *= fold_exp * recovery_eff
                
                fold_increase_pd.append(fold_increase)

            # Create numpy array with distribution of simulated fold increases
            fold_increase_pd = np.array(fold_increase_pd)

            # End while loop if minimum threshold is respected
            if np.percentile(fold_increase_pd, (1 - bio_params["Minimum Threshold"]) * 100) >= self.tfi:
                break

            # Increment number of bioreactor expansion cycles in workflow if minimum threshold is not respected
            required_cycles += 1

        # Optimize the volume of medium used in each bioreactor expansion cycle by avoiding the production
        # of surplus cells (limiting the target fold increase of each cycle), for as long as the minimum
        # threshold continues to be respected

        # Consider the initial optimal fold increase to be the average of the fold increase distribution
        optimal_fold_increase = (self.expansion_simulation["Fold Expansion AVG"].item()
                                 * self.recovery_simulation["Recovery Efficiency AVG"].item())
        
        # Establish a minimum fold increase if a minimum final volume is desired (this is to ensure that
        # the final expansion cycle takes place in a desired bioreactor type, for example)
        min_fold_increase = (self.MIN_FINAL_VOLUME /
                                 self.bioreactors["Min Volume"].min())**(1/(required_cycles-1))
        
        self.fold_increase_pds = []
        while True:
            fold_increase_pd = []
            # Simulate total fold increase for the stipulated number of simulation runs
            for s in range(0, self.SIMULATION_RUNS):
                fold_increase = 1
                for cycle in range(1, required_cycles+1):
                    # For each cycle draw random samples for fold expansion and recovery efficiency and calculate
                    # cycle fold increase (fold expansion x recovery efficiency)
                    fold_exp = nprand.normal(self.expansion_simulation["Fold Expansion AVG"].item(), self.fold_exp_std)
                    recovery_eff = nprand.beta(alpha, beta)
                    fold_increase *= fold_exp * recovery_eff
                    
                    # Limit fold increase of a cycle if it surpasses the established optimal fold increase
                    # (reduce medium usage of next cycle). This does not apply to the last cycle
                    if cycle != required_cycles and fold_increase >= optimal_fold_increase**cycle:
                        fold_increase = optimal_fold_increase**cycle
                
                fold_increase_pd.append(fold_increase)

            # Create numpy array with distribution of simulated fold increases
            fold_increase_pd = np.array(fold_increase_pd)

            # End while loop if minimum threshold is disrespected
            if np.percentile(fold_increase_pd, (1 - bio_params["Minimum Threshold"]) * 100) < self.tfi:
                # Undo decrease of optimal fold increase since this has led to disrespecting the minimum threshold
                optimal_fold_increase *= 1/self.DECREASE_RATIO
                break

            # Save fold increase probability distribution as an object attribute for future reference
            self.fold_increase_pds.append(fold_increase_pd)

            # End while loop once minimum fold increase is simulated
            if optimal_fold_increase == min_fold_increase:
                break

            # Decrease optimal fold increase by the decrease ratio if minimum threshold is still respected (which means
            # that medium usage can be reduced even lower)
            optimal_fold_increase *= self.DECREASE_RATIO

            # Limit optimal fold increase to minimum fold increase
            if optimal_fold_increase < min_fold_increase:
                optimal_fold_increase = min_fold_increase
            

        # Determine optimal medium volume of each cycle
        cycle_medium_volumes = []
        for cycle in range(1, required_cycles+1):
            cycle_medium_volumes.append(self.bioreactors["Min Volume"].min() * optimal_fold_increase**(cycle-1))
        cycle_medium_volumes = np.array(cycle_medium_volumes)

        # Determine bioreactors required for each expansion cycle
        required_bioreactors = pd.DataFrame(0, index=range(1, required_cycles+1),
                                            columns=self.bioreactors.index.tolist())

        for cycle, volume in enumerate(cycle_medium_volumes):
            # Create an auxiliary table with the ratio between the cycle's optimal medium volume
            # and the minimum and maximum volumes of the available bioreactor types
            aux_table = volume / self.bioreactors[["Min Volume", "Max Volume"]]

            # Calculate the max number of bioreactors of a given type which can be properly filled
            aux_table["Min Volume"] = aux_table["Min Volume"].apply(np.floor)

            # Calculate the minimum number of bioreactors of a given type which can contain the medium volume
            aux_table["Max Volume"] = aux_table["Max Volume"].apply(np.ceil)

            # Select bioreactor types where at least 1 bioreactor can be properly filled
            aux_table = aux_table[aux_table["Min Volume"] > 0]

            # Select the bioreactor type where the least number of bioreactors must be used
            aux_table = aux_table[aux_table["Max Volume"] == aux_table["Max Volume"].min()]

            # Save bioreactor required by the current cycle in the respective table
            required_bioreactors.at[cycle+1, aux_table.index[0]] = aux_table["Max Volume"]

            # Calculate minimum volume required to use the required bioreactors
            min_volume = self.bioreactors.loc[aux_table.index[0], "Min Volume"].item() * aux_table["Max Volume"].item()

            # Check if the medium volume is less than the minimum volume. This is done to ensure that if the medium
            # volume does not allow for the proper filling of the required bioreactors an appropriate volume is
            # selected instead. Update optimal medium volumes accordingly
            if volume < min_volume:
                # Correct required fold increase and update each cycle's medium volume accordingly
                corrected_fold_increase = (min_volume / cycle_medium_volumes[0])**(1/cycle)
                for index in range(1, cycle+1):
                    cycle_medium_volumes[index] = cycle_medium_volumes[index-1] * corrected_fold_increase

        return cycle_medium_volumes, required_bioreactors
        

    # Function for determining bioreactor expansion cost
    def determine_bioreactor_expansion_cost(self, db_data, d_facility_cost, d_labor_cost):
        # <>--------------- Consumables Costs ---------------<>
        # Calculate bioreactor cost based on number of bioreactors used        
        N_bioreactors = self.bioreactor_workflow[1].sum()

        bioreactor_use_costs = N_bioreactors * self.bioreactors["Use Cost"]

        T_bioreactor_use_cost = bioreactor_use_costs.sum()

        # <>----------------- Reagent Costs -----------------<>
        # Calculate total dissociation enzyme volume along with corresponding cost
        T_diss_enz_volume = N_bioreactors * self.bioreactors["Max Volume"] * self.DISS_ENZ_VOL_RATIO
        T_diss_enz_cost = T_diss_enz_volume.sum() * db_data["Reagents"].loc[self.dissociation_enzyme].item()

        # Calculate total ROCK Inhibitor cost (ROCKi is added to dissociation enzyme as well as medium)
        T_ROCKi_volume = T_diss_enz_volume.sum() + self.bioreactor_workflow[0].sum()
        T_ROCKi_cost = T_ROCKi_volume * db_data["Reagents"].loc[ROCKI].item()

        # Calculate total supplements cost (!!! when only added at bioreactor inoculation)
        T_supplements_cost = 0
        if self.expansion_simulation["Supplements"].item() != "n/a":
            for supplement in self.expansion_simulation["Supplements"].item().split(";"):
                T_supplements_cost += self.bioreactor_workflow[0].sum() * db_data["Reagents"].loc[supplement].item()
        
        # Calculate total medium volume along with corresponding cost
        T_medium_volume = self.bioreactor_workflow[0].sum() * self.expansion_simulation["Bioreactor Volumes Spent"].item()
        self.T_medium_cost = T_medium_volume * db_data["Reagents"].loc[self.culture_medium].item()

        # <>-------- Quality Control Costs (Subcategory of Reagent Costs)--------<>
        # Calculate flow cytometry cost
        flow_cytometry_cost = db_data["Quality Controls"].loc["Intracellular FC"].item() * (len(INT_FC_ANTIBODIES)+1)
        flow_cytometry_cost += db_data["Quality Controls"].loc["Surface FC"].item() * (len(SUR_FC_ANTIBODIES)+1)
        for antibody in INT_FC_ANTIBODIES:
            flow_cytometry_cost += db_data["Antibodies"].loc[antibody].item()
        for antibody in SUR_FC_ANTIBODIES:
            flow_cytometry_cost += db_data["Antibodies"].loc[antibody].item()

        # Calculate trilineage differentiation cost (cost is calculated for 6 wells of a 12-well plate, 2 wells for
        # each germ layer, and the cost of RT-PCR to analyze differentiation outcome is included)
        diff_planar_platform_data = db_data["2D Platforms"].loc["12-well"]
        diff_total_surfaces = 6

        diff_platform_cost = diff_planar_platform_data["Cost"]

        diff_coating_volume = diff_planar_platform_data["Coating Volume"].item() * diff_total_surfaces
        diff_coating_cost = diff_coating_volume * db_data["Reagents"].loc[self.coating_substrate].item()

        diff_hiPSC_medium_volume = diff_planar_platform_data["Culture Volume"] * diff_total_surfaces
        diff_hiPSC_medium_cost = diff_hiPSC_medium_volume * db_data["Reagents"].loc[self.culture_medium].item()
        diff_ROCKi_cost = diff_hiPSC_medium_volume * db_data["Reagents"].loc[ROCKI].item()

        diff_kit_medium_cost = db_data["Quality Controls"].loc["Trilineage Differentiation"].item()

        diff_RT_PCR_cost = db_data["Quality Controls"].loc["RT-PCR"].item()

        trilineage_differentiation_cost = (diff_platform_cost + diff_coating_cost + diff_hiPSC_medium_cost
                                           + diff_ROCKi_cost + diff_kit_medium_cost + diff_RT_PCR_cost)

        # Calculate immunocytochemistry cost
        immuno_cost = (db_data["Quality Controls"].loc["Intracellular Immunocytochemistry"].item()
                       * (len(INT_IMMUNO_ANTIBODIES)+1))
        immuno_cost += (db_data["Quality Controls"].loc["Surface Immunocytochemistry"].item()
                       * (len(SUR_IMMUNO_ANTIBODIES)+1))
        for antibody in INT_IMMUNO_ANTIBODIES:
            immuno_cost += db_data["Antibodies"].loc[antibody].item()
        for antibody in SUR_IMMUNO_ANTIBODIES:
            immuno_cost += db_data["Antibodies"].loc[antibody].item()        

        # Calculate RT-PCR cost
        RT_PCR_cost = db_data["Quality Controls"].loc["RT-PCR"].item()

        # Calculate karyotyping cost
        karyotyping_cost = db_data["Quality Controls"].loc["Karyotyping"].item()

        # Calculate PCR genomic screening cost
        genetic_analysis_cost = db_data["Quality Controls"].loc["PCR Genomic Screening"].item()

        # Calculate cost of each type of quality control
        ini_qual_cost = flow_cytometry_cost + trilineage_differentiation_cost
        int_qual_cost = flow_cytometry_cost * (len(self.bioreactor_workflow[0])-1)
        fin_qual_cost = (flow_cytometry_cost + trilineage_differentiation_cost + immuno_cost
                         + RT_PCR_cost + karyotyping_cost + genetic_analysis_cost)
        
        # Calculate total quality control cost
        T_qual_ctrl_cost = ini_qual_cost + int_qual_cost + fin_qual_cost

        # <>----------------- Facility Costs ----------------<>
        # Calculate daily bioreactor depreciation
        bioreactor_acquisition_costs = N_bioreactors * self.bioreactors["Acquisition Cost"]

        T_bioreactor_acquisition_cost = bioreactor_acquisition_costs.sum()

        self.d_bioreactor_depreciation = (T_bioreactor_acquisition_cost /
                                          (db_data["Facility Specifications"].loc["Equipment Lifespan"].item() * YEAR_TO_DAYS))
        
        # Calculate total bioreactor energy cost
        bioreactor_energy_consumptions = (N_bioreactors * self.bioreactors["Energy Consumption"]
                                          * self.expansion_simulation["Bioreactor Culture Time"].item())

        T_bioreactor_energy_consumption = bioreactor_energy_consumptions.sum()

        T_bioreactor_energy_cost = (T_bioreactor_energy_consumption
                                    * db_data["Facility Specifications"].loc["Energy Cost"].item())

        # <>--------------- Costs by Category ---------------<>
        # Calculate total category costs of bioreactor expansion
        self.T_consumables_cost = T_bioreactor_use_cost
        self.T_reagents_cost = T_diss_enz_cost + T_ROCKi_cost + self.T_medium_cost + T_qual_ctrl_cost
        self.T_facility_cost = d_facility_cost * self.duration + T_bioreactor_energy_cost
        self.T_labor_cost = d_labor_cost * self.duration

        # Calculate overall bioreactor expansion cost
        self.overall_cost = self.T_consumables_cost + self.T_reagents_cost + self.T_facility_cost + self.T_labor_cost


# Define composite class for bioprocess simulation and computation of costs
class Bioprocess():
    # Initializer of class object
    def __init__(self, db_data, bio_params):
        # <>------------------- Main Body -------------------<>
        # Determine daily facility cost
        self.d_facility_cost = self.determine_daily_facility_cost(db_data)

        # Determine daily labor cost
        self.d_labor_cost = self.determine_daily_labor_cost(db_data)

        # Create instance of Planar Expansion Class
        self.planar_expansion = PlanarExpansion(db_data, bio_params, self.d_facility_cost, self.d_labor_cost)

        # Create instance of Bioreactor Expansion Class
        self.bioreactor_expansion = BioreactorExpansion(db_data, bio_params, self.d_facility_cost, self.d_labor_cost)

        # Adjust the facility cost of both phases by taking into account bioreactor depreciation 
        # (this can only be done after determining the optimal workflow)
        self.planar_expansion.T_facility_cost += (self.bioreactor_expansion.d_bioreactor_depreciation
                                                  * self.planar_expansion.duration)
        self.planar_expansion.overall_cost += (self.bioreactor_expansion.d_bioreactor_depreciation
                                                  * self.planar_expansion.duration)

        self.bioreactor_expansion.T_facility_cost += (self.bioreactor_expansion.d_bioreactor_depreciation
                                                      * self.bioreactor_expansion.duration)
        self.bioreactor_expansion.overall_cost += (self.bioreactor_expansion.d_bioreactor_depreciation
                                                      * self.bioreactor_expansion.duration)

        # <>--------------- Costs by Category ---------------<>
        # Calculate total category costs of entire bioprocess
        self.bioprocess_consumables_cost = self.planar_expansion.T_consumables_cost + self.bioreactor_expansion.T_consumables_cost
        self.bioprocess_reagents_cost = self.planar_expansion.T_reagents_cost + self.bioreactor_expansion.T_reagents_cost
        self.bioprocess_facility_cost = self.planar_expansion.T_facility_cost + self.bioreactor_expansion.T_facility_cost
        self.bioprocess_labor_cost = self.planar_expansion.T_labor_cost + self.bioreactor_expansion.T_labor_cost
        
        # Calculate total medium cost of bioprocess
        self.bioprocess_medium_cost = self.planar_expansion.T_medium_cost + self.bioreactor_expansion.T_medium_cost

        # Calculate overall bioprocess cost
        self.bioprocess_overall_cost = self.planar_expansion.overall_cost + self.bioreactor_expansion.overall_cost


    # Function for determining daily facility cost
    def determine_daily_facility_cost(self, db_data):        
        # Calculate facility construction cost
        standard_rooms_cost = (db_data["Construction Costs"].loc["Standard Rooms", "Area"]
                               * db_data["Construction Costs"].loc["Standard Rooms", "Cost"])
        clean_rooms_cost = (db_data["Construction Costs"].loc["Clean Rooms", "Area"]
                            * db_data["Construction Costs"].loc["Clean Rooms", "Cost"]) 
        facility_construction_cost = standard_rooms_cost.item() + clean_rooms_cost.item()
        
        # Calculate daily facility depreciation (CC / facility lifespan)
        d_facility_depreciation = (facility_construction_cost
                                   / (db_data["Facility Specifications"].loc["Facility Lifespan"].item() * YEAR_TO_DAYS))

        # Calculate equipment acquisition cost
        equipment_cost = (db_data["Equipment"].loc[:, "Amount"] * db_data["Equipment"].loc[:, "Acquisition Cost"]).sum()

        # Calculate daily equipment depreciation (excluding bioreactors)
        d_equipment_depreciation = (equipment_cost
                                    / (db_data["Facility Specifications"].loc["Equipment Lifespan"].item() * YEAR_TO_DAYS))
        
        # Calculate total daily depreciation cost
        d_depreciation_cost = d_facility_depreciation + d_equipment_depreciation

        # Calculate daily facility operating cost
        d_facility_operating_cost = db_data["Operating Costs"].sum().item() / YEAR_TO_DAYS

        # Calculate daily equipment energy cost (excluding bioreactors)
        d_energy_cost = (db_data["Equipment"].loc[:, "Amount"] * db_data["Equipment"].loc[:, "Energy Consumption"]
                         * db_data["Equipment"].loc[:, "Use Factor"]
                         * db_data["Facility Specifications"].loc["Energy Cost"].item()).sum()
        
        # Calculate daily facility cost (taking into account number of parallel bioprocesses)
        d_facility_cost = d_depreciation_cost + d_facility_operating_cost + d_energy_cost
        d_facility_cost /= db_data["Facility Specifications"].loc["Parallel Processes"].item()

        return d_facility_cost


    # Function for determining daily labor cost
    def determine_daily_labor_cost(self, db_data):        
        # Calculate daily labor cost (taking into account number of parallel bioprocesse)
        d_labor_cost = ((db_data["Labor Costs"].loc[:, "Number"] * db_data["Labor Costs"].loc[:, "Salary"]).sum()
                        / YEAR_TO_DAYS / db_data["Facility Specifications"].loc["Parallel Processes"].item())

        return d_labor_cost