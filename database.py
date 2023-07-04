import math
import numpy as np
import pandas as pd
from numpy.random import *

# Global Constants
DOLAR_TO_EURO = 0.8


class DataTable():

    def __init__(self):

        # REAGENTS
        # Reagents for Expansion
        reagents = ['Accutase', 'DMEM_F12', 'DPBS', 'DS', 'EDTA', 'Matrigel', 'Y_27632']

        self.accutase_cost = 468 #[€/L] (Sigma-Aldrich: SCR005)

        self.DMEM_F12_cost = 54.18 #[€/L] (Thermo Fisher Scientific: 31330038)

        self.DPBS_cost = 40.13 #[€/L] (Thermo Fisher Scientific: 14190136)

        self.DS_cost = 2.79 #[€/L] --> concentration of 100 mg/L in culture medium and 0.0248 €/mg
                            # (Sigma-Aldrich: RES2029D-A7)

        self.EDTA_cost = 0.33 #[€/L] --> concentration of 0.5 mM in DPBS and 0.5 M is 305 €/L
                              # (Thermo Fisher Scientific: 15575020)

        self.matrigel_cost = 708.94 #[€/L] --> dilution of 1:100 in DMEM/F12 and 41.3078 €/mL (Corning: 354263)

        self.Y_27632_cost = 139.01 #[€/L] --> concentration of 10 uM, molar mass of 320.3e3 mg/mol
                                   # and 41.2 €/mg (STEMCELL Technologies: 72304)

        self.reagents_costs = pd.DataFrame({'Cost (€/L)': [self.accutase_cost, self.DMEM_F12_cost, self.DPBS_cost,
                  self.DS_cost, self.EDTA_cost, self.matrigel_cost, self.Y_27632_cost]}, index=reagents) #[€]

        # Reagents for Quality Control
        # Flow Cytometry
        self.PBS_cost = 39.3 #[€/L] (Sigma-Aldrich: 806552-1L)

        self.fc_PFA_cost = 40.56 #[€/L] --> 2% solution, which is 20 g/L, and 0.0628 €/g
                                # (Sigma-Aldrich: 158127-500G)

        self.NGS1_cost = 262.3 #[€/L] --> 1% solution in PBS, which is 10 mL/L, and 22.3 €/mL
                               # (Sigma-Aldrich: 566380-10ML)

        self.NGS3_cost = 708.3 #[€/L] --> 3% solution in PBS, which is 30 mL/L, and 22.3 €/mL
                               # (Sigma-Aldrich: 566380-10ML)

        self.saponin_cost = 53.44 #[€/L] --> 1% solution in PBS, which is 10 g/L, and 1.414 €/g
                                  # (Sigma-Aldrich: 47036-50G-F)

        self.fc_OCT4_antibody_cost = 11500 #[€/L] --> 1:300 dilution and 3450 €/mL
                                           # (Thermo Fisher Scientific: MA1-104)

        self.fc_SOX2_antibody_cost = 22600 #[€/L] --> 1:150 dilution and 3390 €/mL
                                           # (Thermo Fisher Scientific: MA1-014)

        self.FBS_cost = 95.38 #[€/L] --> 4% solution in PBS, which is 40 mL/L, and 1.402 €/mL
                              # (Sigma-Aldrich: F2442-50ML)

        self.fc_TRA_1_60_PE_antibody_cost = 38800 #[€/L] --> 1:10 dilution, 388 €/mL
                                                   # (Thermo Fisher Scientific: 12-8863-82)

        # CHANGE!!!!!
        self.fc_SSEA_1_FITC_antibody_cost = 24000 #[€/L] --> 1:20 dilution, 480 €/mL (Thermo Fisher Scientific)

        self.fc_SSEA_1_PE_antibody_cost = 55800 #[€/L] --> 1:10 dilution, 558 €/mL
                                                # (Thermo Fisher Scientific: 12-8813-42)
        # CHANGE!!!!!

        self.fc_secondary_antibody_cost = 1660 #[€/L] --> 1:300 dilution, 498 €/mL
                                               # (Thermo Fisher Scientific: A-11001)

        # qRT-PCR
        self.OCT4_primer_cost = round(120 / 200, 2) #[€/use] --> each primer vial allows for 200 reactions
                                                    # (OriGene: HP206340)

        self.NANOG_primer_cost = round(120 / 200, 2) #[€/use] --> each primer vial allows for 200 reactions
                                                    # (OriGene: HP215086)

        self.SOX1_primer_cost = round(120 / 200, 2) #[€/use] --> each primer vial allows for 200 reactions
                                                    # (OriGene: HP209017)

        self.T_primer_cost = round(120 / 200, 2) #[€/use] --> each primer vial allows for 200 reactions
                                                 # (OriGene: HP206752)

        self.SOX17_primer_cost = round(120 / 200, 2) #[€/use] --> each primer vial allows for 200 reactions
                                                    # (OriGene: HP214460)

        # Immunocytochemistry
        self.ic_PFA_cost = 41.81 #[€/L] --> 4% solution, which is 40 g/L, and 0.0628 €/g
                                 # (Sigma-Aldrich: 158127-500G)

        self.triton_cost = 39.78 #[€/L] --> 0.1% solution, which is 1 mL/L, and 0.482 €/mL
                                 # (Sigma-Aldrich: X100-100ML)

        self.NGS10_cost = 2269.3 #[€/L] --> 10% solution in PBS, which is 100 mL/L, and 22.3 €/mL
                                 #(Sigma-Aldrich: 566380-10ML)

        self.NGS5_cost = 1154.3 #[€/L] --> 5% solution in PBS, which is 50 mL/L, and 22.3 €/mL
                                # (Sigma-Aldrich: 566380-10ML)

        self.DAPI_cost = 55.95 #[€/L] --> 1.5 ug/mL and 11.1 €/mg (Sigma-Aldrich: D9542-1MG)

        self.BSA_cost = 267 #[€/L] --> 3% solution, which is 30 mL/L, and 8.9 €/mL (Sigma-Aldrich: A9576-50ML)

        self.ic_OCT4_antibody_cost = 6900 #[€/L] --> 1:500 dilution and 3450 €/mL
                                           # (Thermo Fisher Scientific: MA1-104)

        self.ic_SOX2_antibody_cost = 16950 #[€/L] --> 1:200 dilution, 3390 €/mL
                                           # (Thermo Fisher Scientific: MA1-014)

        self.ic_TRA_1_60_antibody_cost = 8592.59 #[€/L] --> 1:135 dilution, 1160 €/mL
                                                 # (Thermo Fisher Scientific: 14-8863-82)

        self.ic_SSEA_4_antibody_cost = 14370.37 #[€/L] --> 1:135 dilution, 1940 €/mL
                                                # (Thermo Fisher Scientific: 41-4000)

        self.ic_secondary_antibody_cost = 996 #[€/L] --> 1:500 dilution, 498 €/mL (Thermo Fisher Scientific)

        # Karyotyping (G-Banding)
        self.colcemid_cost = 212.45 #[€/L] --> 0.05 mL/mL and 4.249 €/mL (Thermo Fisher Scientific: 15212012)

        self.KCl_cost = 155.63 #[€/L] (Thermo Fisher Scientific: 10575090)

        self.gurr_buffer_cost = 1.52 #[€/tablet] (Thermo Fisher Scientific: 10582013)

        self.trypsin_cost = 184.4 #[€/L] (Thermo Fisher Scientific: 15090046)

        self.giemsa_stain_cost = 458 #[€/L] (Thermo Fisher Scientific: 10092013)

        # Differentiation
        self.RPMI_cost = 36.06 #[€/L] (Thermo Fisher Scientific: 21875034)

        self.B27_cost = 384 #[€/L] --> 2% solution, which is 20 mL/L, and 19.2 €/mL
                            # (Thermo Fisher Scientific: A1486701)

        self.B27_no_insulin_cost = 398 #[€/L] --> 2% solution, which is 20 mL/L, and 19.9 €/mL
                                       # (Thermo Fisher Scientific: A3695201)

        self.penicillin_streptomycin_cost = 0.80 #[€/L] --> 0.5% solution, which is 5 mL/L, and 0.16 €/mL
                                                 # (Thermo Fisher Scientific)

        self.CHIR99021_cost = 42.52 #[€/L] --> 6 uM, which is 2.79 mg/L, and 15.24 €/mg
                                    # (Sigma-Aldrich: SML1046-25MG)

        self.IWP4_cost = 71.12 #[€/L] --> 5 uM, which is 2.48 mg/L, and 28.64 €/mg (Sigma-Aldrich: SML1114-25MG)

        self.methanol_cost = 30.33 #[€/L] (Sigma-Aldrich: 34860-1L-R)

        self.BSA_05_cost = 83.8 #[€/L] --> 0.5% solution, which is 5 mL/L, and 8.9 €/mL (Sigma-Aldrich: A9576-50ML)

        self.cTnT_antibody_cost = 15600 #[€/L] --> 1:250 dilution, 3900 €/mL (Thermo Fisher Scientific: MA5-17192)

        self.secondary_antibody_cost = 498 #[€/L] --> 1:1000 dilution, 498 €/mL (Thermo Fisher Scientific: A-11001)

        # CONSUMABLES
        # Planar Expansion Platforms (peps)
        peps = ['12-well', '6-well', 'T-25', 'T-75']
            # the values relevant to the 6-well plates are indicated per well

        peps_unit_costs = np.array([0.21, 0.40, 0.90, 1.38]) #[€/unit] (Thermo Fisher Scientific: 150628, 140675)

        peps_matrigel_volumes = np.array([0.5, 1, 3.5, 10.5]) #[mL]
            # the matrigel coating solution is composed of matrigel diluted 1:100 in DMEM/F12

        peps_matrigel_costs = np.round(peps_matrigel_volumes * 1e-3 * (self.DMEM_F12_cost + self.matrigel_cost), 2)
            # [0.23, 0.47, 1.63, 4.89] €

        peps_total_costs = peps_unit_costs + peps_matrigel_costs
            # [0.44, 0.83, 2.53, 6.27] €

        peps_areas = np.array([3.5, 9.6, 25, 75]) #[cm2]

        peps_cells_seeding = np.array([0.1e6, 0.3e6, 0.7e6, 2.1e6]) #[cells]

        peps_cells_confluency = 4 * peps_cells_seeding #[cells]
            # [0.4e6, 1.2e6, 2.8e6, 8.4e6] cells (the fold increase at confluency is approximately 4)

        peps_medium_volumes = np.array([1, 2, 5, 15]) #[mL]

        peps_DPBS_volumes = np.array([1, 2, 5, 15]) #[mL]

        peps_passaging_solution_volumes = np.array([0.5, 1, 3, 8]) #[mL]
            # the passaging solution is composed of 0.5 mM EDTA dissolved in DPBS

        self.peps_costs = pd.DataFrame({'Cost (€)': peps_unit_costs, 'Matrigel Vol (mL)': peps_matrigel_volumes,
               'Matrigel Cost (€)': peps_matrigel_costs, 'Total Cost (€)': peps_total_costs}, index=peps)

        self.peps_parameters = pd.DataFrame({'Area (cm2)': peps_areas, 'Medium Vol (mL)': peps_medium_volumes,
               'DPBS Vol (mL)': peps_medium_volumes, 'PS Vol (mL)': peps_passaging_solution_volumes}, index=peps)

        self.peps_cells = pd.DataFrame({'Cells Seeding': peps_cells_seeding, 'Cells Confluency':
                                        peps_cells_confluency}, index=peps)

        # <6-well> Values from gibco episomal hiPSC line (Thermo Fisher Scientific - MAN0009562)
        self.recovery_medium_per_million = 11 #[mL/1e6 cells]
            # medium used for thawing of frozen hiPSCs

        self.well_culture_duration = 4 #[days]

        self.maximum_passaging_ratio = 4

        self.well_parameters = pd.DataFrame([[self.recovery_medium_per_million, self.well_culture_duration,
               self.maximum_passaging_ratio]], columns=['Recovery Medium (mL/1e6 cells)', 'Culture Duration (d)',
               'Max Passaging Ratio'])

        # Vertical-Wheel Bioreactors (vwbs)
        self.vwbs = ['PBS 0.1MAG', 'PBS 0.5MAG', 'PBS 3MAG']

        vwbs_vessel_costs = np.array([261.75, 327.75, 1248]) #[€/vessel] (PBS Biotech)

        vwbs_min_working_volumes = np.array([60, 300, 1800]) #[mL]

        vwbs_max_working_volumes = np.array([100, 500, 3000]) #[mL]

        self.vwbs_parameters = pd.DataFrame({'Cost (€)': vwbs_vessel_costs, 'Min Vol (mL)':
                         vwbs_min_working_volumes, 'Max Vol (mL)': vwbs_max_working_volumes}, index=self.vwbs)

        vwbs_base_units = ['PBS Mini', 'PBS 3MAG']

        vwbs_base_units_costs = np.array([2473, 65550]) #[€/base] (PBS Biotech)

        self.vwbs_base_units = pd.DataFrame({'Cost (€)': vwbs_base_units_costs}, index=vwbs_base_units)

        # Kits
        self.RNA_purification_kit_cost = round(307 / 50, 2) #[€/use] --> each kit has 50 preps
                                                            # (STEMCELL Technologies: 79040)

        self.cDNA_transcriptase_kit_cost = round(463 / 100, 2) #[€/use] --> each kit allows for 100 reactions
                                                               # (STEMCELL Technologies: 79004)

        self.qPCR_master_mix_kit_cost = round(101 / 200, 2) #[€/use] --> each kit allows for 200 reactions
                                                            # (NZYTech: MB22201)

        self.DNA_purification_kit_cost = round(862 / 250, 2) #[€/use] --> each kit has 250 preps
                                                             # (STEMCELL Technologies: 79020)

        self.hPSC_genetic_kit_cost = round(534 / 20, 2) #[€/use] --> each kit allows the ananlysis of 20 samples
                                                        # (STEMCELL Technologies: 07550)

        # CULTURE MEDIA
        mTeSR1_cost = 654 #[€/L] (STEMCELL Technologies: 85857)

        E8_cost = 584 #[€/L] (Thermo Fisher Scientific: A1517001)

        in_house_E8_cost = 112.29 #[€/L]

        B8_cost = 9.63 #[€/L] --> 9.83 $/L (in-house) !!! Reagent cost only !!!

        culture_media = ['mTeSR1', 'E8', 'ihE8', 'B8']

        self.culture_media_costs = pd.DataFrame({'Cost (€/L)': [mTeSR1_cost, E8_cost, in_house_E8_cost, B8_cost]},
                                                 index=culture_media)

        # LABOR
        self.salary_per_hour = 12 #[€/h] --> online research on the salary of a laboratory technician in Portugal,
                                  # in accordance with Bandeiras et al.

        self.daily_worktime = 8 #[h/day]

        self.weekly_worktime = 40 #[h/week] --> 8 hours of work per day, 5 days per week (Bandeiras Thesis)

        self.yearly_worktime = 48 #[weeks/year] --> 1 month vacation (Bandeiras Thesis)

        # FACILITY
        # Building Costs (Bandeiras Thesis)
        self.cost_clean_room_per_sqmt = 6329 #[€/m2]

        self.cost_generic_room_per_sqmt = 3692 #[€/m2]

        # Depreciation of Facility (Bandeiras Thesis)
        self.facility_depreciation_period = 15 * 365.25 #[days]

        # EQUIPMENT
        equipment = ['Incubator', 'Centrifuge', 'BSC', 'Flow Cytometer', 'Spectrophotometer',
                     'RT-PCR System', 'Microscope', 'Fluorescence Microscope']

        incubator_cost = 9604 #[€] (Thermo Fisher Scientific: 15650667)

        centrifuge_cost = 7163 #[€] (Thermo Fisher Scientific: 15808722)

        bsc_cost = 10865 #[€] (Thermo Fisher Scientific: 10445753)

        flow_cytometer_cost = 94700 #[€] (Thermo Fisher Scientific: 15310677)

        spectrophotometer_cost = 10169.39 #[€] (Thermo Fisher Scientific: 10506405)

        RT_PCR_system_cost = 20560 #[€] (Thermo Fisher Scientific: 15341295)

        microscope_cost = 2471 #[€] (Thermo Fisher Scientific: 11712298)

        fluorescence_microscope_cost = 12245 #[€] (Thermo Fisher Scientific: 11960279)

        self.equipment_costs = pd.DataFrame({'Cost (€)': [incubator_cost, centrifuge_cost, bsc_cost,
                                            flow_cytometer_cost, spectrophotometer_cost, RT_PCR_system_cost,
                                            microscope_cost, fluorescence_microscope_cost]}, index=equipment)

        self.equipment_depreciation_period = 5 * 365.25 # [days]


class Facility():

    def __init__(self, dt):

        # Facility Specifications
        self.facility_area = 400 #[m2]

        self.clean_room_ratio = 0.2

        self.parallel_processes = 6
            # this number cannot exceed the amount of vwb Mini or 3MAG base units

        # Daily Operational Costs (Bandeiras Thesis)  !!!CHANGE TO YEARLY COST FOR LESS ROUNDING ERROR!!!
        co2_daily_cost = round(6000 / 365.25, 2) #[€/day]

        gases_daily_cost = round(15600 / 365.25, 2) #[€/day]

        additional_supplies_daily_cost = round(7900 / 365.25, 2) #[€/day]

        requalification_daily_cost = round(65400 / 365.25, 2) #[€/day]

        maintenance_daily_cost = round(52800 / 365.25, 2) #[€/day]

        cleaning_daily_cost = round(28000 / 365.25, 2) #[€/day]

        garment_daily_cost = round(2000 / 365.25, 2) #[€/day]

        self.operational_daily_costs = round(co2_daily_cost + gases_daily_cost + additional_supplies_daily_cost
               + requalification_daily_cost + maintenance_daily_cost + cleaning_daily_cost + garment_daily_cost, 2)

        # Total Facility Costs
        self.total_facility_costs = round((dt.cost_clean_room_per_sqmt * self.clean_room_ratio
                     + dt.cost_generic_room_per_sqmt * (1 - self.clean_room_ratio)) * self.facility_area, 2) #[€]

        self.daily_facility_costs = round(self.total_facility_costs / dt.facility_depreciation_period
                                     + self.operational_daily_costs, 2)
            # daily facility costs considering depreciation of the facility

        # Available Equipment and Associated Costs
        self.incubator_amount = 8

        self.centrifuge_amount = 4

        self.bsc_amount = 4

        self.flow_cytometer_amount = 1

        self.spectrophotometer_amount = 1

        self.RT_PCR_system_amount = 1

        self.microscope_amount = 2

        self.fluorescence_microscope_amount = 1

        self.base_units_Mini = 18

        self.base_units_3MAG = 6

        self.total_equipment_costs = round(self.incubator_amount * dt.equipment_costs.at['Incubator', 'Cost (€)']
                + self.centrifuge_amount * dt.equipment_costs.at['Centrifuge', 'Cost (€)'] + self.bsc_amount
                * dt.equipment_costs.at['BSC', 'Cost (€)'] + self.flow_cytometer_amount
                * dt.equipment_costs.at['Flow Cytometer', 'Cost (€)'] + self.spectrophotometer_amount
                * dt.equipment_costs.at['Spectrophotometer', 'Cost (€)'] + self.RT_PCR_system_amount
                * dt.equipment_costs.at['RT-PCR System', 'Cost (€)'] + self.microscope_amount
                * dt.equipment_costs.at['Microscope', 'Cost (€)'] + self.fluorescence_microscope_amount
                * dt.equipment_costs.at['Fluorescence Microscope', 'Cost (€)'] + self.base_units_Mini
                * dt.vwbs_base_units.at['PBS Mini', 'Cost (€)'] + self.base_units_3MAG
                * dt.vwbs_base_units.at['PBS 3MAG', 'Cost (€)'], 2)

        # daily equipment costs considering depreciation of the facility
        self.daily_equipment_costs = round(self.total_equipment_costs / dt.equipment_depreciation_period, 2)

        self.overall_daily_facility_costs = self.daily_facility_costs + self.daily_equipment_costs

        # Labor
        self.employee_amount = 4

        self.daily_employee_costs = self.employee_amount * dt.daily_worktime * dt.salary_per_hour


class CultureConditions():

    def __init__(self, **conditions):

        # Default Culture Conditions
        self.initial_cell_number = 1e6
        self.target_cell_number = 2e9
        self.culture_medium = 'mTeSR1'
        self.initial_volume = 60 #[mL]

        # Default Simulation Source
        self.simulation_source = 'Nogueira'
        self.DS_supplementation = False

        # Save Custom Culture Conditions and Simulation Source
        if len(conditions) == 6:
            self.initial_cell_number = conditions['initnum']
            self.target_cell_number = conditions['targnum']
            self.culture_medium = conditions['medium']
            self.initial_volume = conditions['initvol']
            self.simulation_source = conditions['sim']
            self.DS_supplementation = conditions['DS']

        # Set Values Dependent on Simulation Source
        # values from Rodrigues et al. (2019)
        if self.simulation_source == 'Nogueira':
            self.seeding_density = 2.5e5 #[cells/mL]

            if self.DS_supplementation:
                self.fold_increase_average = 9.3
                self.fold_increase_sem = 0.6
                self.fold_increase_sample_size = 3
                self.fold_increase_std = self.fold_increase_sem * math.sqrt(self.fold_increase_sample_size)
                self.expansion_duration = 5 #[days]
            else:
                self.fold_increase_average = 4.8
                self.fold_increase_sem = 0.5
                self.fold_increase_sample_size = 3
                self.fold_increase_std = self.fold_increase_sem * math.sqrt(self.fold_increase_sample_size)
                self.expansion_duration = 7 #[days]

        # values from Borys et al. (2021)
        if self.simulation_source == 'Borys':
            self.seeding_density = 2e4 #[cells/mL]

            self.fold_increase_average = 32.3
            self.fold_increase_sem = 3.2
            self.fold_increase_sample_size = 4
            self.fold_increase_std = self.fold_increase_sem * math.sqrt(self.fold_increase_sample_size)
            self.expansion_duration = 6 #[days]

            self.DS_supplementation = False
                # not tested in the study of Borys et al. (2021)

        self.recovery_efficiency_average = 0.952
        self.recovery_efficiency_sem = 0.040
        self.recovery_efficiency_sample_size = 4
        self.recovery_efficiency_std = (self.recovery_efficiency_sem
                                         * math.sqrt(self.recovery_efficiency_sample_size))

        # Minimum Success Threshold
        self.threshold = 0.95 # 95%

        # Quality Control Variables
        self.final_quality_control_duration = 3 #[days]

        self.qual_ctrl_detection_prob = 0.2 # 20%
        self.bioprocess_failure_prob = 0.1 # 10%

        # Auxiliary Variables
        # determine initial well amount according to the initial cell number
        self.initial_well_amount = math.ceil(self.initial_cell_number / 1e6)

        # determine necessary cell for inoculation of the initial volume
        self.inoculation_cell_number = math.ceil(self.initial_volume * self.seeding_density)


    # Define Feeding Strategy According to Simulation Source and Medium
    def Nogueira_feeding_strategy(self, working_volume):
        return working_volume * 0.8 * (self.expansion_duration - 2)
            # daily exchange of 80% of medium starting 48h after inoculation


    def Nogueira_B8_feeding_strategy(self, working_volume):
        return working_volume * 0.8 * math.ceil((self.expansion_duration - 2) / 2)
            # exchange of 80% of medium starting 48h after inoculation every other day


    def Borys_feeding_strategy(self, working_volume):
        return working_volume * 0.5
            # exchange of 50% of medium on the fourth day after inoculation


    # Define Bioprocess Functions
    # define function for stochastic simulation of fold increase following a normal distribution
    def fold_increase(self):
        avg = self.fold_increase_average
        std = self.fold_increase_std
        # return triangular(avg - std * math.sqrt(6), avg, avg + std * math.sqrt(6))
        return normal(avg, std)


    # define function for stochastic simulation of recovery efficiency following a beta distribution
    # !!!!Calculate alpha and beta outside!!!!
    def recovery_efficiency(self):
        avg = self.recovery_efficiency_average
        std = self.recovery_efficiency_std / 3

        # beta distribution parameters
        precision = ((avg * (1 - avg)) / std**2) - 1

        alpha_parameter = avg * precision
        beta_parameter = (1 - avg) * precision

        # return triangular(avg - std * math.sqrt(6) * 0.33, avg, avg + std * math.sqrt(6) * 0.33)
        # return normal(avg, std * 0.33)
        return beta(alpha_parameter, beta_parameter)

    #define function for stochastic simulation of 6-well cell confluency following a normal distribution
    def well_final_cell_number(self, dt):
        avg = dt.peps_cells.at['6-well', 'Cells Confluency']
        std = avg / 12
        # return triangular(avg - limit, avg, avg + limit)
        return normal(avg, std)
