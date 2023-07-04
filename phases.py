import math
import numpy as np
import pandas as pd
from database import *
from logistics import *


class PreInoculation():

    def __init__(self, dt, ft, cc, necessary_well_amount):

        # set relevant culture medium cost
        self.culture_medium_cost = dt.culture_media_costs.at[cc.culture_medium, 'Cost (€/L)']

        # Recovery of Frozen hiPSCs
        # determine medium volume necessary for recovery of the initial cell number
        self.recovery_medium_volume = dt.recovery_medium_per_million * cc.initial_well_amount

        # Passaging of hiPSCs
        # determine the ideal passaging work flow to achieve the necessary well amount
        (self.passage_number, self.passaging_workflow) = create_2D_workflow(dt, cc.initial_well_amount,
                                                                            necessary_well_amount)

        # determine the total DPBS volume used for washing wells before passaging of hiPSCs
        if self.passage_number == 0:
            self.passaged_well_amount = 0
        else:
            self.passaged_well_amount = (cc.initial_well_amount
                                           + sum(self.passaging_workflow['Passage Well Amount'][:-1]))

        self.pre_DPBS_volume = dt.peps_parameters.at['6-well', 'DPBS Vol (mL)'] * 2 * self.passaged_well_amount

        # determine the total passaging solution (ps) used for passaging of hiPSCs
        self.pre_ps_volume = dt.peps_parameters.at['6-well', 'PS Vol (mL)'] * self.passaged_well_amount

        # Duration [days]
        self.pre_duration = (self.passage_number + 1) * dt.well_culture_duration

        # Worktime [minutes]
        self.total_well_amount = sum(self.passaging_workflow['Passage Well Amount'])

        self.pre_well_coating_worktime = 10 + math.ceil(self.total_well_amount / 6) * 5

        self.pre_recovery_worktime = 20 + 20 * cc.initial_well_amount

        self.pre_medium_exchange_worktime = ((5 * (self.passage_number + 1) + self.total_well_amount * 1) *
                                          dt.well_culture_duration)

        self.pre_passaging_worktime = 20 * self.passage_number + self.passaged_well_amount * 3

        self.pre_total_worktime = (self.pre_well_coating_worktime + self.pre_recovery_worktime +
                                   self.pre_medium_exchange_worktime + self.pre_passaging_worktime)

        # convert the worktime to hours and round it to half hours
        self.pre_total_worktime = (math.ceil(self.pre_total_worktime / 30) * 30) / 60

        # Calculation of Costs
        self.total_well_cost = round(dt.peps_costs.at['6-well', 'Total Cost (€)'] * self.total_well_amount, 2)

        self.recovery_medium_cost = round(self.culture_medium_cost * self.recovery_medium_volume * 1e-3, 2)

        self.pre_DPBS_cost = round(dt.DPBS_cost * self.pre_DPBS_volume * 1e-3, 2)

        self.pre_ps_cost = round((dt.EDTA_cost + dt.DPBS_cost) * self.pre_ps_volume * 1e-3, 2)

        self.pre_daily_medium_volume = (dt.peps_parameters.at['6-well', 'Medium Vol (mL)'] * self.total_well_amount
                                          * dt.well_culture_duration)

        self.pre_daily_medium_cost = round(self.culture_medium_cost * self.pre_daily_medium_volume * 1e-3, 2)

        self.pre_total_medium_cost = round(self.recovery_medium_cost + self.pre_daily_medium_cost, 2)

        self.pre_total_non_medium_cost = round(self.total_well_cost + self.pre_DPBS_cost + self.pre_ps_cost, 2)

        # Organization of Costs into Categories
        self.pre_consumables_cost = self.total_well_cost

        self.pre_reagent_cost = round(self.pre_total_medium_cost + self.pre_DPBS_cost + self.pre_ps_cost, 2)

        self.pre_medium_cost = self.pre_total_medium_cost

        self.pre_facility_cost = round(ft.overall_daily_facility_costs * self.pre_duration
                                       / ft.parallel_processes, 2)

        self.pre_labor_cost = round(self.pre_duration * ft.daily_employee_costs / ft.parallel_processes, 2)

        self.pre_total_cost = round(self.pre_consumables_cost + self.pre_reagent_cost + self.pre_facility_cost
                                    + self.pre_labor_cost, 2)


class PrimaryInoculation():

    def __init__(self, dt, ft, cc, working_volume, vwb_types):

        # Treatment of hiPSCs Prior to Harvesting for Inoculation
        self.final_well_amount = self.passaging_workflow['Passage Well Amount'][-1]

        self.incubation_medium_volume = (dt.peps_parameters.at['6-well', 'Medium Vol (mL)']
                                           * self.final_well_amount)

        self.ino_DPBS_volume = dt.peps_parameters.at['6-well', 'DPBS Vol (mL)'] * self.final_well_amount

        self.ino_accutase_volume = dt.peps_parameters.at['6-well', 'PS Vol (mL)'] * self.final_well_amount

        # Worktime [minutes]
        self.ino_treatment_worktime = 25 + self.final_well_amount * 3

        self.ino_inoculation_worktime = 60 + 10 * vwb_types['Amount'].sum()

        self.ino_total_worktime = self.ino_treatment_worktime + self.ino_inoculation_worktime

        # convert the worktime to hours and round it to half hours
        self.ino_total_worktime = (math.ceil(self.ino_total_worktime / 30) * 30) / 60

        # Calculation of Costs
        self.incubation_medium_cost = round(self.culture_medium_cost * self.incubation_medium_volume * 1e-3, 2)

        self.incubation_Y_27632_cost = round(dt.Y_27632_cost * self.incubation_medium_volume * 1e-3, 2)

        self.ino_DPBS_cost = round(dt.DPBS_cost * self.ino_DPBS_volume * 1e-3, 2)

        self.ino_accutase_cost = round(dt.accutase_cost * self.ino_accutase_volume * 1e-3, 2)

        self.bioreactors_cost = vwb_types['Associated Costs (€)'].sum()

        self.inoculation_medium_cost = round(self.culture_medium_cost * working_volume * 1e-3, 2)

        self.inoculation_Y_27632_cost = round(dt.Y_27632_cost * working_volume * 1e-3, 2)

        if cc.DS_supplementation:
            self.inoculation_DS_cost = round(dt.DS_cost * working_volume * 1e-3, 2)
        else:
            self.inoculation_DS_cost = 0

        self.ino_total_medium_cost = round(self.incubation_medium_cost + self.inoculation_medium_cost, 2)

        self.ino_total_non_medium_cost = round(self.incubation_Y_27632_cost + self.ino_DPBS_cost
                                           + self.ino_accutase_cost + self.bioreactors_cost
                                           + self.inoculation_Y_27632_cost + self.inoculation_DS_cost, 2)

        # Organization of Costs into Categories
        self.ino_consumables_cost = self.bioreactors_cost

        self.ino_reagent_cost = round(self.ino_total_medium_cost + self.incubation_Y_27632_cost
                                       + self.ino_DPBS_cost + self.ino_accutase_cost
                                       + self.inoculation_Y_27632_cost + self.inoculation_DS_cost, 2)

        self.ino_medium_cost = self.ino_total_medium_cost

        self.ino_facility_cost = 0 # included in the expansion phase

        self.ino_labor_cost = 0 # included in the expansion phase

        self.ino_total_cost = round(self.ino_consumables_cost + self.ino_reagent_cost + self.ino_facility_cost
                                    + self.ino_labor_cost, 2)


class SecondaryInoculation():

    def __init__(self, dt, ft, cc, working_volume, vwb_types):

        # Worktime [minutes]
        self.ino_inoculation_worktime = 60 + 10 * vwb_types['Amount'].sum()

        self.ino_total_worktime = self.ino_inoculation_worktime

        # convert the worktime to hours and round it to half hours
        self.ino_total_worktime = (math.ceil(self.ino_total_worktime / 30) * 30) / 60

        # Calculation of Costs
        self.bioreactors_cost = vwb_types['Associated Costs (€)'].sum()

        self.inoculation_medium_cost = round(self.culture_medium_cost * working_volume * 1e-3, 2)

        self.inoculation_Y_27632_cost = round(dt.Y_27632_cost * working_volume * 1e-3, 2)

        if cc.DS_supplementation:
            self.inoculation_DS_cost = round(dt.DS_cost * working_volume * 1e-3, 2)
        else:
            self.inoculation_DS_cost = 0

        self.ino_total_medium_cost = self.inoculation_medium_cost

        self.ino_total_non_medium_cost = round(self.bioreactors_cost + self.inoculation_Y_27632_cost
                                           + self.inoculation_DS_cost, 2)

        # Organization of Costs into Categories
        self.ino_consumables_cost = self.bioreactors_cost

        self.ino_reagent_cost = round(self.ino_total_medium_cost + self.inoculation_Y_27632_cost
                                       + self.inoculation_DS_cost, 2)

        self.ino_medium_cost = self.ino_total_medium_cost

        self.ino_facility_cost = 0 # included in the expansion phase

        self.ino_labor_cost = 0 # included in the expansion phase

        self.ino_total_cost = round(self.ino_consumables_cost + self.ino_reagent_cost + self.ino_facility_cost
                                    + self.ino_labor_cost, 2)


class Expansion():

    def __init__(self, dt, ft, cc, working_volume, vwb_types, min_fold_expansion):

        # Expansion of hiPSCs
        if cc.simulation_source == 'Nogueira':
            if cc.culture_medium == 'B8':
                self.exp_total_medium_volume = cc.Nogueira_B8_feeding_strategy(working_volume)
            else:
                self.exp_total_medium_volume = cc.Nogueira_feeding_strategy(working_volume)
        elif cc.simulation_source == 'Borys':
            self.exp_total_medium_volume = cc.Borys_feeding_strategy(working_volume)

        # determine minimum cell amount at the end of the current expansion cycle
        self.exp_cycle_results = {'Working Volume (mL)': [], 'Final Cell Number': []}

        vwb_index = []
        for vwb_type in vwb_types.index:
            for vwb in range(vwb_types.at[vwb_type, 'Amount']):
                vwb_index.append(vwb_type)
                self.exp_cycle_results['Working Volume (mL)'].append(vwb_types.at[vwb_type, 'Working Volume (mL)'])
                min_cell_number = (cc.seeding_density * vwb_types.at[vwb_type, 'Working Volume (mL)']
                               * min_fold_expansion)
                self.exp_cycle_results['Final Cell Number'].append(math.ceil(min_cell_number))

        self.exp_cycle_results = pd.DataFrame(self.exp_cycle_results, index=vwb_index)

        self.min_cell_number = self.exp_cycle_results['Final Cell Number'].sum()

        # Duration [days]
        self.exp_duration = cc.expansion_duration

        # Worktime [minutes]
        if cc.simulation_source == 'Nogueira':
            self.exp_medium_exchange_worktime = (5 + 10 * vwb_types['Amount'].sum()) * (cc.expansion_duration - 2)
        elif cc.simulation_source == 'Borys':
            self.exp_medium_exchange_worktime = 5 + 10 * vwb_types['Amount'].sum()

        self.exp_total_worktime = self.exp_medium_exchange_worktime

        # convert the worktime to hours and round it to half hours
        self.exp_total_worktime = (math.ceil(self.exp_medium_exchange_worktime / 30) * 30) / 60

        # Calculation of Costs
        self.exp_total_medium_cost = round(self.culture_medium_cost * self.exp_total_medium_volume * 1e-3, 2)

        self.exp_total_non_medium_cost = 0

        # Organization of Costs into Categories
        self.exp_consumables_cost = 0

        self.exp_reagent_cost = self.exp_total_medium_cost

        self.exp_medium_cost = self.exp_total_medium_cost

        self.exp_facility_cost = round(ft.overall_daily_facility_costs * self.exp_duration
                                       / ft.parallel_processes, 2)

        self.exp_labor_cost = round(self.exp_duration * ft.daily_employee_costs / ft.parallel_processes, 2)

        self.exp_total_cost = round(self.exp_consumables_cost + self.exp_reagent_cost + self.exp_facility_cost
                                    + self.exp_labor_cost, 2)


class Harvesting():

    def __init__(self, dt, ft, cc, vwb_types, min_multiplied_distributions):

        # Harvesting of hiPSCs
        self.har_accutase_volume = (vwb_types.at['PBS 0.1MAG', 'Amount']
                * dt.vwbs_parameters.at['PBS 0.1MAG', 'Max Vol (mL)'] + vwb_types.at['PBS 0.5MAG', 'Amount']
                * dt.vwbs_parameters.at['PBS 0.5MAG', 'Max Vol (mL)'] + vwb_types.at['PBS 3MAG', 'Amount']
                * dt.vwbs_parameters.at['PBS 3MAG', 'Max Vol (mL)']) * 0.2
            # the volume of accutase for harvesting is 20% of vwb maximum volume

        # determine minimum cell amount after harvesting
        self.exp_cycle_results['Final Cell Number'] = (cc.seeding_density
               * self.exp_cycle_results['Working Volume (mL)'] * min_multiplied_distributions)

        self.min_cell_number = self.exp_cycle_results['Final Cell Number'].sum()

        # Worktime [minutes]
        self.har_harvesting_worktime = 25 + 15 * vwb_types['Amount'].sum()

        self.har_total_worktime = self.har_harvesting_worktime

        # convert the worktime to hours and round it to half hours
        self.har_total_worktime = (math.ceil(self.har_total_worktime / 30) * 30) / 60

        # Calculation of Costs
        self.har_accutase_cost = round(dt.accutase_cost * self.har_accutase_volume * 1e-3, 2)

        self.har_Y_27632_cost = round(dt.Y_27632_cost * self.har_accutase_volume * 1e-3, 2)

        self.har_total_medium_cost = 0

        self.har_total_non_medium_cost = round(self.har_accutase_cost + self.har_Y_27632_cost, 2)

        # Organization of Costs into Categories
        self.har_consumables_cost = 0

        self.har_reagent_cost = self.har_total_non_medium_cost

        self.har_medium_cost = 0

        self.har_facility_cost = round(self.har_total_worktime * ft.overall_daily_facility_costs / 24, 2)

        self.har_labor_cost = round(self.har_total_worktime * dt.salary_per_hour, 2)

        self.har_total_cost = round(self.har_consumables_cost + self.har_reagent_cost + self.har_facility_cost
                                    + self.har_labor_cost, 2)


class IniQualControl():

    def __init__(self, dt, ft, cc):

        # Flow Cytometry
        # intracellular staining
        self.analyzed_volume = 2 * 0.3
            # 300 uL to test for OCT4 and another to use as a control

        self.iniqual_PBS_volume = 17 * self.analyzed_volume

        self.iniqual_PFA_volume = self.analyzed_volume

        self.iniqual_NGS1_volume = 9 * self.analyzed_volume

        self.iniqual_NGS3_volume = 3 * self.analyzed_volume

        self.iniqual_saponin_volume = self.analyzed_volume

        self.iniqual_OCT4_antibody_volume = self.analyzed_volume / 2

        self.iniqual_secondary_antibody_volume = self.analyzed_volume

        # extracelular staining
        self.analyzed_volume = 3 * 0.1
            # 100 uL to test each of two surface markers and another to use as a control

        self.iniqual_DPBS_volume = 7 * self.analyzed_volume

        self.iniqual_FBS_volume = self.analyzed_volume

        self.iniqual_TRA_1_60_PE_antibody_volume = self.analyzed_volume / 3

        self.iniqual_SSEA_1_FITC_antibody_volume = self.analyzed_volume / 3

        # costs
        self.iniqual_PBS_cost = round(dt.DPBS_cost * self.iniqual_DPBS_volume * 1e-3, 2)

        self.iniqual_PFA_cost = round(dt.fc_PFA_cost * self.iniqual_PFA_volume * 1e-3, 2)

        self.iniqual_NGS1_cost = round(dt.NGS1_cost * self.iniqual_NGS1_volume * 1e-3, 2)

        self.iniqual_NGS3_cost = round(dt.NGS3_cost * self.iniqual_NGS3_volume * 1e-3, 2)

        self.iniqual_saponin_cost = round(dt.saponin_cost * self.iniqual_saponin_volume * 1e-3, 2)

        self.iniqual_OCT4_antibody_cost = round(dt.fc_OCT4_antibody_cost
                                               * self.iniqual_OCT4_antibody_volume * 1e-3, 2)

        self.iniqual_DPBS_cost = round(dt.DPBS_cost * self.iniqual_DPBS_volume * 1e-3, 2)

        self.iniqual_FBS_cost = round(dt.FBS_cost * self.iniqual_FBS_volume * 1e-3, 2)

        self.iniqual_TRA_1_60_PE_antibody_cost = round(dt.fc_TRA_1_60_PE_antibody_cost
                                                   * self.iniqual_TRA_1_60_PE_antibody_volume * 1e-3, 2)

        self.iniqual_SSEA_1_FITC_antibody_cost = round(dt.fc_SSEA_1_FITC_antibody_cost
                                                 * self.iniqual_SSEA_1_FITC_antibody_volume * 1e-3, 2)

        self.iniqual_secondary_antibody_cost = round(dt.fc_secondary_antibody_cost
                                                    * self.iniqual_secondary_antibody_volume * 1e-3, 2)

        # Differentiation
        self.iniqual_well_amount = 3

        self.differentiation_medium_volume_per_well = 1.5

        self.PFA_volume_per_well = 1

        self.methanol_volume_per_well = 1

        self.FCB1_volume_per_well = 6 + 0.6

        self.FCB2_volume_per_well = 0.1 + 8 + 12

        self.cTnT_antibody_volume_per_well = 0.1

        self.secondary_antibody_volume_per_well = 0.2

        self.expansion_medium_volume = (dt.peps_parameters.at['12-well', 'Medium Vol (mL)']
                                        * self.iniqual_well_amount * dt.well_culture_duration)

        self.differentiation_medium_volume = self.differentiation_medium_volume_per_well * self.iniqual_well_amount

        self.PFA_volume = self.PFA_volume_per_well * self.iniqual_well_amount

        self.methanol_volume = self.methanol_volume_per_well * self.iniqual_well_amount

        self.FCB1_volume = self.FCB1_volume_per_well * self.iniqual_well_amount

        self.FCB2_volume = self.FCB2_volume_per_well * self.iniqual_well_amount

        self.cTnT_antibody_volume = self.cTnT_antibody_volume_per_well * self.iniqual_well_amount

        self.secondary_antibody_volume = self.secondary_antibody_volume_per_well * self.iniqual_well_amount

        # costs
        self.expansion_medium_cost =  round(self.culture_medium_cost * self.expansion_medium_volume * 1e-3, 2)

        self.day0_cost = round((dt.RPMI_cost + dt.B27_no_insulin_cost + dt.CHIR99021_cost)
                               * self.differentiation_medium_volume * 1e-3, 2)
            # RPMI/B27-insulin with CHIR99021

        self.day1_cost = round((dt.RPMI_cost + dt.B27_no_insulin_cost)
                               * self.differentiation_medium_volume * 1e-3, 2)
            # RPMI/B27-insulin

        self.day3_cost = round(((dt.RPMI_cost + dt.B27_no_insulin_cost) / 2 + dt.IWP4_cost)
                               * self.differentiation_medium_volume * 1e-3, 2)
            # RPMI/B27-insulin with IWP4

        self.day5_cost = round((dt.RPMI_cost + dt.B27_no_insulin_cost)
                               * self.differentiation_medium_volume * 1e-3, 2)
            # RPMI/B27-insulin

        self.day7_cost = round((dt.RPMI_cost + dt.B27_cost) * self.differentiation_medium_volume * 1e-3, 2)
            # RPMI/B27

        self.day10_cost = self.day7_cost
            # RPMI/B27

        self.day13_cost = self.day7_cost
            # RPMI/B27

        self.PSA_cost = round((dt.PBS_cost + dt.fc_PFA_cost) * self.PFA_volume * 1e-3, 2)

        self.methanol_cost = round(dt.methanol_cost * self.methanol_volume * 1e-3, 2)

        self.FCB1_cost = round((dt.PBS_cost + dt.BSA_05_cost) * self.FCB1_volume * 1e-3, 2)

        self.FCB2_cost = round((dt.PBS_cost + dt.BSA_05_cost + dt.triton_cost) * self.FCB2_volume * 1e-3, 2)

        self.cTnT_antibody_cost = round((dt.PBS_cost + dt.BSA_05_cost + dt.triton_cost + dt.cTnT_antibody_cost)
                                        * self.cTnT_antibody_volume * 1e-3, 2)

        self.secondary_antibody_cost = round((dt.PBS_cost + dt.BSA_05_cost + dt.triton_cost
                                         + dt.secondary_antibody_cost) * self.secondary_antibody_volume * 1e-3, 2)

        self.iniqual_differentiation_cost = round(self.expansion_medium_cost + self.day0_cost + self.day1_cost
               + self.day3_cost + self.day5_cost + self.day7_cost + self.day10_cost + self.day13_cost
               + self.PSA_cost + self.methanol_cost + self.FCB1_cost + self.FCB2_cost + self.cTnT_antibody_cost
               + self.secondary_antibody_cost)

        # Worktime [minutes]
        self.iniqual_flow_cytometry_worktime = 60 + 30 + 2 * 30

        self.iniqual_differentiation_worktime = (5 + self.iniqual_well_amount * 1 * dt.well_culture_duration
               + 10 + math.ceil(self.iniqual_well_amount / 12) * 5 + 60 + 30 + 2 * 30)

        self.iniqual_total_worktime = self.iniqual_flow_cytometry_worktime + self.iniqual_differentiation_worktime

        # convert the worktime to hours and round it to half hours
        self.iniqual_total_worktime = (math.ceil(self.iniqual_total_worktime / 30) * 30) / 60

        # Calculation of Costs
        self.iniqual_total_medium_cost = self.expansion_medium_cost

        self.iniqual_total_non_medium_cost = round(self.iniqual_PBS_cost + self.iniqual_PFA_cost
             + self.iniqual_NGS1_cost + self.iniqual_NGS3_cost + self.iniqual_saponin_cost
             + self.iniqual_OCT4_antibody_cost + self.iniqual_DPBS_cost + self.iniqual_FBS_cost
             + self.iniqual_TRA_1_60_PE_antibody_cost + self.iniqual_SSEA_1_FITC_antibody_cost
             + self.iniqual_secondary_antibody_cost + self.iniqual_differentiation_cost
             - self.expansion_medium_cost, 2)

        # Organization of Costs into Categories
        self.iniqual_consumables_cost = 0

        self.iniqual_reagent_cost = self.iniqual_total_medium_cost + self.iniqual_total_non_medium_cost

        self.iniqual_medium_cost = self.iniqual_total_medium_cost

        self.iniqual_facility_cost = round(self.iniqual_total_worktime * ft.overall_daily_facility_costs / 24, 2)

        self.iniqual_labor_cost = round(self.iniqual_total_worktime * dt.salary_per_hour, 2)

        self.iniqual_total_cost = round(self.iniqual_consumables_cost + self.iniqual_reagent_cost
                                       + self.iniqual_facility_cost + self.iniqual_labor_cost, 2)


class InterQualControl():

    def __init__(self, dt, ft, cc, vwb_types):

        # Flow Cytometry
        # intracellular staining
        self.analyzed_volume = 2 * 0.3 * vwb_types['Amount'].sum()
            # 300 uL to test for OCT4 and another to use as a control

        self.inqual_PBS_volume = 17 * self.analyzed_volume

        self.inqual_PFA_volume = self.analyzed_volume

        self.inqual_NGS1_volume = 9 * self.analyzed_volume

        self.inqual_NGS3_volume = 3 * self.analyzed_volume

        self.inqual_saponin_volume = self.analyzed_volume

        self.inqual_OCT4_antibody_volume = self.analyzed_volume / 2

        self.inqual_secondary_antibody_volume = self.analyzed_volume

        # extracelular staining
        self.analyzed_volume = 3 * 0.1 * vwb_types['Amount'].sum()
            # 100 uL to test each of two surface markers and another to use as a control

        self.inqual_DPBS_volume = 7 * self.analyzed_volume

        self.inqual_FBS_volume = self.analyzed_volume

        self.inqual_TRA_1_60_PE_antibody_volume = self.analyzed_volume / 3

        self.inqual_SSEA_1_FITC_antibody_volume = self.analyzed_volume / 3

        # costs
        self.inqual_PBS_cost = round(dt.DPBS_cost * self.inqual_DPBS_volume * 1e-3, 2)

        self.inqual_PFA_cost = round(dt.fc_PFA_cost * self.inqual_PFA_volume * 1e-3, 2)

        self.inqual_NGS1_cost = round(dt.NGS1_cost * self.inqual_NGS1_volume * 1e-3, 2)

        self.inqual_NGS3_cost = round(dt.NGS3_cost * self.inqual_NGS3_volume * 1e-3, 2)

        self.inqual_saponin_cost = round(dt.saponin_cost * self.inqual_saponin_volume * 1e-3, 2)

        self.inqual_OCT4_antibody_cost = round(dt.fc_OCT4_antibody_cost
                                               * self.inqual_OCT4_antibody_volume * 1e-3, 2)

        self.inqual_DPBS_cost = round(dt.DPBS_cost * self.inqual_DPBS_volume * 1e-3, 2)

        self.inqual_FBS_cost = round(dt.FBS_cost * self.inqual_FBS_volume * 1e-3, 2)

        self.inqual_TRA_1_60_PE_antibody_cost = round(dt.fc_TRA_1_60_PE_antibody_cost
                                                   * self.inqual_TRA_1_60_PE_antibody_volume * 1e-3, 2)

        self.inqual_SSEA_1_FITC_antibody_cost = round(dt.fc_SSEA_1_FITC_antibody_cost
                                                 * self.inqual_SSEA_1_FITC_antibody_volume * 1e-3, 2)

        self.inqual_secondary_antibody_cost = round(dt.fc_secondary_antibody_cost
                                                    * self.inqual_secondary_antibody_volume * 1e-3, 2)

        # Worktime [minutes]
        self.inqual_flow_cytometry_worktime = 60 + 30 + 2 * 30 * vwb_types['Amount'].sum()

        self.inqual_total_worktime = self.inqual_flow_cytometry_worktime

        # convert the worktime to hours and round it to half hours
        self.inqual_total_worktime = (math.ceil(self.inqual_total_worktime / 30) * 30) / 60

        # Calculation of Costs
        self.inqual_total_medium_cost = 0

        self.inqual_total_non_medium_cost = round(self.inqual_PBS_cost + self.inqual_PFA_cost
             + self.inqual_NGS1_cost + self.inqual_NGS3_cost + self.inqual_saponin_cost
             + self.inqual_OCT4_antibody_cost + self.inqual_DPBS_cost + self.inqual_FBS_cost
             + self.inqual_TRA_1_60_PE_antibody_cost + self.inqual_SSEA_1_FITC_antibody_cost
             + self.inqual_secondary_antibody_cost, 2)

        # Organization of Costs into Categories
        self.inqual_consumables_cost = 0

        self.inqual_reagent_cost = self.inqual_total_non_medium_cost

        self.inqual_medium_cost = 0

        self.inqual_facility_cost = round(self.inqual_total_worktime * ft.overall_daily_facility_costs / 24, 2)

        self.inqual_labor_cost = round(self.inqual_total_worktime * dt.salary_per_hour, 2)

        self.inqual_total_cost = round(self.inqual_consumables_cost + self.inqual_reagent_cost
                                       + self.inqual_facility_cost + self.inqual_labor_cost, 2)


class FinalQualControl():

    def __init__(self, dt, ft, cc, vwb_types):

        # qRT-PCR
        self.finqual_qRT_PCR_kits_cost = (dt.RNA_purification_kit_cost + dt.cDNA_transcriptase_kit_cost
                                  + dt.qPCR_master_mix_kit_cost * 3)

        self.finqual_primers_cost = (dt.OCT4_primer_cost + dt.NANOG_primer_cost + dt.SOX1_primer_cost
                                     + dt.T_primer_cost + dt.SOX17_primer_cost) * 3

        self.finqual_qRT_PCR_cost = self.finqual_qRT_PCR_kits_cost + self.finqual_primers_cost

        # Immunocytochemistry
        # intracellular staining
        self.analyzed_volume = 3 * 0.2 * vwb_types['Amount'].sum()
            # 200 uL to test each of two internal markers and another to use as a control

        self.finqual_PBS_volume = 18 * self.analyzed_volume

        self.finqual_PFA_volume = self.analyzed_volume

        self.finqual_triton_volume = 3 * self.analyzed_volume

        self.finqual_NGS10_volume = self.analyzed_volume

        self.finqual_NGS5_volume = 2 * self.analyzed_volume

        self.finqual_OCT4_antibody_volume = self.analyzed_volume / 3

        self.finqual_SOX2_antibody_volume = self.analyzed_volume / 3

        self.finqual_secondary_antibody_volume = self.analyzed_volume

        self.finqual_DAPI_volume = self.analyzed_volume

        # extracelular staining
        self.analyzed_volume = 3 * 0.2 * vwb_types['Amount'].sum()
            # 200 uL to test each of two surface markers and another to use as a control

        self.finqual_DPBS_volume = 11 * self.analyzed_volume

        self.finqual_BSA_volume = 2 * self.analyzed_volume

        self.finqual_TRA_1_60_antibody_volume = self.analyzed_volume / 3

        self.finqual_SSEA_4_antibody_volume = self.analyzed_volume / 3

        self.finqual_secondary_antibody_volume += self.analyzed_volume

        # costs
        self.finqual_PBS_cost = round(dt.PBS_cost * self.finqual_PBS_volume * 1e-3, 2)

        self.finqual_PFA_cost = round(dt.ic_PFA_cost * self.finqual_PFA_volume * 1e-3, 2)

        self.finqual_triton_cost = round(dt.triton_cost * self.finqual_triton_volume * 1e-3, 2)

        self.finqual_NGS10_cost = round(dt.NGS10_cost * self.finqual_NGS10_volume * 1e-3, 2)

        self.finqual_NGS5_cost = round(dt.NGS5_cost * self.finqual_NGS5_volume * 1e-3, 2)

        self.finqual_OCT4_antibody_cost = round(dt.ic_OCT4_antibody_cost
                                                * self.finqual_OCT4_antibody_volume * 1e-3, 2)

        self.finqual_SOX2_antibody_cost = round(dt.ic_SOX2_antibody_cost
                                                * self.finqual_SOX2_antibody_volume * 1e-3, 2)

        self.finqual_DPBS_cost = round(dt.DPBS_cost * self.finqual_DPBS_volume * 1e-3, 2)

        self.finqual_BSA_cost = round(dt.BSA_cost * self.finqual_BSA_volume * 1e-3, 2)

        self.finqual_TRA_1_60_antibody_cost = round(dt.ic_TRA_1_60_antibody_cost
                                                   * self.finqual_TRA_1_60_antibody_volume * 1e-3, 2)

        self.finqual_SSEA_4_antibody_cost = round(dt.ic_SSEA_4_antibody_cost
                                                 * self.finqual_SSEA_4_antibody_volume * 1e-3, 2)

        self.finqual_secondary_antibody_cost = round(dt.ic_secondary_antibody_cost
                                                    * self.finqual_secondary_antibody_volume * 1e-3, 2)

        self.finqual_immunocytochemistry_cost = (self.finqual_PBS_cost + self.finqual_PFA_cost
             + self.finqual_triton_cost + self.finqual_NGS10_cost + self.finqual_NGS5_cost
             + self.finqual_OCT4_antibody_cost + self.finqual_SOX2_antibody_cost + self.finqual_DPBS_cost
             + self.finqual_BSA_cost + self.finqual_TRA_1_60_antibody_cost + self.finqual_SSEA_4_antibody_cost
             + self.finqual_secondary_antibody_cost)

        # Genetic Analysis
        self.finqual_genetic_analysis_cost = dt.DNA_purification_kit_cost + dt.hPSC_genetic_kit_cost

        # G Banding
        self.analyzed_volume = 2 * vwb_types['Amount'].sum()
            # 2 mL to get enough cells for karyotyping

        self.finqual_colcemid_volume = self.analyzed_volume
            # 0.5 ug of colcemid per mL (see concentration of colcemid in database)

        self.finqual_KCl_volume = 10

        self.finqual_gurr_buffer_amount = 2 * vwb_types['Amount'].sum()
            # two tablets of gurr buffer are used

        self.finqual_trypsin_volume = 2.5 * vwb_types['Amount'].sum()

        self.finqual_giemsa_stain_volume = 3 * vwb_types['Amount'].sum()

        # costs
        self.finqual_colcemid_cost = round(dt.colcemid_cost * self.finqual_colcemid_volume * 1e-3, 2)

        self.finqual_KCl_cost = round(dt.KCl_cost * self.finqual_KCl_volume * 1e-3, 2)

        self.finqual_gurr_buffer_cost = round(dt.gurr_buffer_cost * self.finqual_gurr_buffer_amount, 2)

        self.finqual_trypsin_cost = round(dt.trypsin_cost * self.finqual_trypsin_volume * 1e-3, 2)

        self.finqual_giemsa_stain_cost = round(dt.giemsa_stain_cost * self.finqual_giemsa_stain_volume * 1e-3, 2)

        self.finqual_g_banding_cost = (self.finqual_colcemid_cost + self.finqual_KCl_cost
             + self.finqual_gurr_buffer_cost + self.finqual_trypsin_cost + self.finqual_giemsa_stain_cost)

        # Differentiation
        self.finqual_well_amount = 3 * vwb_types['Amount'].sum()

        self.differentiation_medium_volume_per_well = 1.5

        self.PFA_volume_per_well = 1

        self.methanol_volume_per_well = 1

        self.FCB1_volume_per_well = 6 + 0.6

        self.FCB2_volume_per_well = 0.1 + 8 + 12

        self.cTnT_antibody_volume_per_well = 0.1

        self.secondary_antibody_volume_per_well = 0.2

        self.expansion_medium_volume = (dt.peps_parameters.at['12-well', 'Medium Vol (mL)']
                                        * self.finqual_well_amount * dt.well_culture_duration)

        self.differentiation_medium_volume = self.differentiation_medium_volume_per_well * self.finqual_well_amount

        self.PFA_volume = self.PFA_volume_per_well * self.finqual_well_amount

        self.methanol_volume = self.methanol_volume_per_well * self.finqual_well_amount

        self.FCB1_volume = self.FCB1_volume_per_well * self.finqual_well_amount

        self.FCB2_volume = self.FCB2_volume_per_well * self.finqual_well_amount

        self.cTnT_antibody_volume = self.cTnT_antibody_volume_per_well * self.finqual_well_amount

        self.secondary_antibody_volume = self.secondary_antibody_volume_per_well * self.finqual_well_amount

        # costs
        self.expansion_medium_cost =  round(self.culture_medium_cost * self.expansion_medium_volume * 1e-3, 2)

        self.day0_cost = round((dt.RPMI_cost + dt.B27_no_insulin_cost + dt.CHIR99021_cost)
                               * self.differentiation_medium_volume * 1e-3, 2)
            # RPMI/B27-insulin with CHIR99021

        self.day1_cost = round((dt.RPMI_cost + dt.B27_no_insulin_cost)
                               * self.differentiation_medium_volume * 1e-3, 2)
            # RPMI/B27-insulin

        self.day3_cost = round(((dt.RPMI_cost + dt.B27_no_insulin_cost) / 2 + dt.IWP4_cost)
                               * self.differentiation_medium_volume * 1e-3, 2)
            # RPMI/B27-insulin with IWP4

        self.day5_cost = round((dt.RPMI_cost + dt.B27_no_insulin_cost)
                               * self.differentiation_medium_volume * 1e-3, 2)
            # RPMI/B27-insulin

        self.day7_cost = round((dt.RPMI_cost + dt.B27_cost) * self.differentiation_medium_volume * 1e-3, 2)
            # RPMI/B27

        self.day10_cost = self.day7_cost
            # RPMI/B27

        self.day13_cost = self.day7_cost
            # RPMI/B27

        self.PSA_cost = round((dt.PBS_cost + dt.fc_PFA_cost) * self.PFA_volume * 1e-3, 2)

        self.methanol_cost = round(dt.methanol_cost * self.methanol_volume * 1e-3, 2)

        self.FCB1_cost = round((dt.PBS_cost + dt.BSA_05_cost) * self.FCB1_volume * 1e-3, 2)

        self.FCB2_cost = round((dt.PBS_cost + dt.BSA_05_cost + dt.triton_cost) * self.FCB2_volume * 1e-3, 2)

        self.cTnT_antibody_cost = round((dt.PBS_cost + dt.BSA_05_cost + dt.triton_cost + dt.cTnT_antibody_cost)
                                        * self.cTnT_antibody_volume * 1e-3, 2)

        self.secondary_antibody_cost = round((dt.PBS_cost + dt.BSA_05_cost + dt.triton_cost
                                         + dt.secondary_antibody_cost) * self.secondary_antibody_volume * 1e-3, 2)

        self.finqual_differentiation_cost = round(self.expansion_medium_cost + self.day0_cost + self.day1_cost
               + self.day3_cost + self.day5_cost + self.day7_cost + self.day10_cost + self.day13_cost
               + self.PSA_cost + self.methanol_cost + self.FCB1_cost + self.FCB2_cost + self.cTnT_antibody_cost
               + self.secondary_antibody_cost)

        # Duration [days]
        self.finqual_duration = cc.final_quality_control_duration

        # Worktime [minutes]
        self.finqual_qRT_PCR_worktime = 90 * vwb_types['Amount'].sum()

        self.finqual_immunocytochemistry_worktime = 60 + 30 + 2 * 30 * vwb_types['Amount'].sum()

        self.finqual_genetic_analysis_worktime = 120 + 30 * vwb_types['Amount'].sum()

        self.finqual_g_banding_worktime = 110 + 30 + 60 * vwb_types['Amount'].sum()

        self.pre_medium_exchange_worktime = ((5 * (self.passage_number + 1) + self.finqual_well_amount * 1) *
                                          dt.well_culture_duration)

        self.finqual_differentiation_worktime = (5 + self.finqual_well_amount * 1 * dt.well_culture_duration
               + 10 + math.ceil(self.finqual_well_amount / 12) * 5 + 60 + 30 + 2 * 30 * vwb_types['Amount'].sum())

        self.finqual_total_worktime = (self.finqual_qRT_PCR_worktime + self.finqual_immunocytochemistry_worktime
                                       + self.finqual_genetic_analysis_worktime + self.finqual_g_banding_worktime
                                       + self.finqual_differentiation_worktime)

        # convert the worktime to hours and round it to half hours
        self.finqual_total_worktime = (math.ceil(self.finqual_total_worktime / 30) * 30) / 60

        # Calculation of Costs
        self.finqual_total_medium_cost = self.expansion_medium_cost

        self.finqual_total_non_medium_cost = round(self.finqual_qRT_PCR_cost
               + self.finqual_immunocytochemistry_cost + self.finqual_genetic_analysis_cost
               + self.finqual_g_banding_cost + self.finqual_differentiation_cost - self.expansion_medium_cost, 2)

        # Organization of Costs into Categories
        self.finqual_consumables_cost = 0

        self.finqual_reagent_cost = self.finqual_total_medium_cost + self.finqual_total_non_medium_cost

        self.finqual_medium_cost = self.finqual_total_medium_cost

        self.finqual_facility_cost = round(ft.overall_daily_facility_costs * self.finqual_duration
                                       / ft.parallel_processes, 2)

        self.finqual_labor_cost = round(self.finqual_duration * ft.daily_employee_costs / ft.parallel_processes, 2)

        self.finqual_total_cost = round(self.finqual_consumables_cost + self.finqual_reagent_cost
                                       + self.finqual_facility_cost + self.finqual_labor_cost, 2)


class Bioprocess():

    def __init__(self, dt, ft, cc, fold_expansion_distribution, recovery_efficiency_distribution, init_sim_num):

        # define arrays for accumulating values
        cycles = ['2D Exp']

        cycle_durations = []
        cycle_worktimes = []
        cycle_final_cell_number = []
        cycle_medium_costs = []
        cycle_process_costs = []

        cycle_consumables_costs = []
        cycle_reagents_costs = []
        cycle_facility_costs = []
        cycle_labor_costs = []

        cycle_pre_inoculation_costs = []
        cycle_inoculation_costs = []
        cycle_expansion_costs = []
        cycle_qual_control_costs = []
        cycle_harvesting_costs = []

        # carry out optimal workflow calculations based on the culture conditions
        (self.necessary_well_amount, self.average_wells_cell_number, self.optimal_working_volumes,
         self.success_probability, self.average_cycle_cell_numbers,
         self.final_cycle_cell_number_distribution) = workflow_calculations(dt, cc)
        self.average_cycle_cell_numbers.insert(0, self.average_wells_cell_number)

        # organize the bioreactor workflow according to working volume milestones
        (self.expansion_cycles, self.bioreactor_workflow) = create_3D_workflow(dt, self.optimal_working_volumes)

        # run first expansion cycle and calculate associated costs
        PreInoculation.__init__(self, dt, ft, cc, self.necessary_well_amount)

        IniQualControl.__init__(self, dt, ft, cc)

        # detract direct labor and facility cost of the initial quality control from the pre inoculation phase
        self.pre_total_cost -= self.iniqual_facility_cost + self.iniqual_labor_cost
        self.pre_facility_cost -= self.iniqual_facility_cost
        self.pre_labor_cost -= self.iniqual_labor_cost

        # save values for planar expansion phase
        cycle_durations.append(self.pre_duration)
        cycle_worktimes.append(self.pre_total_worktime + self.iniqual_total_worktime)
        cycle_final_cell_number.append(self.average_wells_cell_number)
        cycle_medium_costs.append(self.pre_medium_cost + self.iniqual_medium_cost)
        cycle_process_costs.append(self.pre_total_cost + self.iniqual_total_cost)

        cycle_consumables_costs.append(self.pre_consumables_cost + self.iniqual_consumables_cost)
        cycle_reagents_costs.append(self.pre_reagent_cost + self.iniqual_reagent_cost)
        cycle_facility_costs.append(self.pre_facility_cost + self.iniqual_facility_cost)
        cycle_labor_costs.append(self.pre_labor_cost + self.iniqual_labor_cost)

        cycle_pre_inoculation_costs.append(self.pre_total_cost)
        cycle_inoculation_costs.append(0)
        cycle_expansion_costs.append(0)
        cycle_qual_control_costs.append(self.iniqual_total_cost)
        cycle_harvesting_costs.append(0)

        # determine minimum fold expansion and recovery efficiency according to threshold
        fold_expansion_distribution.sort()
        recovery_efficiency_distribution.sort()

        min_fold_expansion = fold_expansion_distribution[math.ceil(init_sim_num * (1 - cc.threshold))]
        min_recovery_efficiency = recovery_efficiency_distribution[math.ceil(init_sim_num * (1 - cc.threshold))]

        # determine product of the fold expansion and recovery efficiency distributions
        shuffle(fold_expansion_distribution)
        shuffle(recovery_efficiency_distribution)
        fold_expansion_distribution = np.array(fold_expansion_distribution)
        recovery_efficiency_distribution = np.array(recovery_efficiency_distribution)

        multiplied_distributions = fold_expansion_distribution * recovery_efficiency_distribution
        multiplied_distributions.sort()
        min_multiplied_distributions = multiplied_distributions[math.ceil(init_sim_num * (1 - cc.threshold))]

        PrimaryInoculation.__init__(self, dt, ft, cc, self.optimal_working_volumes[0], self.bioreactor_workflow[0])

        # run the following expansion cycles, according to the predetermined optimal workflow
        for expansion_cycle, current_workflow in enumerate(self.bioreactor_workflow):
            # carry out bioreactor-bioreactor inoculation if currently not the first expansion phase
            if expansion_cycle != 0:
                SecondaryInoculation.__init__(self, dt, ft, cc, self.optimal_working_volumes[expansion_cycle],
                                              current_workflow)

            Expansion.__init__(self, dt, ft, cc, self.optimal_working_volumes[expansion_cycle], current_workflow,
                               min_fold_expansion)

            InterQualControl.__init__(self, dt, ft, cc, current_workflow)

            # detract direct labor and facility cost of the intermediate quality control from the expansion phase
            self.exp_total_cost -= self.inqual_facility_cost + self.inqual_labor_cost
            self.exp_facility_cost -= self.inqual_facility_cost
            self.exp_labor_cost -= self.inqual_labor_cost

            # carry out the harvesting phase if currently not in the last expansion phase
            if expansion_cycle < len(self.bioreactor_workflow) - 1:
                Harvesting.__init__(self, dt, ft, cc, current_workflow, min_multiplied_distributions)

                # detract direct labor and facility cost of the harvesting phase from the expansion phase
                self.exp_total_cost -= self.har_facility_cost + self.har_labor_cost
                self.exp_facility_cost -= self.har_facility_cost
                self.exp_labor_cost -= self.har_labor_cost

            cycles.append(f'C{expansion_cycle}')

            # consider cost of the harvesting phase if currently not in the last expansion phase
            if expansion_cycle < len(self.bioreactor_workflow) - 1:
                # save values for current expansion cycle (Ci)
                cycle_durations.append(self.exp_duration)
                cycle_worktimes.append(self.har_total_worktime + self.ino_total_worktime + self.exp_total_worktime
                                       + self.inqual_total_worktime)
                cycle_final_cell_number.append(self.min_cell_number)
                cycle_medium_costs.append(self.har_medium_cost + self.ino_medium_cost + self.exp_medium_cost
                                          + self.inqual_medium_cost)
                cycle_process_costs.append(self.har_total_cost + self.ino_total_cost + self.exp_total_cost
                                           + self.inqual_total_cost)

                cycle_consumables_costs.append(self.har_consumables_cost + self.ino_consumables_cost
                                               + self.exp_consumables_cost + self.inqual_consumables_cost)
                cycle_reagents_costs.append(self.har_reagent_cost + self.ino_reagent_cost + self.exp_reagent_cost
                                            + self.inqual_reagent_cost)
                cycle_facility_costs.append(self.har_facility_cost + self.ino_facility_cost + self.exp_facility_cost
                                            + self.inqual_facility_cost)
                cycle_labor_costs.append(self.har_labor_cost + self.ino_labor_cost + self.exp_labor_cost
                                         + self.inqual_labor_cost)

                cycle_pre_inoculation_costs.append(0)
                cycle_inoculation_costs.append(self.ino_total_cost)
                cycle_expansion_costs.append(self.exp_total_cost)
                cycle_qual_control_costs.append(self.inqual_total_cost)
                cycle_harvesting_costs.append(self.har_total_cost)

            else:
                # save values for current expansion cycle (Ci)
                cycle_durations.append(self.exp_duration)
                cycle_worktimes.append(self.ino_total_worktime + self.exp_total_worktime + self.inqual_total_worktime)
                cycle_final_cell_number.append(self.min_cell_number)
                cycle_medium_costs.append(self.ino_medium_cost + self.exp_medium_cost + self.inqual_medium_cost)
                cycle_process_costs.append(self.ino_total_cost + self.exp_total_cost + self.inqual_total_cost)

                cycle_consumables_costs.append(self.ino_consumables_cost + self.exp_consumables_cost
                                               + self.inqual_consumables_cost)
                cycle_reagents_costs.append(self.ino_reagent_cost + self.exp_reagent_cost + self.inqual_reagent_cost)
                cycle_facility_costs.append(self.ino_facility_cost + self.exp_facility_cost
                                            + self.inqual_facility_cost)
                cycle_labor_costs.append(self.ino_labor_cost + self.exp_labor_cost + self.inqual_labor_cost)

                cycle_pre_inoculation_costs.append(0)
                cycle_inoculation_costs.append(self.ino_total_cost)
                cycle_expansion_costs.append(self.exp_total_cost)
                cycle_qual_control_costs.append(self.inqual_total_cost)

            # move on to the next expansion cycle
            expansion_cycle += 1

        # no harvesting after last expansion cycle
        cycle_harvesting_costs.append(0)

        # take into account the final quality control phase
        FinalQualControl.__init__(self, dt, ft, cc, self.bioreactor_workflow[-1])

        cycle_durations[-1] += self.finqual_duration
        cycle_worktimes[-1] += self.finqual_total_worktime
        cycle_medium_costs[-1] += self.finqual_medium_cost
        cycle_process_costs[-1] += self.finqual_total_cost

        cycle_consumables_costs[-1] += self.finqual_consumables_cost
        cycle_reagents_costs[-1] += self.finqual_reagent_cost
        cycle_facility_costs[-1] += self.finqual_facility_cost
        cycle_labor_costs[-1] += self.finqual_labor_cost

        cycle_qual_control_costs[-1] += self.finqual_total_cost

        # create dataframes for presentation of data to user
        self.cycle_info = pd.DataFrame({'Duration (d)': cycle_durations, 'Worktime (h)': cycle_worktimes,
                'Cell Number': self.average_cycle_cell_numbers, 'Medium Cost (€)': cycle_medium_costs,
                'Phase Cost (€)': cycle_process_costs}, index=cycles)

        self.cycle_cost_categories = pd.DataFrame({'Consumables': cycle_consumables_costs, 'Reagents':
                cycle_reagents_costs, 'Facility': cycle_facility_costs, 'Labor': cycle_labor_costs}, index=cycles)

        self.cycle_cost_stages = pd.DataFrame({'Pre-Inoculation': cycle_pre_inoculation_costs,
                'Expansion': np.add(cycle_inoculation_costs, cycle_expansion_costs), 'Quality Control':
                cycle_qual_control_costs, 'Harvesting': cycle_harvesting_costs}, index=cycles)

        # determine accumulated values over the entire bioprocess
        self.total_duration = self.cycle_info['Duration (d)'].sum()
        self.total_worktime = math.ceil(self.cycle_info['Worktime (h)'].sum())
        self.final_cell_number = self.average_cycle_cell_numbers[-1]
        self.total_medium_cost = self.cycle_info['Medium Cost (€)'].sum()

        self.total_consumables_cost = self.cycle_cost_categories['Consumables'].sum()
        self.total_reagents_cost = self.cycle_cost_categories['Reagents'].sum()
        self.total_facility_cost = self.cycle_cost_categories['Facility'].sum()
        self.total_labor_cost = self.cycle_cost_categories['Labor'].sum()

        self.total_pre_inoculation_cost = self.cycle_cost_stages['Pre-Inoculation'].sum()
        self.total_expansion_cost = self.cycle_cost_stages['Expansion'].sum()
        self.total_quality_control_cost = self.cycle_cost_stages['Quality Control'].sum()
        self.total_harvesting_cost = self.cycle_cost_stages['Harvesting'].sum()

        # calculate direct, indirect and total cost per bioprocess
        self.total_direct_cost = round(self.total_consumables_cost + self.total_reagents_cost, 2)
        self.total_indirect_cost = round(self.total_facility_cost + self.total_labor_cost, 2)

        self.total_cost = round(self.total_direct_cost + self.total_indirect_cost, 2)

        # take into account potential failures and distribute failed bioprocess costs among successful bioprocesses
        failure_detection_probs = []
        qual_ctrl_checkpoints = len(cycles)
        failure_occurence_prob = 1 / qual_ctrl_checkpoints

        for index in range(qual_ctrl_checkpoints - 1):
            i = 0
            failure_detection_prob = 0
            while i <= index:
                failure_detection_prob += ((cc.bioprocess_failure_prob * failure_occurence_prob)
                                           * (1 - cc.qual_ctrl_detection_prob)**i)
                i += 1
            failure_detection_prob *= cc.qual_ctrl_detection_prob
            failure_detection_probs.append(failure_detection_prob)

        failure_detection_probs.append(cc.bioprocess_failure_prob - sum(failure_detection_probs))
        #print(f"Failure: {failure_detection_probs}")

        self.cumulative_cycle_info = self.cycle_info.cumsum()
        self.cumulative_cycle_cost_categories = self.cycle_cost_categories.cumsum()
        #print(self.cumulative_cycle_cost_categories)
        self.cumulative_cycle_cost_stages = self.cycle_cost_stages.cumsum()
        for index, cycle in enumerate(cycles):
            #print(self.total_consumables_cost)
            self.total_medium_cost += (self.cumulative_cycle_info.at[cycle, 'Medium Cost (€)']
                                       * failure_detection_probs[index] / (1 - cc.bioprocess_failure_prob))

            self.total_consumables_cost += (self.cumulative_cycle_cost_categories.at[cycle, 'Consumables']
                                            * failure_detection_probs[index] / (1 - cc.bioprocess_failure_prob))
            self.total_reagents_cost += (self.cumulative_cycle_cost_categories.at[cycle, 'Reagents']
                                         * failure_detection_probs[index] / (1 - cc.bioprocess_failure_prob))
            self.total_facility_cost += (self.cumulative_cycle_cost_categories.at[cycle, 'Facility']
                                         * failure_detection_probs[index] / (1 - cc.bioprocess_failure_prob))
            self.total_labor_cost += (self.cumulative_cycle_cost_categories.at[cycle, 'Labor']
                                      * failure_detection_probs[index] / (1 - cc.bioprocess_failure_prob))

            self.total_pre_inoculation_cost += (self.cumulative_cycle_cost_stages.at[cycle, 'Pre-Inoculation']
                                        * failure_detection_probs[index] / (1 - cc.bioprocess_failure_prob))
            self.total_expansion_cost += (self.cumulative_cycle_cost_stages.at[cycle, 'Expansion']
                                  * failure_detection_probs[index] / (1 - cc.bioprocess_failure_prob))
            self.total_quality_control_cost += (self.cumulative_cycle_cost_stages.at[cycle, 'Quality Control']
                                        * failure_detection_probs[index] / (1 - cc.bioprocess_failure_prob))
            self.total_harvesting_cost += (self.cumulative_cycle_cost_stages.at[cycle, 'Harvesting']
                                   * failure_detection_probs[index] / (1 - cc.bioprocess_failure_prob))

        #print(self.total_consumables_cost)
        # consider cases where the target cell number was not reached as failures
        self.total_medium_cost /= self.success_probability

        self.total_consumables_cost /= self.success_probability
        self.total_reagents_cost /= self.success_probability
        self.total_facility_cost /= self.success_probability
        self.total_labor_cost /= self.success_probability
        #print(self.total_consumables_cost, self.total_reagents_cost, self.total_facility_cost, self.total_labor_cost)

        self.total_pre_inoculation_cost /= self.success_probability
        self.total_expansion_cost /= self.success_probability
        self.total_quality_control_cost /= self.success_probability
        self.total_harvesting_cost /= self.success_probability
        #print(self.total_pre_inoculation_cost, self.total_expansion_cost, self.total_quality_control_cost, self.total_harvesting_cost)

        self.total_failure_cost = round((self.total_consumables_cost + self.total_reagents_cost
                                   + self.total_facility_cost + self.total_labor_cost) - self.total_cost, 2)

        self.total_cost = round(self.total_cost + self.total_failure_cost, 2)
