import os
import database
import utils
import bioprocess
import outputs

# -----------------------------------------------------------------------------
#    FUNCTIONS
# -----------------------------------------------------------------------------
# Function for executing tasks related to "compare" command
def compare_command():
    # Ask user which sets of bioprocess parameters should be compared (accept only valid input)
    print('\nSelect bioprocess parameters from the list below for comparison study (type them in one at a time):\n')
    print(database_data["Bioprocess Parameters"].index.tolist())
    selection = []
    sets = input('\n>>> ')
    while sets.capitalize() != 'Done':
        if sets not in database_data["Bioprocess Parameters"].index.tolist():    
            print('\nInvalid input. Please select bioprocess parameters from the list above (case sensitive):')
            sets = input('\n>>> ')
        else:
            selection.append(sets)
            print(f'\nSelected bioprocess parameters: {selection}')
            print(' \nSelect more bioprocess parameters from the list above or type "Done":')
            sets = input('\n>>> ')

    # Ask user which set of bioprocess parameters should be used as the base-case and determine its repective index
    print('\nWhich of the selected bioprocess parameters should be used as the base-case?')
    base_case = input('\n>>> ')
    while base_case not in selection:
        print('\nInvalid input. Please select bioprocess parameters from those previously selected (case sensitive):')
        base_case = input('\n>>> ')

    base_case_index = selection.index(base_case)

    # Ask user what the name of the comparison study should be
    print('\nDefine name of comparison study:')
    study_name = input('\n>>> ')

    # Ask user to choose a color palett for graphical output
    print('\nChoose color palette from the list below to be applied to output graphs:\n')
    print(database_data["Color Palettes"].index.tolist())
    color_palette = input('\n>>> ')
    while color_palette not in database_data["Color Palettes"].index.tolist():
        print('\nInvalid input. Please select color palette from the list above (case sensitive):')
        color_palette = input('\n>>> ')

    # Ask user to define the label to be used for each set of bioprocess parameters when making graphs
    print('\nDefine graph labels for each set of bioprocess parameters, in order (type them in one at a time):')
    labels = []
    label = input('\n>>> ')
    labels.append(label)
    while True:
        print(f'\nDefined labels: {labels} -> {len(selection) - len(labels)} still left to define')
        label = input('\n>>> ')
        labels.append(label)
        if (len(selection) - len(labels) == 0):
            break    

    # Define output customization array
    output_customization = [base_case_index, study_name, color_palette, labels]

    # Simulate bioprocess using each set of bioprocess parameters
    simulation_results = []
    for sets in selection:
        simulation_results.append(bioprocess.Bioprocess(database_data, database_data["Bioprocess Parameters"].loc[[sets]]))

    # Present outputs of comparison study
    outputs.compare_output(database_data, selection, simulation_results, output_customization)


# Function for executing tasks related to "help" command
def help_command():
    print('\nTo simulate a specific set of bioprocess parameters type "Simulate".')
    print('To execute a comparison study with more than one set of bioprocess parameters type "Compare".')
    print('To terminate BEMSCA type "Quit".')


# Function for executing tasks related to "simulate" command
def simulate_command():
    # Ask user which set of bioprocess parameter he wishes to simulate (accept only valid input)
    print('\nSelect which set of bioprocess parameters should be simulated from the list below:\n')
    print(database_data["Bioprocess Parameters"].index.tolist())
    set = input('\n>>> ')
    while set not in database_data["Bioprocess Parameters"].index.tolist():
        print('\nInvalid input. Please select preset bioprocess parameters from the list below (case sensitive):')
        set = input('\n>>> ')

    # Simulate bioprocess using the selected preset bioprocess parameters
    simulation_results = bioprocess.Bioprocess(database_data, database_data["Bioprocess Parameters"].loc[[set]])

    # Present outputs of simulated bioprocess
    outputs.simulate_output(database_data, set, simulation_results)


# -----------------------------------------------------------------------------
#    INPUT LOOP
# -----------------------------------------------------------------------------
# Check if database exists and is accessible, if it does not initialize it
if not os.path.exists("database.db"):
    database.initialize_database()

# Load stored database data so BEMSCA has access to it (stored in pandas dataframes)
database_data = utils.get_database_data()

# Initialize variable to store user command
command = ''

# Print welcome message
print('\n<>------------------------------------------------------------------------------------------------------<>')
print('  Welcome to BEMSCA. For a list of available commands, type "help". Alternatively, type desired command.  ')
print('<>------------------------------------------------------------------------------------------------------<>')

# Start input loop (terminate when 'quit' command is received)
while command != 'Quit':

    command = input('\n>>> ').capitalize()

    if command == 'Compare':
        compare_command()
    elif command == 'Help':
        help_command()
    elif command == 'Simulate':
        simulate_command()
    elif command != 'Quit':
        print('\nInvalid command.')