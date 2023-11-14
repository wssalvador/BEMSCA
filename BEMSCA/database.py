import sqlite3


# -----------------------------------------------------------------------------
#    IMPORTANT NOTES
# -----------------------------------------------------------------------------
# All costs are in €

# Energy cost is in €/kWh

# All volumes are in L

# All lifespans are in years

# All energies are in kWh/day

# All costs in the Antibodies Table are in €/use

# All costs of the Reagents Table are in €/L (for solids the working concentration is taken into account)

# All areas in the Construction Costs Table are in m2 and its costs in €/m2

# All salaries in the Labor Costs Table and costs in the Operating Costs Table are in €/year

# The costs presented in the Quality Controls Table represent the cumulative cost of reagents exclusive to
# the respective quality assay, they DO NOT necessarily represent the total cost of that assay. The total
# cost is calculated during BEMSCA's workflow simulations (see "bioprocess.py" -> Quality Control Costs)

# Use factor denotes the percentage of a day during which an equipment is consuming energy (being used)


# -----------------------------------------------------------------------------
#    FUNCTIONS
# -----------------------------------------------------------------------------
# Function for initializing BEMSCA's default database (if database does not already exist)
def initialize_database():
    
    # Connect to database (creates database if it does not exist)
    connection = sqlite3.connect('database.db')

    # Create cursor
    cursor = connection.cursor()

    # ||---------------------[ PROGRAM VARIABLES TABLES ]--------------------||
    # Program Variables Tables Index (in order of appearence in code)
    # [Color Palettes]

    # <>----------------- Color Palletes ----------------<>    
    # Create Color Palletes Table    
    cursor.execute('''CREATE TABLE "Color Palettes" (
                    "Name" text,
                    "Primary Colors" text,
                    "Secondary Colors" text 
                   )''')
    
    # Insert default values into Color Palletes Table
    cursor.execute('''INSERT INTO "Color Palettes" VALUES
                    ("Gold/Red (2 media)", "goldenrod;firebrick", "gold;salmon"),
                    ("Gold/Brown (2 media)", "goldenrod;peru", "gold;peachpuff"),
                    ("Gold/Warm (2 media)", "rosybrown;goldenrod;indianred;lightcoral;lightsalmon",
                     "peachpuff;gold;peachpuff;peachpuff;peachpuff")
                   ''')
    # ||-----------------[ END OF PROGRAM VARIABLES TABLES ]-----------------||

    # ||--------------------[ ECONOMIC VARIABLES TABLES ]--------------------||
    # Economic Variables Tables Index (in order of appearence in code)
    # [2D Platforms, Antibodies, Bioreactors, Construction Costs, Equipment,
    #  Facility Specifications, Labor Costs, Operating Costs, Quality Controls, Reagents]

    # <>------------------ 2D Platforms -----------------<>    
    # Create 2D Platforms Table
    cursor.execute('''CREATE TABLE "2D Platforms" (
                    "Name" text,
                    "Surfaces" integer,
                    "Coating Volume" real,
                    "Washing Volume" real,
                    "Culture Volume" real,
                    "Surface Confluency" integer,
                    "Confluency STD" real,
                    "Cost" real
                   )''')

    # Insert default values into 2D Platforms Table
    cursor.execute('''INSERT INTO "2D Platforms" VALUES
                    ("12-well", 12, 0.0005, 0.0005, 0.001, 400000, 33333, 2.51),
                    ("6-well", 6, 0.001, 0.001, 0.002, 1200000, 100000, 2.40)
                   ''')

    # <>---------------- Antibodies Table ---------------<>    
    # Create Antibodies Table
    cursor.execute('''CREATE TABLE "Antibodies" (
                    "Name" text,
                    "Cost" real
                   )''')

    # Insert default values into Antibodies Table
    cursor.execute('''INSERT INTO "Antibodies" VALUES
                    ("OCT4", 3.85),
                    ("Secondary Antibody", 0.58),
                    ("SSEA-4", 0.79),
                    ("SSEA-4-PE", 0.58),
                    ("SOX2", 3.74),
                    ("TRA-1-60", 0.90),
                    ("TRA-1-60-PE", 0.58)
                   ''')
    
    # <>-------------- Bioreactors Table -------------<>
    # Create Bioreactors Table
    cursor.execute('''CREATE TABLE "Bioreactors" (
                    "Name" text,
                    "Min Volume" real,
                    "Max Volume" real,
                    "Acquisition Cost" real,
                    "Use Cost" real,
                    "Energy Consumption" real
                   )''')
    
    # Insert default values into Bioreactors Table
    cursor.execute('''INSERT INTO "Bioreactors" VALUES
                    ("PBS 0.1MAG", 0.060, 0.100, 2473.00, 261.75, 0.07),
                    ("PBS 0.5MAG", 0.300, 0.500, 2473.00, 327.75, 0.07),
                    ("PBS 3MAG", 1.800, 3.000, 65550.00, 1248.00, 8.64)
                   ''')

    # <>--------------- Construction Costs --------------<>    
    # Create Construction Costs Table
    cursor.execute('''CREATE TABLE "Construction Costs" (
                    "Name" text,
                    "Area" integer,
                    "Cost" real
                   )''')
    
    # Insert default values into Construction Costs Table
    cursor.execute('''INSERT INTO "Construction Costs" VALUES
                    ("Standard Rooms", 268, 3692.00),
                    ("Clean Rooms", 132, 6329.00)
                   ''')    

    # <>------------------- Equipment -------------------<>    
    # Create Equipment Table
    cursor.execute('''CREATE TABLE "Equipment" (
                    "Name" text,
                    "Amount" integer,
                    "Acquisition Cost" real,
                    "Energy Consumption" real,
                    "Use Factor" real
                   )''')

    # Insert default values into Equipment Table
    cursor.execute('''INSERT INTO "Equipment" VALUES
                    ("Autoclave", 2, 18330.00, 240.00, 0.50),
                    ("Biosafety Cabinet", 6, 10865.00, 8.16, 0.50),
                    ("Centrifuge", 6, 11995.00, 40.80, 0.10),
                    ("Cryo Freezer", 1, 33491.00, 64.80, 0.05),
                    ("Flow Cytometer", 1, 94700.00, 3.60, 0.20),
                    ("Fluorescence Microscope", 1, 12245.00, 2.40, 0.10),
                    ("Freezer (-20 °C)", 2, 9075.00, 82.80, 0.33),
                    ("Freezer (-80 °C)", 1, 16760.00, 33.67, 0.33),
                    ("Fridge", 4, 6500.00, 82.80, 0.33), 
                    ("Incubator", 12, 9604.00, 13.92, 0.50),
                    ("Microscope", 6, 2471.00, 0.72, 0.10),
                    ("RT-PCR System", 1, 20560.00, 10.61, 0.10),
                    ("Spectrophotometer", 1, 10169.39, 0.06, 0.10)
                   ''')
    
    # <>------------ Facility Specifications ------------<>    
    # Create Facility Specifications Table
    cursor.execute('''CREATE TABLE "Facility Specifications" (
                    "Name" text,
                    "Value" real
                   )''')

    # Insert default values into Facility Specifications Table
    cursor.execute('''INSERT INTO "Facility Specifications" VALUES
                    ("Energy Cost", 0.12266),
                    ("Equipment Lifespan", 5),
                    ("Facility Lifespan", 15),
                    ("Parallel Processes", 6)
                   ''')    
    
    # <>------------------ Labor Costs ------------------<>    
    # Create Labor Costs Table
    cursor.execute('''CREATE TABLE "Labor Costs" (
                    "Name" text,
                    "Number" integer,
                    "Salary" real
                   )''')
    
    # Insert default values into Labor Costs Table
    cursor.execute('''INSERT INTO "Labor Costs" VALUES
                    ("Lab Technician", 6, 30000),
                    ("Supervisors", 3, 50000)
                   ''')

    # <>---------------- Operating Costs ----------------<>    
    # Create Operating Costs Table
    cursor.execute('''CREATE TABLE "Operating Costs" (
                    "Name" text,
                    "Cost" real
                   )''')

    # Insert default values into Operating Costs Table
    cursor.execute('''INSERT INTO "Operating Costs" VALUES
                    ("Additional Supplies", 7900.00),
                    ("Cleaning", 28000.00),
                    ("Energy", 54461.00),
                    ("Garments", 2000.00), 
                    ("Gases", 21600.00),
                    ("Maintenance", 52800.00),  
                    ("Requalification", 65400.00)
                   ''')
    
    # <>------------- Quality Controls Table ------------<>    
    # Create Quality Controls Table
    cursor.execute('''CREATE TABLE "Quality Controls" (
                    "Name" text,
                    "Cost" real
                   )''')

    # Insert default values into Quality Controls Table
    cursor.execute('''INSERT INTO "Quality Controls" VALUES
                    ("Trilineage Differentiation", 85.00),
                    ("PCR Genomic Screening", 30.00),
                    ("Intracellular FC", 1.50),
                    ("Intracellular Immunocytochemistry", 1.10),
                    ("Karyotyping", 7.00),
                    ("RT-PCR", 22.50),
                    ("Surface FC", 0.10),
                    ("Surface Immunocytochemistry", 0.20)                  
                   ''')

    # <>----------------- Reagents Table ----------------<>    
    # Create Reagents Table
    cursor.execute('''CREATE TABLE "Reagents" (
                    "Name" text,
                    "Cost" real
                   )''')

    # Insert default values into Reagents Table
    cursor.execute('''INSERT INTO "Reagents" VALUES
                    ("Accutase", 468.00),
                    ("B8", 9.63),
                    ("DMEM/F-12", 54.18),
                    ("DPBS", 40.13),
                    ("Dextran Sulfate", 2.79),
                    ("E8", 584.00),
                    ("EDTA", 40.46),
                    ("Matrigel", 708.94),
                    ("mTeSR1", 654.00),
                    ("Y_27632", 139.01)
                   ''')
    # ||-----------------[ END OF ECONOMIC VARIABLES TABLES ]----------------||

    # ||-------------------[ BIOPROCESS VARIABLES TABLES ]-------------------||
    # Biologic Variables Tables Index (in order of appearence in code)
    # [Bioprocess Parameters, Expansion Simulations, Recovery Simulations]

    # <>---------- Bioprocess Parameters Table ----------<>    
    # Create Bioprocess Parameters Table
    cursor.execute('''CREATE TABLE "Bioprocess Parameters" (
                    "Name" text,
                    "Initial Cell Number" integer,
                    "Target Cell Number" integer,
                    "2D Platform" text,
                    "Coating Substrate" text,
                    "Bioreactors" text,
                    "Dissociation Enzyme" text,
                    "Expansion Simulation" text,
                    "Recovery Simulation" text,
                    "Minimum Threshold" real
                   )''')
    
    # Insert default values into Bioprocess Parameters Table
    cursor.execute('''INSERT INTO "Bioprocess Parameters" VALUES
                    ("Default", 1000000, 2000000000, "6-well", "Matrigel", "PBS 0.1MAG;PBS 0.5MAG;PBS 3MAG",
                     "Accutase", "Nogueira2019", "Borys2021", 0.95),
                    ("DS Supplementation", 1000000, 2000000000, "6-well", "Matrigel", "PBS 0.1MAG;PBS 0.5MAG;PBS 3MAG",
                     "Accutase", "Nogueira2019 (+DS)", "Borys2021", 0.95),
                    ("Low-Density Inoculation", 1000000, 2000000000, "6-well", "Matrigel", "PBS 0.1MAG;PBS 0.5MAG;PBS 3MAG",
                     "Accutase", "Borys2021", "Borys2021", 0.95),
                    ("B8 (0.75x)", 1000000, 2000000000, "6-well", "Matrigel", "PBS 0.1MAG;PBS 0.5MAG;PBS 3MAG",
                     "Accutase", "Nogueira2019 (B8-0.75x)", "Borys2021", 0.95),
                    ("B8 (1x)", 1000000, 2000000000, "6-well", "Matrigel", "PBS 0.1MAG;PBS 0.5MAG;PBS 3MAG",
                     "Accutase", "Nogueira2019 (B8-1x)", "Borys2021", 0.95),
                    ("B8 (1.5x)", 1000000, 2000000000, "6-well", "Matrigel", "PBS 0.1MAG;PBS 0.5MAG;PBS 3MAG",
                     "Accutase", "Nogueira2019 (B8-1.5x)", "Borys2021", 0.95),
                    ("B8 (2x)", 1000000, 2000000000, "6-well", "Matrigel", "PBS 0.1MAG;PBS 0.5MAG;PBS 3MAG",
                     "Accutase", "Nogueira2019 (B8-2x)", "Borys2021", 0.95)
                   ''')

    # <>---------- Expansion Simulations Table ----------<>    
    # Create Expansion Simulations Table
    cursor.execute('''CREATE TABLE "Expansion Simulations" (
                    "Name" text,
                    "Culture Medium" text,
                    "Supplements" text,
                    "Seeding Density" integer,
                    "Fold Expansion AVG" real,
                    "Fold Expansion SEM" real,
                    "Experiment Sample Size" integer,
                    "Bioreactor Culture Time" integer,
                    "Bioreactor Volumes Spent" real
                   )''')
    
    # Insert default values into Expansion Simulations Table
    cursor.execute('''INSERT INTO "Expansion Simulations" VALUES
                    ("Nogueira2019", "mTeSR1", "n/a", 250000, 4.8, 0.5, 3, 7, 5.0),
                    ("Nogueira2019 (B8-0.75x)", "B8", "n/a", 250000, 3.6, 0.375, 3, 8, 5.8),
                    ("Nogueira2019 (B8-1x)", "B8", "n/a", 250000, 4.8, 0.5, 3, 7, 5.0),
                    ("Nogueira2019 (B8-1.5x)", "B8", "n/a", 250000, 7.2, 0.75, 3, 6, 4.2),
                    ("Nogueira2019 (B8-2x)", "B8", "n/a", 250000, 9.6, 1.0, 3, 5, 3.4),
                    ("Nogueira2019 (+DS)", "mTeSR1", "Dextran Sulfate", 250000, 9.3, 0.6, 3, 5, 3.4),
                    ("Borys2021", "mTeSR1", "n/a", 20000, 32.3, 3.2, 4, 6, 1.5)
                   ''')
    
    # <>----------- Recovery Simulations Table ----------<>    
    # Create Recovery Simulations Table
    cursor.execute('''CREATE TABLE "Recovery Simulations" (
                    "Name" text,
                    "Recovery Efficiency AVG" real,
                    "Recovery Efficiency SEM" real,
                    "Experiment Sample Size" integer
                   )''')
    
    # Insert default values into Recovery Simulations Table
    cursor.execute('''INSERT INTO "Recovery Simulations" VALUES
                    ("Borys2021", 0.952, 0.040, 4)
                   ''')

    # ||----------------[ END OF BIOPROCESS VARIABLES TABLES ]---------------||

    # Commit database changes (save initialized database)
    connection.commit()

    # Close database
    connection.close()