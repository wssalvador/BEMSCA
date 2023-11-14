# BEMSCA
BEMSCA - Bioprocess Economic Model for Stem Cell Applications

Stem Cell Engineering Research Group, Instituto Superior Tecnico, Universidade de Lisboa, Portugal

## Introduction
BEMSCA is a stochastic decision support tool for the optimal design of stem cell bioprocesses. It is capable of providing an optimized workflow and detailed economic evaluation of a user-designed bioprocess based on relevant inputs. With the source code of this repository, BEMSCA was applied to the optimization of the large-scale expansion of human induced pluripotent stem cells, as a demonstration of its functionality. For a more detailed explanation concerning the importance of bioprocess economic modeling, as well as the structure and potential applications of BEMSCA, please refer to the following research paper:

## Using BEMSCA
In order to use BEMSCA, the user must have a version of Python 3 installed on their computer. If the user does not have Python 3 installed, they may downloaded it from: https://www.python.org/downloads/.

Aside from a Python 3 source, the numpy, pandas and matplotlib packages should be installed, since these are essential for the proper functioning of BEMSCA's source code. The easiest way to do this is through Python's pip package installer, which should already come with the downloaded Python 3 source. The user should execute the following 3 commands in their terminal (or command prompt):

```
pip3 install numpy
pip3 install pandas
pip3 install matplotlib
```

Note: pip may have to be used instead of pip3, or whatever alias has been defined in the user's operating system.

The folder provided in this GitHub repository includes all the necessary files for BEMSCA to function, along with a default SQLite database (database.db) and example graphs in the "results" subfolder. All new graphs created by the user will also be saved to the "results" subfolder.

To run BEMSCA, the user should open their terminal (or command prompt) within the BEMSCA folder and execute the following command:

```
python3 main.py
```

Note: python may have to be used instead of python3, or whatever alias has been defined in the user's operating system.

The user should then follow the instructions presented to interact with BEMSCA. As an alternative, BEMSCA can be run using a code editor of the user's choice (e.g., Visual Studio Code).

The user is encouraged to alter BEMSCA's source code according to his specific production scenarios. If the user wishes to alter BEMSCA's database, they must first remove the existing database from the "BEMSCA" folder. They can then modify the database.py file according to their preferences, but must take care to respect the existing organization of the tables present in this file. The user can change values, create new table entries, or even create entirely new tables, but the user may need to execute additional modifications to the rest of BEMSCA's source code. When the user next runs BEMSCA, a new database.db file will be created reflecting the modifications to database.py.

Programming in Python is required to employ BEMSCA to its full potential, but by following the existing structures of the source code a basic level of understanding is sufficient.

Thank you for taking an interest in BEMSCA, we hope it may prove useful for your work. If you have any queries related to BEMSCA, fell free to send them to: william.salvador@tecnico.ulisboa.pt.
