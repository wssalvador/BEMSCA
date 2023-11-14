# BEMSCA
BEMSCA - Bioprocess Economic Model for Stem Cell Applications

Stem Cell Engineering Research Group, Instituto Superior Tecnico, Universidade de Lisboa, Portugal

# Introduction
BEMSCA is a stochastic decision support tool for the optimal design of stem cell bioprocesses. It is capable of providing an optimized workflow and detailed economic evaluation of a user-designed bioprocess based on relevant inputs. For a more detailed explanation concerning the importance of bioprocess economic modeling, as well as the structure and potential applications of BEMSCA, please refer to the following research paper:

# Using BEMSCA
In order to use BEMSCA, a user must have a version of Python 3 installed on their computer. If a user does not have Python 3 installed, they may downloaded it from: https://www.python.org/downloads/.

Aside from a Python 3 source, the numpy, pandas and matplotlib packages should be installed, since these are essential for the proper functioning of BEMSCA's source code. The easiest way to do this is through Python's pip package installer, which should already come with the downloaded Python 3 source. A user should execute the following 3 commands in their command prompt (or terminal):

```
pip3 install numpy
pip3 install pandas
pip3 install matplotlib
```

Note: pip may have to be used instead of pip3, or whatever alias has been defined in the user's operating system.

When altering BEMSCA's source code to adapt it to specific production scenarios, the use of a code editor is recommended. Visual Studio Code is a good candidate for this purpose, but there are many other alternatives that can be employed according to the user's preference.

The folder provided in this GitHub repository includes all the necessary files for BEMSCA to function, along with a default SQLite database (database.db) and example graphs in the "results" subfolder.
