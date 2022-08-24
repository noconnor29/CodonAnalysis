import csv
import tkinter as tk       # sudo apt install python3-tk
from tkinter.filedialog import askopenfilename
tk.Tk().withdraw()      # part of the import if not using other tkinter f(x)

filename = askopenfilename()
#inputfile ="\path\to\DNA_sequence_placeholder.csv" # use for testing
file = open(filename, "r")
seq = file.read()
