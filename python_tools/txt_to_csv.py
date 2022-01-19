import numpy as np
import csv
import pandas as pd

with open("results.csv", "w", newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter='\t')

    f = open("results.txt", "r")
    lines = f.readlines()
    print(lines[0])
    # count = 0
    for line in lines:
        # count += 1

        if line != lines[0]:
            elements = line.split(' ')
            elements = [x for x in elements if x]
            writer.writerow(elements[:5])