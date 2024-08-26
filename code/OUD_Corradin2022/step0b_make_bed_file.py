import pandas as pd

# read xlsx file
dflost = pd.read_excel("data/NIHMS1778910-supplement-Table_S4.xlsx", sheet_name="LostVELs")
print(dflost)
dfgain = pd.read_excel("data/NIHMS1778910-supplement-Table_S4.xlsx", sheet_name="GainVELs")
print(dfgain)
# output bed file
gain = "data/NIHMS1778910-supplement-Table_S4_gain.bed"
lost = "data/NIHMS1778910-supplement-Table_S4_lost.bed"
dflost.to_csv(lost, sep='\t', header=False, index=False)
dfgain.to_csv(gain, sep='\t', header=False, index=False)
