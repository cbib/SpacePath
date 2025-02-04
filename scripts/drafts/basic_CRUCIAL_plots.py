import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

"""
completely unrelated note of Bash use: 
change a substring on filename ("" to make it work even if spaces within)
then copy to new filename
for i in *Batch*;do echo $i; loo="${i//Batch/Grps}"; cp "$i" "$loo" ; done
finally I delete the old files
"""

pixels_number = [11677, 16943, 17334, 10433, 13622,
                  19714, 28992, 4947, 92476, 16875,
                 103874, 133444 ]

total_mz = [857, 857, 857, 857, 857, 857,
      1060, 1110, 2279, 1179, 3622, 2930]

annotated_mz = [118, 118, 118, 118, 118, 118,
          70, 42, 109, 62, 160, 163]

kegg_compounds = [208, 208,208,208,208,208,
       103, 71, 162, 93, 232, 340]

samples_group = ['Ravi', 'Ravi', 'Ravi', 'Ravi', 'Ravi', 'Ravi',
                 'UB', 'UB', 'UB', 'UB', 'UB', 'Ellis']

df = pd.DataFrame({
    "pixels_number":  pixels_number,
    "total_mz": total_mz,
    "annotated_mz":  annotated_mz,
    "kegg_compounds": kegg_compounds ,
    "samples_group": samples_group
})


sns.scatterplot(data=df,
                x="pixels_number", y="total_mz",
                hue="annotated_mz", size="annotated_mz")
plt.show()


g = sns.scatterplot(data=df,
                x="pixels_number", y="total_mz",
                hue="samples_group", s=100, alpha=0.5,
                )
sns.move_legend(g, "upper left")
plt.show()


sns.scatterplot(data=df,
                x="pixels_number", y="annotated_mz",
                hue="kegg_compounds")
plt.show()


g = sns.scatterplot(data=df,
                x="pixels_number", y="annotated_mz",
                hue="samples_group", s=100, alpha=0.5,
                )
sns.move_legend(g, "upper left")
plt.show()

