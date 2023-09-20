
import pandas as pd
import numpy as np

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 100)
pd.set_option('display.width', 1000)

df_ori = pd.read_excel('./input/ashrae_db2.01.xlsx')

columns_to_keep = ['Koppen climate classification', 'Cooling startegy_building level',
                   'Thermal sensation','Clo','Met','Air temperature ()','Radiant temperature ()','Relative humidity (%)',
                   'Air velocity (m/s)']

df = df_ori[columns_to_keep]

df.rename(columns={'Thermal sensation': 'TSV',
                   'Air temperature ()': 'Ta',
                   'Radiant temperature ()': 'Tr',
                   'Relative humidity (%)':'RH',
                   'Air velocity (m/s)':'Vel',
                   'Koppen climate classification':'Climate',
                   'Cooling startegy_building level': 'Operation'}, inplace=True)

df = df.dropna(subset=['Ta','Clo','RH','Vel','Met','TSV'])

df = df[df['Operation'] == 'Naturally Ventilated']

# upper and lower limits
upper_fences = df.quantile(0.75) + 1.5 * (df.quantile(0.75) - df.quantile(0.25))
lower_fences = df.quantile(0.25) - 1.5 * (df.quantile(0.75) - df.quantile(0.25))

print('upper_fences:\n',upper_fences)
print('lower_fences:\n',lower_fences)

# replace empty Tr with Ta
df['Tr'] = df['Tr'].fillna(df['Ta'])

# remove outliers
for column in ['Ta','Clo','RH','Vel','Met','Tr']:
    if pd.api.types.is_numeric_dtype(df[column]):
        df = df[(df[column] <= upper_fences[column]) & (df[column] >= lower_fences[column])]

print()

from pythermalcomfort.models import pmv_ppd_optimized

ta = df["Ta"].values
tr = df["Tr"].values
vel = df["Vel"].values
rh = df["RH"].values
met = df["Met"].values
clo = df["Clo"].values
df["PMV"] = None
df["PPD"] = None

results = np.vectorize(pmv_ppd_optimized)(ta, tr, vel, rh, met, clo, 0)

for index, row in df.iterrows():
    results = pmv_ppd_optimized(
        tdb=row["Ta"],
        tr=row["Tr"],
        vr=row["Vel"],
        rh=row["RH"],
        met=row["Met"],
        clo=row["Clo"],
        wme=0,
    )
    df.loc[index, "PMV"] = results["pmv"]
    df.loc[index, "PPD"] = results["ppd"]

df = df[(df['PMV'] >= -3) & (df['PMV'] <= 3)]

df.to_csv('./input/ashrae_db2.01_cleaned.csv')

print(df.describe())


