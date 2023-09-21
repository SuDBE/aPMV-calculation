
import pandas as pd
import numpy as np

df = pd.read_csv('./aPMV.csv')

print(df.describe())

"""
    it was copied from: pythermalcomfort.models import pmv_ppd_optimized
"""
import math

def pmv_ppd_optimized(tdb, tr, vr, rh, met, clo, wme):

    pa = rh * 10 * math.exp(16.6536 - 4030.183 / (tdb + 235))

    icl = 0.155 * clo  # thermal insulation of the clothing in M2K/W
    m = met * 58.15  # metabolic rate in W/M2
    w = wme * 58.15  # external work in W/M2
    mw = m - w  # internal heat production in the human body
    # calculation of the clothing area factor
    if icl <= 0.078:
        f_cl = 1 + (1.29 * icl)  # ratio of surface clothed body over nude body
    else:
        f_cl = 1.05 + (0.645 * icl)

    # heat transfer coefficient by forced convection
    hcf = 12.1 * math.sqrt(vr)
    hc = hcf  # initialize variable
    taa = tdb + 273
    tra = tr + 273
    t_cla = taa + (35.5 - tdb) / (3.5 * icl + 0.1)

    p1 = icl * f_cl
    p2 = p1 * 3.96
    p3 = p1 * 100
    p4 = p1 * taa
    p5 = (308.7 - 0.028 * mw) + (p2 * (tra / 100.0) ** 4)
    xn = t_cla / 100
    xf = t_cla / 50
    eps = 0.00015

    n = 0
    while abs(xn - xf) > eps:
        xf = (xf + xn) / 2
        hcn = 2.38 * abs(100.0 * xf - taa) ** 0.25
        if hcf > hcn:
            hc = hcf
        else:
            hc = hcn
        xn = (p5 + p4 * hc - p2 * xf ** 4) / (100 + p3 * hc)
        n += 1
        if n > 150:
            raise StopIteration("Max iterations exceeded")

    tcl = 100 * xn - 273

    # heat loss diff. through skin
    hl1 = 3.05 * 0.001 * (5733 - (6.99 * mw) - pa)
    # heat loss by sweating
    if mw > 58.15:
        hl2 = 0.42 * (mw - 58.15)
    else:
        hl2 = 0
    # latent respiration heat loss
    hl3 = 1.7 * 0.00001 * m * (5867 - pa)
    # dry respiration heat loss
    hl4 = 0.0014 * m * (34 - tdb)
    # heat loss by radiation
    hl5 = 3.96 * f_cl * (xn ** 4 - (tra / 100.0) ** 4)
    # heat loss by convection
    hl6 = f_cl * hc * (tcl - tdb)

    ts = 0.303 * math.exp(-0.036 * m) + 0.028
    _pmv = ts * (mw - hl1 - hl2 - hl3 - hl4 - hl5 - hl6)
    _ppd = 100.0 - 95.0 * math.exp(-0.03353 * pow(_pmv, 4.0) - 0.2179 * pow(_pmv, 2.0))

    return {"pmv": round(_pmv, 2), "ppd": round(_ppd, 2)}
    # return _pmv

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

print(df)

# temperature range, checked with database, may need to be changed in practice
temperature_bins = range(0, 50)

result_df = pd.DataFrame(columns=['Ta', 'TSV', 'PMV'])

climate_data = df

"""
    bin data
"""

for temp_range in zip(temperature_bins[:-1], temperature_bins[1:]):
    temp_min, temp_max = temp_range
    temp_data = climate_data[(climate_data['Ta'] >= temp_min) & (climate_data['Ta'] < temp_max)]

    avg_tsv = temp_data['TSV'].mean()
    avg_pmv = temp_data['PMV'].mean()

    result_df = result_df.append({
        'Ta': temp_min,
        'TSV': avg_tsv,
        'PMV': avg_pmv
    }, ignore_index=True)

df_bin = result_df
print(df_bin)

"""
    calculate lambda
"""
from sklearn.metrics import mean_squared_error

climate_data = df_bin
y_true = []
y_hat1 = []
y_hat2 = []

# PMV < 0
sum_div = 0
sum_multi = 0
lambda_value_neg = None
for index, row in climate_data.iterrows():
    if row['PMV'] < 0:
        div = row['PMV'] - row['TSV']
        multi = row['PMV'] * row['TSV']
        sum_div += div
        sum_multi += multi
        y_true.append(row['TSV'])
        y_hat1.append(row['PMV'] / (1 + 0.66 * row['PMV']))
        y_hat2.append(row['PMV'] / (1 - 0.66 * row['PMV']))

        error1 = mean_squared_error(y_true, y_hat1, squared=False)
        error2 = mean_squared_error(y_true, y_hat2, squared=False)

if sum_multi != 0:
    lambda_value_neg = sum_div / sum_multi
    if abs(lambda_value_neg) > 0.66:
        if error1 < error2:
            lambda_value_neg = 0.66
        else:
            lambda_value_neg = -0.66

# PMV > 0
sum_div = 0
sum_multi = 0
lambda_value_pos = None
for index, row in climate_data.iterrows():
    if row['PMV'] >= 0:
        div = row['PMV'] - row['TSV']
        multi = row['PMV'] * row['TSV']
        sum_div += div
        sum_multi += multi
        y_true.append(row['TSV'])
        y_hat1.append(row['PMV'] / (1 + 0.66 * row['PMV']))
        y_hat2.append(row['PMV'] / (1 - 0.66 * row['PMV']))

        error1 = mean_squared_error(y_true, y_hat1, squared=False)
        error2 = mean_squared_error(y_true, y_hat2, squared=False)

if sum_multi != 0:
    lambda_value_pos = sum_div / sum_multi

    if abs(lambda_value_pos) > 0.66:
        if error1 < error2:
            lambda_value_pos = 0.66
        else:
            lambda_value_pos = -0.66

if lambda_value_neg is not None:
    lambda_value_neg = "{:.2f}".format(lambda_value_neg)
else:
    lambda_value_neg = 'NaN'

if lambda_value_pos is not None:
    lambda_value_pos = "{:.2f}".format(lambda_value_pos)
else:
    lambda_value_pos = 'NaN'

print('λ (PMV<0): ', lambda_value_neg)
print('λ (PMV>0): ', lambda_value_pos)

"""
    Draw aPMV curve
"""
x_neg_P = np.arange(-3, 0, 0.01)
y_neg_P = list()
if lambda_value_neg != 'NaN':
    lambda_value_neg = float(lambda_value_neg)
    for t in x_neg_P:
        yi = t / (1 + lambda_value_neg * t)  # y轴具体数值
        y_neg_P.append(yi)

if lambda_value_pos != 'NaN':
    x_pos_P = np.arange(0, 3, 0.01)
    y_pos_P = list()
    lambda_value_pos = float(lambda_value_pos)
    for t in x_pos_P:
        yi = t / (1 + lambda_value_pos * t)  # y轴具体数值
        y_pos_P.append(yi)

import time
if lambda_value_pos == 'NaN' and lambda_value_neg == 'NaN':
    print("Please use valid data!")
    time.sleep(5)

import matplotlib.pyplot as plt
import seaborn as sns

climate_data['PMV category'] = ['PMV > 0' if pmv >= 0 else 'PMV < 0' for pmv in climate_data['PMV']]

plt.figure(figsize=(6, 5))

alpha = 0.85
width1 = 2.1

if lambda_value_neg != 'NaN':
    label_neg = f'PMV<0, λ={lambda_value_neg:.2f}'

if lambda_value_pos != 'NaN':
    label_pos = f'PMV>0, λ={lambda_value_pos:.2f}'

if lambda_value_neg != 'NaN':
    sns.lineplot(x=x_neg_P, y=y_neg_P, label=label_neg, color='royalblue',alpha = alpha, linestyle='dashed', linewidth=width1)
if lambda_value_pos != 'NaN':
    sns.lineplot(x=x_pos_P, y=y_pos_P, label=label_pos, color='orangered',alpha =alpha, linestyle='dashed', linewidth=width1)

sns.scatterplot(x='PMV', y='TSV', data=climate_data, hue='PMV category', legend=False)

plt.grid(axis='y', linestyle='--', alpha=0.4)

plt.xlim(-3, 3)
plt.ylim(-3, 3)
plt.legend(loc='upper left')
plt.ylabel('TSV/aPMV')
plt.title('aPMV curve')
plt.savefig('pic-aPMV curve.png', format='png', dpi=400)
# plt.show()

def calculate_pred(row):
    if row['PMV category'] == 'PMV < 0' and lambda_value_neg != 'NaN':
        return row['PMV'] / (1 + lambda_value_neg * row['PMV'])
    elif row['PMV category'] == 'PMV > 0' and lambda_value_pos != 'NaN':
        return row['PMV'] / (1 + lambda_value_pos * row['PMV'])
    else:
        return None

climate_data['aPMV'] = climate_data.apply(calculate_pred, axis=1)

print(climate_data)

fig, ax = plt.subplots(figsize=(8, 5))

df = climate_data

ax.scatter(df['Ta'], df['TSV'], label='TSV', marker='o', color='black', facecolors='none', edgecolors='black')
ax.scatter(df['Ta'], df['PMV'], label='PMV', marker='s', color='royalblue', facecolors='none', edgecolors='royalblue')
ax.scatter(df['Ta'], df['aPMV'], label='λ_N', marker='^', color='forestgreen')

ax.axhline(y=0, color='gray', linestyle='--', label=None)
plt.grid(axis='y', linestyle='--', alpha=0.4)

ax.set_title('aPMV prediction')
ax.set_xlabel('Air temperature, °C')
ax.set_ylabel('TSV/Prediction')
ax.set_ylim(-3.5, 3.5)

ax.legend()
plt.legend(loc='upper left')
plt.savefig('pic-aPMV points.png', format='png', dpi=400)
# plt.show()