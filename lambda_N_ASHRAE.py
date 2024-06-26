
import pandas as pd
import numpy as np

df = pd.read_csv('./input/ashrae_db2.01_cleaned.csv')

print(df.describe())

# temperature range, checked with database, may need to be changed in practice
temperature_bins = range(0, 50)

result_df = pd.DataFrame(columns=['Climate', 'Ta', 'TSV', 'PMV'])

# using climate Csb as an example
Climate_sample = 'Aw'
climate_data = df[df['Climate'] == Climate_sample]

"""
    bin data
"""

for temp_range in zip(temperature_bins[:-1], temperature_bins[1:]):
    temp_min, temp_max = temp_range
    temp_data = climate_data[(climate_data['Ta'] >= temp_min) & (climate_data['Ta'] < temp_max)]

    avg_tsv = temp_data['TSV'].mean()
    avg_pmv = temp_data['PMV'].mean()

    result_df = result_df.append({
        'Climate': Climate_sample,
        'Ta': temp_min,
        'TSV': avg_tsv,
        'PMV': avg_pmv
    }, ignore_index=True)

df_bin = result_df
print(df_bin)

"""
    calculate lambda: using climate Csb as an example
"""
from sklearn.metrics import mean_squared_error

climate_data = df_bin[df_bin['Climate'] == Climate_sample]
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

print('New lambda (PMV<0): ', Climate_sample,'', lambda_value_neg)
print('New lambda (PMV>0): ', Climate_sample,'', lambda_value_pos)


