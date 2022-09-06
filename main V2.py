import pandas as pd
import numpy as np

# read data .csv files
cell_data = pd.read_csv('ParameterData_Main.csv')
pex_data = pd.read_csv('ParameterData_peroxisome.csv')

# remove calls out of gate
cell_data = cell_data[cell_data['R01'] == 1].reset_index(drop=True)
pex_data = pex_data[pex_data['MO.R01'] == 1].reset_index(drop=True)

# add strain column
pex_data['Strain'] = pex_data['Well '].replace({7: 'Control', 16: 'MutD', 25: 'MutS'})
cell_data['Strain'] = cell_data['Well '].replace({7: 'Control', 16: 'MutD', 25: 'MutS'})

# create new dataframe - for each cells, find corresponding peroxisomes
new_df = cell_data[['Object ID', 'Well ', 'Strain', 'Mean Intensity GFP',
                    'Total Intensity GFP', 'Area ']].copy(deep=True)
new_df = new_df.rename(columns={'Object ID': 'Cell ID', 'Mean Intensity GFP': 'Cell Mean Intensity GFP',
                                'Total Intensity GFP': 'Cell Total Intensity GFP', 'Area ': 'Cell Area'})
# add empty (nan values) columns (to be filled later)
new_df['Pex Mean Intensity GFP'] = np.nan
new_df['Pex Total Intensity GFP'] = np.nan
new_df['Pex Area'] = np.nan
new_df = new_df.reset_index(drop=True)
# replace nan values with actual values - mean intensity and area of the corresponding peroxisome(s):
for idx in range(len(new_df)):
    cell_ID = new_df['Cell ID'][idx]
    cur_pex_data = pex_data.loc[pex_data['Parent Object ID (MO)'] == cell_ID]
    if cur_pex_data.empty:
        new_df.drop(idx, inplace=True)
    else:
        new_df.loc[idx, 'Pex Mean Intensity GFP'] = cur_pex_data['Mean Intensity GFP'].mean()
        new_df.loc[idx, 'Pex Total Intensity GFP'] = cur_pex_data['Total Intensity GFP'].mean()
        new_df.loc[idx, 'Pex Area'] = cur_pex_data['Area '].mean()

# save new_df to csv file
new_df.to_csv('organized data (V2).csv', index=False)

# sample 207 objects from each well
sampled_well_7 = new_df[new_df['Well '] == 7].sample(n=207)
sampled_well_16 = new_df[new_df['Well '] == 16].sample(n=207)
sampled_well_25 = new_df[new_df['Well '] == 25].sample(n=207)
# concatenate sampled data to one dataframe
sampled_data = pd.concat([sampled_well_7, sampled_well_16, sampled_well_25], ignore_index=True)

# save sampled_data to csv file
sampled_data.to_csv('sampled data (V2).csv', index=False)
