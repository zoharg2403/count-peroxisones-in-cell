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

# create new dataframe - for each pex mean intensity, find the parent cell mean intensity
new_df = pex_data[['Object ID', 'Well ', 'Strain', 'Mean Intensity GFP', 'Total Intensity GFP', 'Area ',
                   'Parent Object ID (MO)']].copy(deep=True)
new_df = new_df.rename(columns={'Object ID': 'Object ID (pex)', 'Mean Intensity GFP': 'Pex Mean Intensity GFP',
                                'Total Intensity GFP': 'Pex Total Intensity GFP', 'Area ': 'Pex Area'})
# add empty (nan values) columns (to be filled later)
new_df['Cell Mean Intensity GFP (MO)'] = np.nan
new_df['Cell Total Intensity GFP (MO)'] = np.nan
new_df['Cell Area (MO)'] = np.nan
new_df = new_df.reset_index(drop=True)
# replace nan values with actual values - Mean GFP intensity of the corresponding cell:
# for each pex 'MO_ID' is the Parent Object ID (MO = Main object) (ID of the corresponding cell)
for idx in range(len(new_df)):
    MO_ID = new_df['Parent Object ID (MO)'][idx]
    cur_MO_data = cell_data.loc[cell_data['Object ID'] == MO_ID]
    new_df.loc[idx, 'Cell Mean Intensity GFP (MO)'] = cur_MO_data['Mean Intensity GFP'].values[0]
    new_df.loc[idx, 'Cell Total Intensity GFP (MO)'] = cur_MO_data['Total Intensity GFP'].values[0]
    new_df.loc[idx, 'Cell Area (MO)'] = cur_MO_data['Area '].values[0]

# save new_df to csv file
new_df.to_csv('organized data.csv', index=False)

# sample 271 objects from each well
sampled_well_7 = new_df[new_df['Well '] == 7].sample(n=271, replace=True)
sampled_well_16 = new_df[new_df['Well '] == 16].sample(n=271, replace=True)
sampled_well_25 = new_df[new_df['Well '] == 25].sample(n=271, replace=True)
# concatenate sampled data to one dataframe
sampled_data = pd.concat([sampled_well_7, sampled_well_16, sampled_well_25], ignore_index=True)

# save sampled_data to csv file
sampled_data.to_csv('sampled data.csv', index=False)
