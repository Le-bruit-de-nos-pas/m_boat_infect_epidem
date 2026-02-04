import pandas as pd
import numpy as np

from scipy.stats import fisher_exact, mannwhitneyu

import shap
import xgboost
import inspect

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import sklearn

import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv('../data/neg_bacilli_turin_2024.txt', delimiter='\t')
print(df)
for name in df.columns:
    print(name)



print("Total missing # :", df.isna().sum().sum())
print("Total missing % :", 100* df.isna().sum().sum() / (df.shape[0]*df.shape[1]) )
df.fillna(0, inplace=True)
print("Total missing # :", df.isna().sum().sum())
print("Total missing % :", 100* df.isna().sum().sum() / (df.shape[0]*df.shape[1]) )



def calculate_se_proportion(p, n):
    return np.sqrt((p * (1 - p)) / n).round(2)



def create_summary_table(df):
    """ Generate summary table.
    Args:
    df (DataFrame): The DataFrame for which to generate the summary table.
    Returns:
    DataFrame: The summary table."""

    # Generate summary statistics for ALL columns
    summary_table = df.describe()
    # Transpose the summary table
    summary_table = summary_table.T
    # Round all values to two decimal points
    summary_table = summary_table.round(2)
    # Create 'Mean ± Std' column
    summary_table['Mean ± Std'] = summary_table.apply(lambda row: f"{row['mean']} ± {row['std']}", axis=1)
    # Create 'Median [IQR]' column
    summary_table['Median [IQR]'] = summary_table.apply(lambda row: f"{row['50%']} [{row['25%']}-{row['75%']}]", axis=1)
    # Create 'Proportion ± SE' column
    summary_table['Proportion ± SE'] = summary_table.apply(lambda row: f"{row['mean']} ± {calculate_se_proportion(row['mean'], row['count'])}", axis=1)
    # Create 'SE_Error_Flag' column
    summary_table['SE_Error_Flag'] = summary_table['Proportion ± SE'].apply(lambda row: 1 if 'nan' in str(row) else 0)
    # Reset the index
    summary_table.reset_index(inplace=True)
    # Create 'Row_Number' column
    summary_table['Row_Number'] = summary_table.index + 1
    # Drop unnecessary columns
    summary_table.drop(['count', 'mean', 'std', 'min', '25%', '50%', '75%', 'max'], axis=1, inplace=True)

    return summary_table


import pandas as pd

def calculate_sum_and_percentage(df):
    """ Calculate sum and percentage for each column.
    Args:
    df (DataFrame): The DataFrame for which to calculate the sum and percentage.
    Returns:
    DataFrame: A DataFrame with two columns: 'Column' for column names and 'Sum/Percentage' for their corresponding sum or percentage.
    """

    result = []
    
    total_rows = len(df)
    
    for column in df.columns:
        column_sum = df[column].sum()
        percentage = (column_sum / total_rows) * 100
        result.append([column, f"{column_sum} ({percentage:.2f}%)"])

    return pd.DataFrame(result, columns=['Column', 'Sum/Percentage'])


summary_table = create_summary_table(df)
summary_table.to_csv('../data/summary_table.csv')
summary_table

summary_table = calculate_sum_and_percentage(df)
summary_table.to_csv('../data/summary_table.csv')
summary_table





df_cult_pos = df[df['Positive_2culture_50perc'] == 1] # 43 out of 63
summary_table = create_summary_table(df_cult_pos)
summary_table.to_csv('../data/summary_table_cult_pos.csv')

df_cult_pos = df[df['Positive_2culture_50perc'] == 1] # 43 out of 63
summary_table = calculate_sum_and_percentage(df_cult_pos)
summary_table.to_csv('../data/summary_table_cult_pos.csv')

df_cult_neg = df[df['Positive_2culture_50perc'] == 0] # 20 out of 63
summary_table = create_summary_table(df_cult_neg)
summary_table.to_csv('../data/summary_table_cult_neg.csv')

df_cult_neg = df[df['Positive_2culture_50perc'] == 0] # 20 out of 63
summary_table = calculate_sum_and_percentage(df_cult_neg)
summary_table.to_csv('../data/summary_table_cult_neg.csv')

exclude_cols = ['Positive_2culture_50perc', 'pat_id']
binary_results = []
numeric_results = []

for col in df.columns:
    if col not in exclude_cols:
        if df[col].nunique() <= 2:
            filtered_col = df[col].dropna().astype(int)
            if filtered_col.nunique() == 2:
                contingency_table = pd.crosstab(df['Positive_2culture_50perc'], filtered_col)
                odds_ratio , p_value = fisher_exact(contingency_table)
                binary_results.append((col, odds_ratio, p_value))
            else:
                print(f"Warning: Colunm '{col}' does not exactly 2 unique values. Dropping it..." )
        else:
            group_1 = df[df['Positive_2culture_50perc'] == 1][col]
            group_2 = df[df['Positive_2culture_50perc'] == 0][col]
            _ , p_value = mannwhitneyu(group_1, group_2, alternative='two-sided')
            numeric_results.append((col, p_value))



print("Numeric Variables (Wilcoxon/Mann-Whitney Test): ")
for result in numeric_results:
    print(f"{result[0]} , p-value: {result[1]}" )

print("Binary Variables (Fisher's Exact Test): ")
for result in binary_results:
    print(f"{result[0]} - Odds Ratio: {result[1]}, p-value: {result[2]}" )



df_hosp_mort_pos = df[df['Mortality_hospital'] == 1] # 10 out of 63
summary_table = create_summary_table(df_hosp_mort_pos)
summary_table.to_csv('../data/summary_table_hosp_mort_pos.csv')

df_hosp_mort_pos = df[df['Mortality_hospital'] == 1] # 10 out of 63
summary_table = calculate_sum_and_percentage(df_hosp_mort_pos)
summary_table.to_csv('../data/summary_table_hosp_mort_pos.csv')

df_hosp_mort_neg = df[df['Mortality_hospital'] == 0] # 53 out of 63
summary_table = create_summary_table(df_hosp_mort_neg)
summary_table.to_csv('../data/summary_table_hosp_mort_neg.csv')

df_hosp_mort_neg = df[df['Mortality_hospital'] == 0] # 53 out of 63
summary_table = calculate_sum_and_percentage(df_hosp_mort_neg)
summary_table.to_csv('../data/summary_table_hosp_mort_neg.csv')


exclude_cols = ['Mortality_hospital', 'pat_id', 'Amoxicillin_clavulanate_including', 'Metronidazole_including', 'Mortality_30_days', 'Mortality_14_days']
binary_results = []
numeric_results = []

for col in df.columns:
    if col not in exclude_cols:
        if df[col].nunique() <= 2:
            filtered_col = df[col].dropna().astype(int)
            if filtered_col.nunique() == 2:
                contingency_table = pd.crosstab(df['Mortality_hospital'], filtered_col)
                odds_ratio , p_value = fisher_exact(contingency_table)
                binary_results.append((col, odds_ratio, p_value))
            else:
                print(f"Warning: Colunm '{col}' does not exactly 2 unique values. Dropping it..." )
        else:
            group_1 = df[df['Mortality_hospital'] == 1][col]
            group_2 = df[df['Mortality_hospital'] == 0][col]
            _ , p_value = mannwhitneyu(group_1, group_2, alternative='two-sided')
            numeric_results.append((col, p_value))

print("Numeric Variables (Wilcoxon/Mann-Whitney Test): ")
for result in numeric_results:
    print(f"{result[0]} , p-value: {result[1]}" )

print("Binary Variables (Fisher's Exact Test): ")
for result in binary_results:
    print(f"{result[0]} - Odds Ratio: {result[1]}, p-value: {result[2]}" )

df = pd.read_csv('../data/neg_bacilli_turin_2024.txt', delimiter='\t')
df.fillna(0, inplace=True)


# columns_to_drop = ['pat_id', 'Amoxicillin_clavulanate_including', 'Metronidazole_including', 'Mortality_30_days', 'Mortality_14_days']
columns_to_drop = ['pat_id', 'Amoxicillin_clavulanate_including', 'Metronidazole_including', 'Mortality_30_days', 'Mortality_14_days', 'LOS_days']
df = df.drop(columns=columns_to_drop)

y = df['Mortality_hospital']
X = df.drop(columns=['Mortality_hospital'])


X100 = shap.utils.sample(X, 100)

model = xgboost.XGBClassifier(nestimators=100, max_depth=2, objective = "binary:logistic").fit(X, y)
explainer = shap.Explainer(model, X)
shap_values = explainer(X)

shap_values.display_data = X.values

shap.initjs()

shap.force_plot(shap_values[10])

shap.force_plot(shap_values[11])

shap.force_plot(shap_values[14])

shap.force_plot(shap_values[41])

shap.force_plot(shap_values[42])

shap.force_plot(shap_values[36])

shap.plots.bar(shap_values)

shap.plots.bar(shap_values.abs.max(0))

shap.plots.beeswarm(shap_values)

shap.plots.beeswarm(shap_values.abs, color="shap_red")

shap.plots.heatmap(shap_values, max_display=10, plot_width=8)

df = pd.read_csv('../data/neg_bacilli_turin_2024.txt', delimiter='\t')
df.fillna(0, inplace=True)


columns_to_drop = ['pat_id', 'Amoxicillin_clavulanate_including', 'Metronidazole_including', 'Mortality_14_days', 'Mortality_30_days', 'LOS_days']
df = df.drop(columns=columns_to_drop)

numerical_features = df.select_dtypes(include=['float64', 'int64'])
numerical_features = numerical_features.drop(columns=['Mortality_hospital'])

binary_numerical_features = [col for col in numerical_features.columns if numerical_features[col].nunique() == 2]

non_binary_numerical_features = numerical_features.drop(columns=binary_numerical_features)

scaler = StandardScaler()
scaled_numerical_features = scaler.fit_transform(numerical_features)

df[numerical_features.columns] = scaled_numerical_features

numerical_features = df.drop(columns=['Mortality_hospital']).select_dtypes(include=['float64', 'int64'])
numerical_features.shape

pca = PCA(n_components=0.90) 
pca.fit(numerical_features)


principal_components = pca.transform(numerical_features)
principal_components_df = pd.DataFrame(data=principal_components, 
                                       columns=[f'PC{i}' for i in range(1, pca.n_components_ + 1)])


loadings = pca.components_
loadings_df = pd.DataFrame(loadings, columns=numerical_features.columns)


for i, component in enumerate(loadings_df.iterrows(), start=1):
    print(f"Principal Component {i} Loadings:")
    print(component)
    print("\n")

loadings_df.index = loadings_df.index + 1

plt.figure(figsize=(9, 16))
sns.heatmap(loadings_df.T, cmap='RdGy', annot=False, fmt=".1f", cbar=True)
plt.title('Principal Component Loadings')
plt.xlabel(None)  
plt.ylabel(None)  
plt.show()


cumulative_variance_ratio = np.cumsum(pca.explained_variance_ratio_)


plt.figure(figsize=(10, 6))
plt.plot(range(len(cumulative_variance_ratio)), cumulative_variance_ratio, marker='o', linestyle='-', color='black', markerfacecolor='none')
plt.title('Cumulative Explained Variance Ratio')
plt.xlabel('Number of Principal Components')
plt.ylabel('Cumulative Explained Variance Ratio')
plt.xticks(range(len(cumulative_variance_ratio)), range(1, len(cumulative_variance_ratio) + 1))
plt.yticks(np.arange(0, 1.1, 0.1))  # Start y-axis from 0 and include increments of 0.1
plt.grid(False)
plt.show()


y = df['Mortality_hospital']

correlation_with_y = principal_components_df.apply(lambda col: col.corr(y))

correlation_with_y_sorted = correlation_with_y.abs().sort_values(ascending=False)

print("Correlation with Output Variable (Absolute Values):")
print(correlation_with_y_sorted)

pca_with_y = pd.concat([principal_components_df, y], axis=1)

mean_loadings_by_y = pca_with_y.groupby('Mortality_hospital').mean()

print("Mean Loadings by Output Variable Group:")
print(mean_loadings_by_y)

mean_loadings_by_y.index = mean_loadings_by_y.index.map({0: 'No', 1: 'Yes'})


plt.figure(figsize=(25, 1))
sns.heatmap(mean_loadings_by_y, cmap='RdGy', annot=True, fmt=".1f", cbar=True)
plt.title('Mean Loadings by In-hospital Mortality Status')
plt.xlabel('Principal Components')
plt.ylabel('In-hospital Mortality')
plt.show()

top_features_by_pc = {}

for pc in loadings_df.index:
    top_features_by_pc[pc] = loadings_df.loc[pc].abs().nlargest(5).index.tolist()

for pc, features in top_features_by_pc.items():
    print(f"Principal Component {pc}: {features}")
