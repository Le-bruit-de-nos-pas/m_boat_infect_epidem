# Paulo Bastos   May 3rd 2024 

```python
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
````


```python
df = pd.read_csv('../data/neg_bacilli_turin_2024.txt', delimiter='\t')
print(df)
for name in df.columns:
    print(name)
````


## Missing Data

```python
print("Total missing # :", df.isna().sum().sum())
print("Total missing % :", 100* df.isna().sum().sum() / (df.shape[0]*df.shape[1]) )
df.fillna(0, inplace=True)
print("Total missing # :", df.isna().sum().sum())
print("Total missing % :", 100* df.isna().sum().sum() / (df.shape[0]*df.shape[1]) )
````

## Summary Stats

```python
def calculate_se_proportion(p, n):
    return np.sqrt((p * (1 - p)) / n).round(2)
````


```python
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
````



```python
summary_table = create_summary_table(df)
summary_table.to_csv('../data/summary_table.csv')
summary_table
````



# Summary Stats + Fisher's Exat Test | Wilcoxon/Mann-Whitney Test



`
Summary statistics including mean, standard deviation (SD), median, interquartile range (IQR), proportion, and standard error (SE) were calculated to provide insights into the distribution and variability of the data.
For binary variables, proportions of patients with specific events were calculated, along with their corresponding standard errors. We employed Fisher's exact test to assess the association between binary variables and outcomes. 
For continuous variables, we used the Mann-Whitney U test to compare distributions between groups based on outcome status.
`

## Positive >2 cultures  || +50%

```python
df_cult_pos = df[df['Positive_2culture_50perc'] == 1] # 43 out of 63
summary_table = create_summary_table(df_cult_pos)
summary_table.to_csv('../data/summary_table_cult_pos.csv')
````



```python
df_cult_neg = df[df['Positive_2culture_50perc'] == 0] # 20 out of 63
summary_table = create_summary_table(df_cult_neg)
summary_table.to_csv('../data/summary_table_cult_neg.csv')
````




```python
exclude_cols = ['Positive_2culture_50perc', 'pat_id']
binary_results = []
numeric_results = []
````




```python
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
````



```python
print("Numeric Variables (Wilcoxon/Mann-Whitney Test): ")
for result in numeric_results:
    print(f"{result[0]} , p-value: {result[1]}" )
````



```python
print("Binary Variables (Fisher's Exact Test): ")
for result in binary_results:
    print(f"{result[0]} - Odds Ratio: {result[1]}, p-value: {result[2]}" )
````


# Summary Stats + Fisher's Exat Test | Wilcoxon/Mann-Whitney Test


## In-Hospital Mortality


```python
df_hosp_mort_pos = df[df['Mortality_hospital'] == 1] # 10 out of 63
summary_table = create_summary_table(df_hosp_mort_pos)
summary_table.to_csv('../data/summary_table_hosp_mort_pos.csv')
````


```python
df_hosp_mort_neg = df[df['Mortality_hospital'] == 0] # 53 out of 63
summary_table = create_summary_table(df_hosp_mort_neg)
summary_table.to_csv('../data/summary_table_hosp_mort_neg.csv')
````

```python
exclude_cols = ['Mortality_hospital', 'pat_id', 'Amoxicillin_clavulanate_including', 'Metronidazole_including', 'Mortality_30_days', 'Mortality_14_days']
binary_results = []
numeric_results = []
````

`
Features containing only zeros have been removed as they provided no addiitonal information.
`
`
Mortality at 14 and 30 days post hospitalization has been ignored.
`



```python
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
````


```python
print("Numeric Variables (Wilcoxon/Mann-Whitney Test): ")
for result in numeric_results:
    print(f"{result[0]} , p-value: {result[1]}" )
````



```python
print("Binary Variables (Fisher's Exact Test): ")
for result in binary_results:
    print(f"{result[0]} - Odds Ratio: {result[1]}, p-value: {result[2]}" )
````


# Model: XGBoost

```python
df = pd.read_csv('../data/neg_bacilli_turin_2024.txt', delimiter='\t')
df.fillna(0, inplace=True)
````


```python
columns_to_drop = ['pat_id', 'Amoxicillin_clavulanate_including', 'Metronidazole_including', 'Mortality_30_days', 'Mortality_14_days']
df = df.drop(columns=columns_to_drop)
````


```python
y = df['Mortality_hospital']
X = df.drop(columns=['Mortality_hospital'])
````


```python
X100 = shap.utils.sample(X, 100)
df.shape
````


```python
model = xgboost.XGBClassifier(nestimators=100, max_depth=2).fit(X, y)
explainer = shap.Explainer(model, X)
shap_values = explainer(X)
````


`
We employed XGBoost as an explanatory model to rank and identify relative feature importance, using  SHAP values as a proxy for feature importance measured at the local data point level. SHAP values are reported as log odds.
`

```python
shap_values.display_data = X.values
````


```python
shap.plots.bar(shap_values)
````

```python
shap.plots.bar(shap_values.abs.max(0))
````

```python
shap.plots.beeswarm(shap_values)
````

```python
shap.plots.heatmap(shap_values)
````

```python
shap.initjs()
shap.force_plot(shap_values[10])
````


` 
(a) Beeswarm Plot of SHAP Values: Beeswarm plot illustrateing the distribution of SHAP (SHapley Additive exPlanations) values for each feature in the dataset. Each dot represents a feature value for a specific patient, with the position along the x-axis indicating the magnitude of the SHAP value. The color of each dot indicates the corresponding feature value, providing insight into the relationship between feature values and their impact on model predictions. Features with wider distributions and greater dispersion of SHAP values suggest higher variability and importance in the model's decision-making process. 
Redder dots on the right side of the plot indicate that being positive for that feature or experiencing higher values on it increases the predicted probability of death. Redder dots on the left side indicate the opposite (lower chances of dying).
(b) Heatmap of SHAP Values:  Heatmap illustrating the SHAP (SHapley Additive exPlanations) values for the top features impacting the model's mortality prediction scores.  Each row corresponds to a feature in the dataset, and each column represents a patient sample. 
Patient samples are ordered using hierarchical clustering by their explanation similarity, resulting in samples with closer model outputs for the same reason getting grouped together.
The color intensity indicates the magnitude and direction of the feature's impact on the model output: red indicates positive impact (increasing the predicted likelihood of death), while blue indicates negative impact (decreasing the predicted likelihood of death).  
A feature's importance can be inferred from the range and variability of its SHAP values across samples.
Features with higher absolute SHAP values exert a greater influence on the model predictions. Only the top 9 features are herein individually depicted.
The output of the model is shown above the heatmap matrix as a line plot centered around the explaination’s base valu and the global importance of each model input shown as a bar plot on the right hand side of the plot.
(c) Force plots with individual patient examples, breaking down the contribution of each feature to the prediction of a given patient (3 random patients with a high/low predicted score shown). Scores are on a log odds scale. Probabilities can be easily infered as probability = exp(log-odds)/(1+exp(log-odds)).
`


`
There are some trends, whereby patients who did die tend to have received combination therapy as oposed to monotherapy, were more comorbid (most notably in regard to diabetes). However, by far, the only features strongly predictive of in-hospital mortality seem to be time from admission to BSI and total length of stay. But this is, of course, not very helpful. The longer a patient stays on the hospital, the more likely he/she is to die at the hospital.
`


# Principal Component Analysis

```python
df = pd.read_csv('../data/neg_bacilli_turin_2024.txt', delimiter='\t')
df.fillna(0, inplace=True)
````


```python
columns_to_drop = ['pat_id', 'Amoxicillin_clavulanate_including', 'Metronidazole_including', 'Mortality_14_days', 'Mortality_30_days', 'LOS_days']
df = df.drop(columns=columns_to_drop)
````


```python
numerical_features = df.select_dtypes(include=['float64', 'int64'])
numerical_features = numerical_features.drop(columns=['Mortality_hospital'])

binary_numerical_features = [col for col in numerical_features.columns if numerical_features[col].nunique() == 2]

non_binary_numerical_features = numerical_features.drop(columns=binary_numerical_features)

scaler = StandardScaler()

scaled_numerical_features = scaler.fit_transform(numerical_features)
````


```python
df[numerical_features.columns] = scaled_numerical_features
df.shape
````


```python
numerical_features = df.drop(columns=['Mortality_hospital']).select_dtypes(include=['float64', 'int64'])
numerical_features.shape
````

```python
pca = PCA(n_components=0.75) 
pca.fit(numerical_features)
````

```python
principal_components = pca.transform(numerical_features)
principal_components_df = pd.DataFrame(data=principal_components, 
                                       columns=[f'PC{i}' for i in range(1, pca.n_components_ + 1)])
````

```python
loadings = pca.components_
loadings_df = pd.DataFrame(loadings, columns=numerical_features.columns)
````

```python
for i, component in enumerate(loadings_df.iterrows(), start=1):
    print(f"Principal Component {i} Loadings:")
    print(component)
    print("\n")
````


```python
loadings_df.index = loadings_df.index + 1

plt.figure(figsize=(9, 16))
sns.heatmap(loadings_df.T, cmap='RdGy', annot=False, fmt=".1f", cbar=True)
plt.title('Principal Component Loadings')
plt.xlabel(None)  
plt.ylabel(None)  
plt.show()
````


```python
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
````


```python
principal_components_df

y = df['Mortality_hospital']

correlation_with_y = principal_components_df.apply(lambda col: col.corr(y))

````


```python
correlation_with_y_sorted = correlation_with_y.abs().sort_values(ascending=False)

print("Correlation with Output Variable (Absolute Values):")
print(correlation_with_y_sorted)
````



```python
pca_with_y = pd.concat([principal_components_df, y], axis=1)

mean_loadings_by_y = pca_with_y.groupby('Mortality_hospital').mean()

print("Mean Loadings by Output Variable Group:")
print(mean_loadings_by_y)
````



```python
mean_loadings_by_y.index = mean_loadings_by_y.index.map({0: 'No', 1: 'Yes'})


plt.figure(figsize=(25, 1))
sns.heatmap(mean_loadings_by_y, cmap='RdGy', annot=True, fmt=".1f", cbar=True)
plt.title('Mean Loadings by In-hospital Mortality Status')
plt.xlabel('Principal Components')
plt.ylabel('In-hospital Mortality')
plt.show()
````




```python
top_features_by_pc = {}

for pc in loadings_df.index:
    top_features_by_pc[pc] = loadings_df.loc[pc].abs().nlargest(5).index.tolist()

for pc, features in top_features_by_pc.items():
    print(f"Principal Component {pc}: {features}")

````


`
Variance is spread too thin throughout. Upon varimax roation the loadings are very weak and even the more modest ones (0.3 - 0.4 max) have nothing to do with the mortality-related features which don't appear until PC +30. Which reinforces that the patients are not significantly different at baseline. It is probably just a bias due to longer  hospitalization times.
`
