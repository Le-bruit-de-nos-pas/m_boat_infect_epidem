import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from pypalettes import load_cmap
from matplotlib.font_manager import FontProperties
from drawarrow import fig_arrow
from pyfonts import load_font

from highlight_text import fig_text, ax_text


url = "https://raw.githubusercontent.com/holtzy/the-python-graph-gallery/master/static/data/europe.geojson"
world = gpd.read_file(url)
world.head()


#url = "https://raw.githubusercontent.com/holtzy/the-python-graph-gallery/master/static/data/co2PerCapita.csv"
#df = pd.read_csv(url)


df = pd.read_csv('Countries_Isolates.csv')
df.head()

df.Country.unique()

data = world.merge(df, how='left', left_on='name', right_on='Country')

data.name.unique()

data = data[data['continent'] == 'Europe']

data = data[data['name'].isin(['Slovakia', 'Romania', 'Bulgaria', 'Czechia', 'Denmark',
       'Hungary', 'Greece', 'Belgium', 'Ireland', 'Austria', 'Switzerland',
       'Italy', 'Netherlands', 'Poland', 'Slovenia', 'Portugal', 'Spain',
       'Norway', 'Germany', 'France', 'Sweden', 'Croatia'])]

data = data[['name', 'Peds', 'geometry']]
data = data.dropna()
data

adjustments = {
    'France': (10, 3),
    'Italy': (-2.4, 2.5),
    'Ireland': (0, -1),
    'Germany': (-0.2, 0),
    'Poland': (0, 0.2),
    'Sweden': (-1.2, -2.8),
    'Norway': (-4, -5.5),
}




cmap = load_cmap('BrwnYl', cmap_type='continuous')
background_color = 'white'
text_color = 'black'


font = load_font(
   'https://github.com/dharmatype/Bebas-Neue/blob/master/fonts/BebasNeue(2018)ByDhamraType/ttf/BebasNeue-Regular.ttf?raw=true'
)
other_font = load_font(
   'https://github.com/bBoxType/FiraSans/blob/master/Fira_Sans_4_3/Fonts/Fira_Sans_TTF_4301/Normal/Roman/FiraSans-Light.ttf?raw=true'
)
other_bold_font = load_font(
   'https://github.com/bBoxType/FiraSans/blob/master/Fira_Sans_4_3/Fonts/Fira_Sans_TTF_4301/Normal/Roman/FiraSans-Medium.ttf?raw=true'
)



fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
fig.set_facecolor(background_color)
ax.set_facecolor(background_color)

data.plot(ax=ax, column='Peds', cmap=cmap, edgecolor='black', linewidth=0.5, alpha=0.7)

ax.set_xlim(-11, 41)
ax.set_ylim(32, 73)
ax.set_axis_off()



data_projected = data.to_crs(epsg=3035)
data_projected['centroid'] = data_projected.geometry.centroid
data['centroid'] = data_projected['centroid'].to_crs(data.crs)

countries_to_annotate = [
       'Slovakia', 'Romania', 'Bulgaria', 'Czechia', 'Denmark',
       'Hungary', 'Greece', 'Belgium', 'Ireland', 'Austria', 'Switzerland',
       'Italy', 'Netherlands', 'Poland', 'Slovenia', 'Portugal', 'Spain',
       'Norway', 'Germany', 'France', 'Sweden', 'Croatia'
]

for country in countries_to_annotate:
   centroid = data.loc[data['name'] == country, 'centroid'].values[0]
   x, y = centroid.coords[0]
   try:
      x += adjustments[country][0]
      y += adjustments[country][1]
   except KeyError:
      pass
   rate = data.loc[data['name'] == country, 'Peds'].values[0]
   if country=='United Kingdom':
      country = 'UK'
   if rate > 70000:
      color_text = 'white'
   else:
      color_text = text_color # 'black'
   ax_text(
      x=x, y=y, s=f"<{country.upper()}>: {rate:.0f}", fontsize=6, font=other_bold_font, color=color_text,
      ha='center', va='center', ax=ax, highlight_textprops=[{'font': other_bold_font}]
   )

plt.savefig('isolates.png', dpi=300)
plt.show()