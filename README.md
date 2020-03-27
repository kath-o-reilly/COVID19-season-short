# COVID19-season-short

## Repo outline

This repo is provided to replicate the short analysis we've been working on. Please see https://cmmid.github.io/topics/covid19/current-patterns-transmission/role-of-climate.html for the comment. Please note that this is merely a comment and is not peer-reviewed.

### Data

All the data used to carry out the analysis is provided in the _Data_ folder. Please see _WHO_SITREP_COVID_24032020_incTemp_v2.csv_ for the details of each country.

### Code

The analysis was carried out in R. The analysis is very simple; we combine data together from several sources;
- The WHO 'sitreps' from the 24th March 2020, where the transmission status of each country is provided.
- ERA5 hourly 2-metre temperature and dewpoint temperature data was extracted from the Copernicus Climate Change Service Data Store. Relative humidity was calculated using the August–Roche–Magnus formula. Seasonal averages for each country are population adjusted using data from worldpop and publicly available shapefiles.

Any errors in the code should be notified to kathleen.oreilly@lshtm.ac.uk