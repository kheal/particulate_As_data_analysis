# Read me

This github repository is the data analysis portion associated with the manuscript by Heal *et. al* "Diverse arsenic-containing lipids in the surface ocean", published in Limnology and Oceanography, DOI: 10.1002/lol2.10216.


#### **Control.Rmd**
The **Control.Rmd** file is the data analysis notebook, which creates arsenoliplid database, cleans ICPMS and QE data for further analysis, and calculates bulk and individual arsenic concentrations.

* Sources scripts from 'SourceCode' folder
* Sources files from 'MetaData' and 'RawOutput' folders
* Writes products into the 'Intermediates' folder

#### **Exploration.Rmd**
The **Exploration.Rmd** file is the data exploration notebook, which uses the arsenolipid databsae and ESI data to search for possibile matches in MS1 and MS2 space for arsenolipids.  After these explorations, I manually inspected the output and curated a list of compounds saved in 'MetaData/Targeted_visually_IDd_2.csv' 

* Sources scripts from 'SourceCode' folder
* Sources files from 'MetaData' and 'RawOutput' folders
* Writes products into the 'Intermediates' folder 


#### **Figures.Rmd**
The **Figures.Rmd** file is the figure creation notebook. Main text and supplemental figures are both made and not separated in the output.

* Sources scripts from 'Figures/Code' folder
* Sources files from 'Intermediates'
* Writes products figures into the 'Figures/Manuscript_figures' 

#### **Tables.Rmd**
The **Tables.Rmd** is the table creation notebook. Tables are made into csv tables and supplemental tables are combined into a .xcel file.  

* Sources scripts from 'Tables/Code' folder
* Sources files from 'Intermediates' folder 
* Writes products figures into the 'Manuscript_tables' and 'Manuscript_tables/SuppTables' folders
