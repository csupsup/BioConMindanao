# Caught in statistical noise: Pitfalls of a unidimensional approach to understanding biodiversity-conflict in the southern Philippines
Authors: Kier Mitchel E. Pitogo, Camila G. Meneses, Andrie Bon A. Flores, Aljohn Jay L. Saavedra, Ace Kevin S. Amarga, Marjorie Delos Angeles, Cristian C. Lucañas, Syrus Cesar Decena, Russell Evan L. Venturina, Jay S. Fidelino, Dexcem Pantinople, Mark W. Herr, Justin M. Bernstein, Kin Onn Chan, Marites B. Sanguila, Neil Aldrin Mallari, Rafe M. Brown, and Christian E. Supsup

This repository documents our reanalysis of the dataset used to examine biodiversity and conflict relationships in the Southern Philippines presented in this paper https://doi.org/10.1038/s44185-024-00044-8

Our critique can be accessed here: https://www.nature.com/articles/s44185-025-00088-4

Their response is available here: https://www.nature.com/articles/s44185-025-00089-3

Some key points raised in their response:

1. They argue that our findings on the relationships between observed species richness and distance differ due to varying analytical choices and decisions. However, they overlook that our critique originated from their failure to provide sufficient details on how they processed the data and conducted the analyses. For instance, they challenged our decision to log-transform independent variables, such as conflict frequency values. We opted for this transformation on the assumption that their values might have been log-transformed as well, given the presence of negative values in Figure 4a of their original paper, which can arise from data transformation. We requested clarification about these values in our correspondence, but they did not respond. Once again, the discrepancies observed could have been avoided had they offered sufficient information about their methods, allowing us to replicate their work accurately.

2. They also raised concerns about our findings regarding the decreasing trend of observed species richness with increasing distance. We acknowledge that we did not provide an in-depth explanation on this matter in our critique. However, to us, this trend primarily stems from the inherent bias in the methodology they employed to link conflict data with biodiversity points using the QGIS plugin "Join attributes by nearest." This approach utilizes an algorithm that connects conflict data to the nearest biodiversity points, resulting in a greater amount of biodiversity data at closer distances compared to those farther away. Consequently, the decline in observed species richness with distance is merely a byproduct of this method. Any interpretations or conclusions relating sociopolitical conflict to biodiversity should be drawn with considerable caution, particularly without a thorough understanding of the method. 

Despite our call to provide adequate information regarding their data treatment and analysis, they proceeded to present another analysis, this time at the municipal level, without supplying any supporting materials, such as cleaned datasets or documentation of their analytical procedures. Therefore, the analyses presented here were conducted to the best of our ability, given the limited information provided in their methods.

Supplementary Files:

FS1: https://csupsup.github.io/BioConMindanao/supplementary_files/Supplementary_FS1.xlsx

FS2: https://csupsup.github.io/BioConMindanao/r_scripts/Supplementary_FS2.html

FS3: https://csupsup.github.io/BioConMindanao/r_scripts/Supplementary_FS3.html
