# Network-based gene prioritization analysis with NetGPA
it has been demonstrated that genes who function in the same pathway tend to distribute in a coherent sub-network. NetGPA implements a network-based gene prioritization method. Briefly, a sub-network is built using a given set of seed genes that are assumed to function in the same pathway or similar pathways. To predict whether a query gene is functionally related to the seed genes, we project this gene to a global network, and test whether connection of this gene to the subnetwork is random or statistically significant.

## Installation
library("devtools")
install_github("JiantaoShi/NetGPA")

## Vignettes
Please read this vignettes [vignettes](https://jiantaoshi.github.io/Package/NetGPA_vignettes.html) for introduction and application.

## Contacts
If you have any questions or feedback, please contact us at:
Email: jshi@jimmy.harvard.edu

