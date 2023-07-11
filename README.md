# plume_modification

These scripts relate to a plume and mixing model for glacial fjords, created as part of a the paper:

Cowton T, D Slater & M Inall (2023). Subglacial-discharge plumes drive widespread subsurface warming in northwest Greenland’s fjords. Geophysical Research Letters (accepted July 2023).

Together, the scripts run a numerical plume model (forced by subglacial discharge and ambient water profiles), then use the plume and ambient water properties to calculate the fraction of plume modified water present in the fjord. For further information, please see the above publication.

The script workflow_example.m runs through the functions required to calculate subglacial-discharge plume
fluxes and water fractions for a given fjord and year based CTD and
runoff inputs

It requires the scripts gsw_sigma0.m and gsw_rho.m available as part of
the Gibbs-SeaWater (GSW) Oceanographic Toolbox
(https://www.teos-10.org/software.htm) (McDougall, T.J. and P.M. Barker,
2011)

It includes example CTD data from Kangerlussuup Sermia in 2019, obtained
from NASA's Ocean's Melting Greenland project archive (OMG 2020)
(https://podaac.jpl.nasa.gov/omg), and timeseries of subglacial discharge
based on daily ice sheet surface runoff outputs from the regional climate
model RACMO2.3, statistically downscaled to 1 km resolution (Noël et al.,
2016)

References:

McDougall, T.J. and P.M. Barker,
2011: Getting started with TEOS-10 and the Gibbs Seawater (GSW)
Oceanographic Toolbox, 28pp., SCOR/IAPSO WG127, ISBN 978-0-646-55621-5.

Noël, B., van de Berg, W. J., Machguth, H., Lhermitte, S., Howat, I.,
Fettweis, X., & Van Den Broeke, M. R. (2016). A daily, 1 km resolution
data set of downscaled Greenland ice sheet surface mass balance
(1958–2015). The Cryosphere, 10(5), 2361-2377.

OMG. (2020). OMG CTD Conductivity Temperature Depth. Ver. 1. PO.DAAC, CA,
USA. https://doi.org/https://doi.org/10.5067/OMGEV-CTDS1
