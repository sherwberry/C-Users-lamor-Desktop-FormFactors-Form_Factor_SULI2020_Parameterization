# Deuteron Form Factor Parameterization

## Introducton

This work was done as part of the 2020 DOE SULI program hosted remotely by Jefferson Lab.

The figures were made using the python matplotlib library in Jupyter notebook and the Abbott parameterizations [1].

## New Parameteriztion

Over the range of the data, the orginal Abbot I and II parameterizations are, in general reasonable representations of the data.   As people tried to use these parameterizatons for the EIC, going beyond the original range of validity of the fits problems were discovered.   Most notibliy an non-physical signularity in the charge form factor was noted.   For this summer project, we redid the fits but adding in constains so that the functions were more physical beyond the range of the data. 

## Results

Interestingly, just adding this constraint, makes a new function that, over the range of the experimental data, nicely falls between Abbot I and II.  The Kiwi in the legends refers to our new fit and is shown plotted with Abbott I and II parameterizations

![G_C](https://github.com/sherwberry/C-Users-lamor-Desktop-FormFactors-Form_Factor_SULI2020_Parameterization/blob/master/g_c_parameterization_weighted.png)
![G_Q](https://github.com/sherwberry/C-Users-lamor-Desktop-FormFactors-Form_Factor_SULI2020_Parameterization/blob/master/g_q_parameterization_weighted.png)
![G_M](https://github.com/sherwberry/C-Users-lamor-Desktop-FormFactors-Form_Factor_SULI2020_Parameterization/blob/master/g_m_parameterization_weighted.png)




## References

[1] D. Abbott (Jefferson Lab) et al., Eur.Phys.J.A. 7 (2000) 421-427. <a href="http://doi.org/10.1007/PL00013629">DOI: 10.1007/PL00013629</a>

## Zenodo

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4074280.svg)](https://doi.org/10.5281/zenodo.4074280)
