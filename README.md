#  Introduction

This repository presents a Bayesian model designed by Richard McElreath allowing to evaluate the size of signal repertoire (vocal/gestural/other) at the individual level.
The model accounts for 
  - The duration individuals have been observed
  - The number of recordings collected
  - The frequency at which each individual vocalizes
  - The commonness of each signal type

#  Applicability
This model can thus be used for samples with unbalanced datasets on each individual and can be used to compare individuals across age, sex, populations and even species. 
It can be apply to single signals but also, like here, to signal combinations (call sequences or gesture sequences). 
A detail description of the model can be found at:  https://doi.org/10.1038/s42003-025-07922-2 

#  Citation
Sigmundson, R., Girard-Buttoz, C., Le Floch, A., Souha Azaiez, T., McElreath, R., Zuberbühler, K., Wittig, R. M., & Crockford, C. (2025). Vocal sequence diversity and length remain stable across ontogeny in a catarrhine monkey (Cercocebus atys). Communications Biology, 8:465. https://doi.org/10.1038/s42003-025-07922-2

#  Content
The script: "Mangabey ontogeny Model Rethinking version.R" provides the code for the model.
The "input" folder comprises all the datasets used in the paper that apply the same model structure on them.  

#  Resources needed
Library needed for running the model: Rethinking Package: https://github.com/rmcelreath/stat_rethinking_2024
