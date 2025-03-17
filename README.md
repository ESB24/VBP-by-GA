# VBP


## TODO miscelaneous

## TODO Presentation

 - la poste pic: https://www.gettyimages.fr/search/2/image?phrase=mail+sorting+machine

## TODO article

Read:
  - `Mommessin2023`: "Affinity-aware resource provisioning for long-running applications in shared clusters" (l. 106 of .bib)
  - `Bansal2010`: "A New Approximation Method for Set Covering Problems, with Applications to Multidimensional Bin Packing" (l. 179 of .bib)

## TODO algorithm

### Bounds
 - Upper bound :
   - Tests : `biggest mail`, `total volume`, `route number`
   - Heuristic : `best-fit-decreasing`, sort lexico by = (`biggest minus smallest mail`, `|mail|`)
   - Insertion : `smooth assigned`, ?`left alligned`?

### Partitioning
 - create a real 1D-bin-packing solver (model or (meta)heuristic)
 - rely only on the partitionning
 - column generation ?
 - 

### Rebuild session
 - minimize loss obj fct ?
 - dinamic weighting of the objective function! (rely on prevision, backtrack, ...)
 - change model update in back-tracking! (volume cst)
 -  

### Miscelaneous
 - local search ?
 - make the whole paper readable even if we skip/don't understand the formula/model.
