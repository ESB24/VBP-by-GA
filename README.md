# VBP


## TODO

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
