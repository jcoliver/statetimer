# statetimer

States through time plotting package (incomplete)

Contact [Jeff Oliver](mailto:jcoliver@arizona.edu?subject=statetimer%20package)
if you have any questions.

## TODO:
+ `joint_probs()` documentation says that x and y are ancestral state 
reconstructions, indicating they are the states element (a matrix) of a rayDISC 
object. But the function itself tries to extract the states element from x and 
y. The documentation is probably wrong, because other elements (e.g. 
`tip.states`) are also extracted from x and y by `joint_probs()`
