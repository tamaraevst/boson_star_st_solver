In the source folder, you may have noticed two versions for the code that computes the 1D solution for the massive gravitational scalar field. Each one of them is used for:
* ***massive_no_match.c*** - shoot for bosonic frequency and then compute the solution straight away with no matching (at the radius where the solutions grow too quickly, the respective fields are set to zero so that end profiles can be recorded).
* ***massive.c*** - the code to compute teh solution using a 'smoothness' criteria that uses Newton-Raphson method to obtain a physical solution.
