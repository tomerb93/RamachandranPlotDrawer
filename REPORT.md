# RamachandranPlotDrawer

## For computer science course Intro to BioInformatics:

### What I did
I initially read the PDB file, storing the protein in a variable and 
making sure to skip any atomic acid different from N,CA,C in the chains.
After reading the file, in a nested for loop over each chain in the loaded protein, I calculated each 
phi angle by using the previous amino acid's C value and the current amino acid's N,CA and C values.
I calculated each psi angle by using the current AA's N,CA,C values and next AA's N value.
Then, I calculated the vectors of these coordinates, and their cross product vectors n1,n2 (like shown in class).
The acos of the dot product of n1,n2 (normalized) gave me the resulting phi/psi angle in radians.
At this point I was only getting positive numbers between 0 and PI, so in order to determine the sign
needed I calculated the dot product of the (cross product of n1,n2) and r3. If the result was greater than 0, 
the sign would be negated.
Once I calculated the values I added each phi/psi tuple into list which I then scattered on a plot using matplotlib

### Differences
The main differences I saw between different proteins of the same type (A,B,A+B or A/B) is the amount of vertices describing them,
apart from that they were usually pretty similar.
Different types varied in the amount of alphas (A types had more concentration around alpha area) and betas (B types had more concentration around beta area),
but also usually all concentrating around same 2 areas, and loops also appearing in generally similar areas between all types.

### Results
Images of the ramachandran plot drawn can be found under the directory "results" for the proteins found in the "proteins" directory
