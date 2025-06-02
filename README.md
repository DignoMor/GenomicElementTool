# GenomicElementTool
Tool to work with genomic elements

## Genomic Element

Genomic Element is the key data structure this tool is based on. 
A Genomic Element contains 3 key components:

- A reference genome file
- A region file
- Annotations

### Reference Genome file

Reference genome file of a Genomic Element is of fasta format. 

### Region file

A region file defines the coordinates of the "elements". 
Region files are one of the few supported variants of 
bed format. The coordinates follow the same half open 
convension as bed files.

### Annotations

Annotations are numpy arrays stored in npy files. 
The first dimension of the array always equal to 
the number of regions. Depending on the second dimension 
size, an annotation can be of different types.

- track: store signal tracks for each element. 
  The second dimension is larger than 1.
- stat: store statistics for each element, the 
  second dimension is of size 1.
