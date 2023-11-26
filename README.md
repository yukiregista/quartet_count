# quartet_count

Given a set of trees, it returns the most occuring quartet topologies. Can be used to create an input to qstar algorithm.

I only tested this code with mac OS and linux.

## Inputs/Outputs
Input: (pure) newick file containing multiple trees. It should only contain newick trees (possibly with branch lengths, although we don't account for branch lengths now), without any annotaions or rooting information (like [\&U]) that is usually added for nexus files. See examples/example.tre for an example.

Output: Most occuring quartet topologies. See examples/example.out for the structure of the output. In case of a tie, no quartet topology is returned for that quartet.


## Usage

`make` will compile `quartet_count.cpp` and creates an executable `quartet_count`.

To run the program, you need to pass --input and --output.
```
    quartet_count --input examples/example.tre --output examples/example.out
```
