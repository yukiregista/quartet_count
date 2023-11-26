# quartet_count

Given a set of trees, it returns the most occuring quartet topologies. Can be used to create an input to qstar algorithm.

I only tested this code with mac OS and linux.

## Inputs/Outputs
Input: (pure) newick file containing multiple trees. It should only contain newick trees (possibly with branch lengths, although we don't account for branch lengths now), without any annotaions or rooting information (like [\&U]) that is usually added for nexus files. See `examples/example.tre` for an example.

Output: Most occuring quartet topologies. See `examples/example_quartfile` for the structure of the output. In case of a tie, no quartet topology is returned for that quartet.


## Usage
### Download/Compiles
```
git clone https://github.com/yukiregista/quartet_count.git
cd quartet_count
make
```

### Run the program
You need to pass --input (input file path) and --output (output file path).
```
./quartet_count --input examples/example.tre --output examples/example_quartfile
```
