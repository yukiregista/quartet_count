# quartet_count

Given a set of trees, 

- `quartet_counts` counts number of each quartet topologies.
- `pick_quartets` picks the set of quartet topologies based on counts and picking rule.

This program can be used to create an input to qstar algorithm.

I have only tested this code with mac OS and linux.

## Inputs/Outputs
### quartet_count
Input: (pure) newick file containing multiple trees. It should only contain newick trees (possibly with branch lengths, although we don't account for branch lengths now), without any annotaions or rooting information (like [\&U]) that is usually added for nexus files. See `examples/example.tre` for an example.

Output: Counts of quartet topologies. See `examples/quartcount` for the structure of the output.

### pick_quartets
Input: `quartcount` file produced by `quartet_count`.

Output: The set of quartets picked by user-specified picking rule.


## Usage
### Download/Compile
```
git clone https://github.com/yukiregista/quartet_count.git
cd quartet_count
make
```

### Run the program
For `quartet_count`, you need to pass `--input (input file path)` and `--output (output file path)`.

For `pick_quartets`, you need to pass `--input (input file path)` and `--output (output file path)`.
You can optionally pass `--picking-rule (picking rule)` to specify the rule to pick quartets.
Currently, following picking rules can be specified.
- "most-present": picks quartet topolgies that are present at most. In case of a tie, no topolgoies are included.
- "majority": picks quartet topologies that are present in more than 0.5 proportions of the input trees. This is more conservative than "most-present" rule.
By default, "most-present" rule is applied.

### Example
```
## create quartcount file
./quartet_count --input examples/example.tre --output examples/quartcount
## pick quartet topologies by most-present rule
./quartet_count --input examples/quartcount --output examples/quartfile1 --picking-rule most-present
## This will also pick by most-present rule.
./quartet_count --input examples/quartcount --output examples/quartfile2 --picking-rule most-present
## pick quartet topologies by majority rule
./quartet_count --input examples/quartcount --output examples/quartfile2 --picking-rule majority
```
