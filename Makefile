all: quartet_count pick_quartets
	echo "Compiled successfully"

quartet_count: quartet_count.cpp
	g++ --std=c++14 -Ofast -o quartet_count quartet_count.cpp 

pick_quartets: pick_quartets.cpp
	g++ --std=c++14 -Ofast -o pick_quartets pick_quartets.cpp 

clean: 
	rm quartet_count pick_quartets