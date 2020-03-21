## Classes: 

#### Grid_class(caseFile,customCaseName); 
	
```  
	Constructor of an object needs:
	caseFile - MATPOWER case; 
	customCaseName - name your case
```

1. **dcopf();** runs MATPOWER optimal power flow on the loaded case. Saves results to the structure

	```
	somecode
	```

2. **[number_of_violations, absolute_margin, relative_margin, top10] = N_0_analysis();** runs N-0 analysis

	```  
	number_of_violations - is number of N-0 violation  
	absolute_margin      - array of (Limit-Flow) for each line  
	relative_margin      - array of (Limit-Flow)/Limit for each line  
	top10                - top10 the most closest lines
	```

3. **N_1_analysis();**
runs N-1 analysis and saves results to `'/results/N_1_analysis_dangerous_lines_{{custom name}}.mat'`

	```  
	among results are:
	L       - matrix of LODFs
	margins - array of maximum relative margins for all lines 
	lines   - array of numbers, lines(i) = k, tells that i line was violated k times 
	```

4. **N_2_analysis(approach);** 
runs N-2 analysis and saves results to `'/results/'` folder. `Approach = 'fast' or 'bruteforce'`

	```
	among results are:
	cont_fast_algoritm   - massive of dangerous contingencies found by fast algoritm
	cont_brute_force_algorithm - massive of dangerous pairs found by brute force approach'
	```

#### Reduction(); 
class contains various useful reductions of the grid 

1. **new_rnc = remove_parallel(rnc);**
removes parallel lines from the grid `rnc` is a MATPOWER grid case

	```
	new_rnc - case without parallel lines
	```

2. **[area] = connectivity_analysis(E);** 
runs connectivity analysis 

	```
	area - array. area(bus_number) - shows area for 'bus_number' bus. 
	If all areas are the same - grid is connected
	```

3. **[rnc,map] = remap_grid(rnc);**  
remaps odd numbered buses in grid

	``` 
	rnc - remapped case
	map - map for implemented remapping
	```

#### Sz(); 
class for fast and easy calculation of rows, columns and emptiness of any matrix

```
	rows    = Sz.r(matrix); // calculates rows of the matrix
	columns = Sz.c(matrix); // calculates columns of the matrix
	boolean = Sz.z(matrix); // gives true if matrix is empty and falst otherwise  
```