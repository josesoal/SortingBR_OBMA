#ifndef OPERATORS_H_
#define OPERATORS_H_

#include "structs_ga.h"

void generate_initial_population(population *pob, permutation *perm);
void restart_population_LS(population *pob, permutation *perm, int *num_eval_fit_function);
void restart_population_OBL_LS(population *pob, permutation *perm, int *num_eval_fit_function, int *flag_obl);
int opposite_permutation(population *pob, int i, int start, int end, int *num_eval_fit_function);
void calculate_fitness_population(population *pob, int generation_number, int *num_eval_fit_function);
void selection(population *pob, int generation_number, int *best_solution);
void crossover(population *pob, permutation *perm);
void mutation(population *pob, permutation *perm);
void replacement(population *pob);
void LS_for_population(population *pob, int generation_number,int *num_eval_fit_function);

void OBL_LS_for_population(population *pob,int generation_number, int *num_eval_fit_function,int *flag_obl);
int local_search(population *pob, int i, int *num_eval_fit_function);
int obma_reached_degenerate_state(population *pob);

void read_permutation(permutation *perm);
void read_permutation_from_file(char* nombre_archivo, permutation *perm);
void show_population(population *pop, permutation *perm, int generation_number);

#endif /* OPERADORES_H_ */
