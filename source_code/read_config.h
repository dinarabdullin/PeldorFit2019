#ifndef READ_CONFIG_H
#define READ_CONFIG_H

#include "definitions.h"
#include "spin.h"
#include <vector>

void read_config(const char* path_config, 
	             std::vector<experiment>& experiments,
	             std::vector<Spin>& spin_system,
                 optimization_parameters& opt_param,
				 genetic_parameters& genetic_param,
				 output_parameters& output_param,
				 errorplot_parameters& errorplot_param);
#endif