/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 


	setup_stalk(); 
	setup_tip(); 

	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 
	
	// create some of each type of cell 
	
	Cell* pC;
	
	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 
	
	// load cells from your CSV file (if enabled)
	load_cells_from_pugixml(); 	

	
	// make any cells near boundary "rigid"
	for( int n=0; n < (*all_cells).size() ; n++ )
	{
		Cell* pC = (*all_cells)[n]; 
		if( fabs(pC->position[0]) > 450 || fabs(pC->position[1]) > 450 )
		{
			// set_single_behavior( pC, "is movable" , 0); 
			pC->is_movable = false; 
		} 
	}

	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ 
	std::vector<std::string> out = paint_by_number_cell_coloring(pCell); 
	return out; 

//	if( pCell->type_name == "stalk")
	{
		int R = (int)( get_single_signal( pCell , "custom:tip_signal" ) * 255.0 );
		if( R > 255 )
		{ R = 255; }

		std::string color = "rgb(" + std::to_string(R) + ",128,128)"; 
		out[0] = color;  
	}


	return out; 
}

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 

void setup_stalk( void )
{
	Cell_Definition* pCD = find_cell_definition( "stalk"); 

	pCD->phenotype.mechanics.maximum_number_of_attachments = 2; 
	pCD->phenotype.mechanics.attachment_rate = 10; 
	pCD->phenotype.mechanics.attachment_elastic_constant = 0.01; 

	pCD->functions.update_phenotype = stalk_phenotype; 
	pCD->functions.contact_function = flux_contact; 
	pCD->functions.custom_cell_rule = custom_helper; 

	return; 
}

void stalk_phenotype( Cell* pCell , Phenotype& phenotype , double dt )
{
	tip_signaling(pCell,phenotype,dt);

	double tip = get_single_signal( pCell, "custom:tip_signal"); 
	double c = get_single_signal( pCell, "chemoattractant"); 
	bool movable = (bool) get_single_behavior( pCell, "is movable"); 

	if( movable == false )
	{ return; }

	// if c is high and tip is low, convert to tip 
	double r_transform = 0; 
	double c_threshold = 0.003 ; 0.001; // 0.01 
	if( c > c_threshold && tip < 0.01 )
	{ r_transform = 0.005; } 
	set_single_behavior( pCell, "transform to tip" , r_transform) ; 

	int n_adhesions = pCell->state.attached_cells.size(); 
	set_single_behavior( pCell, "custom:n_adhesions" , n_adhesions ); 

	return; 
}

void flux_contact( Cell* pMe, Phenotype& phenoMe, Cell* pOther, Phenotype& phenoOther, double dt )
{
	// flux constant across the connected cells 
	double k = 0.5*( get_single_signal( pMe , "custom:tip_flux" ) + get_single_signal( pOther, "custom:tip_flux") ); 

	// if no flux, exit
	if( k < 1e-10 )
	{ return; }

	double signal_me = get_single_signal( pMe, "custom:tip_signal" ); 
	double signal_other = get_single_signal( pOther, "custom:tip_signal" ); 

	// one-directional flux: from other into me
	if( signal_other > signal_me )
	{
		double d_signal = dt * k * (signal_other - signal_me); 
		signal_me += d_signal;
		signal_other -= d_signal; 

		set_single_behavior( pMe , "custom:tip_signal" , signal_me ); 
		#pragma omp critical 
		{ set_single_behavior( pOther , "custom:tip_signal" , signal_other); }
	}

	return; 
}

	/* tip signaling */
void tip_signaling( Cell* pCell, Phenotype& phenotype, double dt )
{
	double r_create = get_single_signal( pCell , "custom:tip_creation");
	double r_decay = get_single_signal( pCell, "custom:tip_decay"); 

	double tip = get_single_signal( pCell , "custom:tip_signal"); 
	// std::cout << tip << std::endl; 
	double d_tip = dt * r_create * ( 1.0 - tip ) - dt * r_decay * tip; 
	tip += d_tip; 
	// std::cout << "\t" << d_tip << " " << tip << std::endl; 

	set_single_behavior( pCell, "custom:tip_signal" , tip ); 

	return; 
}

void custom_helper( Cell* pCell, Phenotype& phenotype , double dt )
{
	// copy current attachmetns 
	pCell->state.attached_cells = pCell->state.spring_attachments; 
}

void setup_tip( void )
{
	Cell_Definition* pCD = find_cell_definition( "tip"); 

	pCD->phenotype.mechanics.maximum_number_of_attachments = 2; 
	pCD->phenotype.mechanics.attachment_rate = 0; 
	// pCD->phenotype.mechanics.detachment_rate = 10; 
	pCD->phenotype.mechanics.attachment_elastic_constant = 0.01; 

	pCD->functions.update_phenotype = tip_phenotype; 
	pCD->functions.contact_function = flux_contact; 
	pCD->functions.custom_cell_rule = tip_custom; // custom_helper; 

	return; 
}

void tip_phenotype( Cell* pCell , Phenotype& phenotype , double dt )
{
	/* tip signaling */
	tip_signaling(pCell,phenotype,dt); 

	int n_adhesions = pCell->state.attached_cells.size(); 

	set_single_behavior( pCell, "custom:n_adhesions" , n_adhesions ); 

	return; 
}

void tip_custom( Cell* pCell , Phenotype& phenotype , double dt )
{
	custom_helper(pCell,phenotype,dt);

	// check for intersection with source 
	for( int n=0; n < pCell->state.neighbors.size() ; n++ )
	{
		Cell* pC = pCell->state.neighbors[n]; 
		if( pC->type_name == "source")
		{
			// if found, attach as spring
			// and rapidly transform back to stalk 
			attach_cells_as_spring( pCell, pC ); 
			set_single_behavior( pCell , "transform to stalk" , 10); 

			// also, turn off the source
			#pragma omp critical 
			{ set_single_behavior( pC, "chemoattractant secretion", 0); }
		}
	}

	return;
}

Cell* clone_cell( Cell* pSource, std::vector<double>& direction )
{
	static Cell_Definition* pCD_tip = find_cell_definition( "tip"); 

	Cell* pNew = create_cell( *pCD_tip ); 

	std::vector<double> x = pSource->position; 
	double distance = pSource->phenotype.geometry.radius + pNew->phenotype.geometry.radius; 
	direction *= distance; 
	x += direction; 
	pNew->assign_position( x );

	return pNew; 
}

void custom_extension( double dt )
{
	static double t = 0; 
	static double interval = 60; 
	static double t_next = 0.0;

	static int n_chemo = microenvironment.find_density_index( "chemoattractant");
	static int n_inhib = microenvironment.find_density_index( "inhibitory"); 
	// static int n_inhib = 0;

	t += dt; 
	if( t < t_next )
	{ return; }

	double extension_rate = 0.01; 

	int n_last = (*all_cells).size(); 

	for( int n=0 ; n < n_last ; n++ )
	{
		Cell* pC = (*all_cells)[n];
		if( pC->type_name == "tip" )
		{
			// extension direction 
			std::vector<double> x = pC->position; 
			std::vector<double> direction = microenvironment.nearest_gradient_vector(x)[n_chemo];
			std::vector<double> direction_inhib = microenvironment.nearest_gradient_vector(x)[n_inhib]; 
			normalize( &direction ); 
			normalize( &direction_inhib);

			// std::vector<double> rand = UniformOnUnitCircle(); 
			double bias = 0.70;
			double neg = -1.0;
			direction = bias * direction + (1-bias) * (neg * direction_inhib); 
			normalize( &direction ); 
			
			// extend IF chemokine is high enough 
			double c_threshold = 0.001; 
			double in_threshold = 0.15;
			if( get_single_signal(pC,"chemoattractant") > c_threshold )
			{ 
				if( get_single_signal(pC,"inhibitory") < in_threshold )
				{
				// clone (and position)
				Cell* pNew = clone_cell( pC , direction ); 

				// spring-link 
				attach_cells_as_spring( pC , pNew ); 

				// force fast conversion of the prior tip 
				set_single_behavior( pC, "transform to stalk" , 10); 
				}
			}
		}

	}	

	t_next += interval; 
	return; 
}
