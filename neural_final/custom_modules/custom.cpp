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
# BSD 3-Clause License (see https://opentarget.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in target and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of target code must retain the above copyright notice,   #
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


	setup_axon(); 
	setup_gcone(); 
   setup_start();
	setup_bridge();
   setup_inhib();

	
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

//	if( pCell->type_name == "axon")
	{
		int R = (int)( get_single_signal( pCell , "custom:gcone_signal" ) * 255.0 );
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

void setup_axon( void )
{
	Cell_Definition* pCD = find_cell_definition( "axon"); 

	pCD->phenotype.mechanics.maximum_number_of_attachments = 2; 
	pCD->phenotype.mechanics.attachment_rate = 10; 
	pCD->phenotype.mechanics.attachment_elastic_constant = 0.01; 

	pCD->functions.update_phenotype = axon_phenotype; 
	pCD->functions.contact_function = flux_contact; 
	pCD->functions.custom_cell_rule = axon_custom; 

	return; 
}

void axon_phenotype( Cell* pCell , Phenotype& phenotype , double dt )
{
	gcone_signaling(pCell,phenotype,dt);
   alive_signaling(pCell,phenotype,dt);

	double gcone = get_single_signal( pCell, "custom:gcone_signal"); 
	double c = get_single_signal( pCell, "chemoattractant"); 
	bool movable = (bool) get_single_behavior( pCell, "is movable"); 

	if( movable == false )
	{ return; }

	// if c is high and gcone is low, convert to gcone 
	double r_transform = 0; 
	double c_threshold = 0.003 ; 0.001; // 0.01 
	if( c > c_threshold && gcone < 0.01 )
	{ r_transform = 0.005; } 
	set_single_behavior( pCell, "transform to gcone" , r_transform) ; 

	int n_adhesions = pCell->state.attached_cells.size(); 
	set_single_behavior( pCell, "custom:n_adhesions" , n_adhesions ); 

	return; 
}

void setup_bridge( void ){

	Cell_Definition* pCD = find_cell_definition( "bridge");
	pCD->functions.update_phenotype = bridge_phenotype;


}

void bridge_phenotype( Cell* pCell , Phenotype& phenotype , double dt ){

   double c = get_single_signal( pCell, "chemoattractant");
   double i = get_single_signal( pCell, "inhibitory");
   double c_thresh = 0.01;
   double inhib_thresh = 0.05;
   if ( c > c_thresh){
      double c_out = c * 0.5;

      set_single_behavior( pCell, "chemoattractant secretion", c_out);
   }

   if (i > inhib_thresh){

      // pCell->phenotype.cycle.advance_cycle(pCell, phenotype, dt);
      double vol = get_single_signal(pCell, "volume");
      double vol_thresh = 2400;
      if (vol > vol_thresh){
         set_single_behavior( pCell, "cycle entry", 0.0023);
         //turn motility speed to 0
         pCell->phenotype.motility.migration_speed = 0;
      }
      else{
         set_single_behavior( pCell, "cycle entry", 0.0);
      }


   }
   else{
         set_single_behavior( pCell, "cycle entry", 0.0);
   }

}

void flux_contact( Cell* pMe, Phenotype& phenoMe, Cell* pOther, Phenotype& phenoOther, double dt )
{
   // flux constant across the connected cells
   double k = 0.5*( get_single_signal( pMe , "custom:gcone_flux" ) + get_single_signal( pOther, "custom:gcone_flux") );

   double k2 = 0.5*( get_single_signal( pMe , "custom:alive_flux" ) + get_single_signal( pOther, "custom:alive_flux") );


   // if no flux, exit
   if( k < 1e-10 || k2 < 1e-10 )
   { return; }

   double gcone_me = get_single_signal( pMe, "custom:gcone_signal" );
   double gcone_other = get_single_signal( pOther, "custom:gcone_signal" );

   double alive_me = get_single_signal( pMe, "custom:alive_signal" );
   double alive_other = get_single_signal( pOther, "custom:alive_signal" );

   // one-directional flux: from other into me
   if( gcone_other > gcone_me )
   {
      double d_signal = dt * k * (gcone_other - gcone_me);
      gcone_me += d_signal;
      gcone_other -= d_signal;

      set_single_behavior( pMe , "custom:gcone_signal" , gcone_me );
      #pragma omp critical
      { set_single_behavior( pOther , "custom:gcone_signal" , gcone_other); }
   }

   if( alive_other > alive_me )
   {
      double d_alive = dt * k2 * (alive_other - alive_me);
      alive_me += d_alive;
      alive_other -= d_alive;

      set_single_behavior( pMe , "custom:alive_signal" , alive_me );
      #pragma omp critical
      { set_single_behavior( pOther , "custom:alive_signal" , alive_other); }
   }

   return;
}

	/* gcone signaling */
void gcone_signaling( Cell* pCell, Phenotype& phenotype, double dt )
{
	double r_create = get_single_signal( pCell , "custom:gcone_creation");
	double r_decay = get_single_signal( pCell, "custom:gcone_decay"); 

	double gcone = get_single_signal( pCell , "custom:gcone_signal"); 
	// std::cout << gcone << std::endl; 
	double d_gcone = dt * r_create * ( 1.0 - gcone ) - dt * r_decay * gcone; 
	gcone += d_gcone; 
	// std::cout << "\t" << d_gcone << " " << gcone << std::endl; 

	set_single_behavior( pCell, "custom:gcone_signal" , gcone ); 

	return; 
}

void custom_helper( Cell* pCell, Phenotype& phenotype , double dt )
{
	// copy current attachmetns 
	pCell->state.attached_cells = pCell->state.spring_attachments; 
}

void setup_gcone( void )
{
	Cell_Definition* pCD = find_cell_definition( "gcone"); 

	pCD->phenotype.mechanics.maximum_number_of_attachments = 2; 
	pCD->phenotype.mechanics.attachment_rate = 0; 
	// pCD->phenotype.mechanics.detachment_rate = 10; 
	pCD->phenotype.mechanics.attachment_elastic_constant = 0.01; 

	pCD->functions.update_phenotype = gcone_phenotype; 
	pCD->functions.contact_function = flux_contact; 
	pCD->functions.custom_cell_rule = gcone_custom; // custom_helper; 

	return; 
}

void gcone_phenotype( Cell* pCell , Phenotype& phenotype , double dt )
{
	/* gcone signaling */
	gcone_signaling(pCell,phenotype,dt); 

	int n_adhesions = pCell->state.attached_cells.size(); 

	set_single_behavior( pCell, "custom:n_adhesions" , n_adhesions ); 

	return; 
}

void gcone_custom( Cell* pCell , Phenotype& phenotype , double dt )
{
	custom_helper(pCell,phenotype,dt);

	// check for intersection with target 
	for( int n=0; n < pCell->state.neighbors.size() ; n++ )
	{
		Cell* pC = pCell->state.neighbors[n]; 
		if( pC->type_name == "target")
		{
			// if found, attach as spring
			// and rapidly transform back to axon 
			attach_cells_as_spring( pCell, pC ); 
			set_single_behavior( pCell , "transform to axon" , 10); 

			// also, turn off the target
			#pragma omp critical 
			{ set_single_behavior( pC, "chemoattractant secretion", 0); }
		}
	}

	return;
}

Cell* clone_cell( Cell* pSource, std::vector<double>& direction )
{
	static Cell_Definition* pCD_gcone = find_cell_definition( "gcone"); 

	Cell* pNew = create_cell( *pCD_gcone ); 

	std::vector<double> x = pSource->position; 
	double distance = pSource->phenotype.geometry.radius + pNew->phenotype.geometry.radius; 
	direction *= distance; 
	x += direction; 
	pNew->assign_position( x );

	return pNew; 
}

void axon_custom( Cell* pCell , Phenotype& phenotype , double dt )
{
   custom_helper(pCell,phenotype,dt);

   bool movable = (bool) get_single_behavior( pCell, "is movable");

   if(pCell->state.attached_cells.size() > 0)
   {
      if(movable == true)
      {
         phenotype.motility.migration_speed = 0;
      }
   }

   // check for intersection with starting neuron
   if(pCell->state.spring_attachments.size() == 0)
   {
      for( int n=0; n < pCell->state.neighbors.size() ; n++ )
      {
         Cell* pC = pCell->state.neighbors[n];
         if( pC->type_name == "start")
         {
            // if found, attach as spring
            attach_cells_as_spring( pCell, pC );
            // also, turn off the starter's recruiter signal
            #pragma omp critical
            { set_single_behavior( pC, "recruiter secretion", 0); }
         }
      }
   }
   return;
}

/* keep-alive signaling */
void alive_signaling( Cell* pCell, Phenotype& phenotype, double dt )
{
   double r_create = get_single_signal( pCell , "custom:alive_creation");
   double r_decay = get_single_signal( pCell, "custom:alive_decay");

   double alive = get_single_signal( pCell , "custom:alive_signal");
   // std::cout << alive << std::endl;
   double d_alive = dt * r_create * ( 1.0 - alive ) - dt * r_decay * alive;
   alive += d_alive;
   // std::cout << "\t" << d_alive << " " << alive << std::endl;

   set_single_behavior( pCell, "custom:alive_signal" , alive );

   return;
}

void setup_start( void )
{
   Cell_Definition* pCD = find_cell_definition( "start");

   pCD->phenotype.mechanics.maximum_number_of_attachments = 2;
   pCD->phenotype.mechanics.attachment_rate = 10;
   pCD->phenotype.mechanics.attachment_elastic_constant = 0.01;

   pCD->functions.update_phenotype = start_phenotype;
   pCD->functions.contact_function = flux_contact;
   pCD->functions.custom_cell_rule = custom_helper;

   return;
}

void start_phenotype( Cell* pCell , Phenotype& phenotype , double dt )
{
   /* alive signaling */
   alive_signaling(pCell,phenotype,dt);

   int n_adhesions = pCell->state.attached_cells.size();

   set_single_behavior( pCell, "custom:n_adhesions" , n_adhesions );

   return;
}

void setup_inhib( void ){

   Cell_Definition* pCD = find_cell_definition( "inhib");
   pCD->functions.update_phenotype = inhib_phenotype;
}

void inhib_phenotype( Cell* pCell , Phenotype& phenotype , double dt ){

   double c = get_single_signal( pCell, "inhibitory");
   double secrete_thresh = 0.01;
   double chemotaxis_thresh = 0.001;
   if (c > chemotaxis_thresh){
      pCell->phenotype.motility.migration_speed = 5;
   }
   if ( c > secrete_thresh){
      // double c_out = c * 0.8;

      set_single_behavior( pCell, "inhibitory secretion", 10);
   }

}

void sever_axon(double dt)
{
   static bool done = false;
   static double t = 0;
   static double activation_time = 6000;
   static double new_time = 6100;
   static bool really_done = false;
   static int severed_id;
   int n;

   t += dt;
   if(t>=activation_time && done==false)
   {
      //std::cout << "entering first if" << std::endl;
      //instantiate prng
      std::mt19937 engine{ static_cast<unsigned int>( 
         std::chrono::high_resolution_clock::now().time_since_epoch().count()) }; 
      std::uniform_real_distribution<double> dist{0.0, 1.0};
      std::vector<int> axon_indices;
      //find all the axons
      //std::cout << "Check all cells" << std::endl;
      int n_last = (*all_cells).size();
      for(n=0 ; n < n_last ; n++ )
      {
         Cell* pC = (*all_cells)[n];
         //std::cout << pC->ID << std::endl;
         if(pC->type_name == "axon" && pC->state.spring_attachments.size() == 2)
            axon_indices.push_back(n);
      }
      //iterate through all the axons
      //std::cout << "run through the axons" << std::endl;
      for(n=1; n < axon_indices.size()-1; n++)
      {
         double p = dist(engine);
         int axon_index = axon_indices[n];
         Cell* pC = (*all_cells)[axon_index];
         //std::cout << pC->ID << std::endl;
         if(p <= 0.075)
         {
            //std::cout << pC->state.spring_attachments.size() << std::endl;
            
            Cell* pLeft = pC->state.spring_attachments[0];
            Cell* pRight = pC->state.spring_attachments[1];
            detach_cells_as_spring(pC, pLeft);
            detach_cells_as_spring(pC, pRight);
            set_single_behavior( pC , "transform to inhib" , 10);
            //std::cout << pC->type_name << std::endl;
            //std::cout << pC->state.spring_attachments.size() << std::endl;
            //std::cout << pLeft->state.spring_attachments.size() << std::endl;
            //std::cout << pRight->state.spring_attachments.size() << std::endl;
            severed_id = pC->ID;
            //std::cout << pC->ID << std::endl;
            done = true;
            return;
         }

        /*Cell* pPrev = (*all_cells)[n-1];
         if( pC->type_name == "axon" && pPrev->type_name == "axon" && p <= 0.15)
         {
            detach_cells_as_spring(pC, pPrev);
            set_single_behavior( pC , "transform to inhib" , 10);
            done = true;
         } */
      }
   }
   if(t >= new_time && done == true && really_done == false)
   {
      //std::cout << "entering second if" << std::endl;
      int n_last = (*all_cells).size();
      for(n=0 ; n < n_last ; n++ )
      {
         Cell* pC = (*all_cells)[n];
         
         //std::cout << pC->ID << std::endl;
         if(pC->ID == severed_id)
         {
            //std::cout << "found the severed guy" << std::endl;
            //std::cout << pC->type_name << std::endl;
            set_single_behavior( pC, "inhibitory secretion", 100);
            really_done = true;
            return;
         }
      }
   }


   return;
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
		if( pC->type_name == "gcone" )
		{
			// extension direction 
			std::vector<double> x = pC->position; 
			std::vector<double> direction = microenvironment.nearest_gradient_vector(x)[n_chemo];
			std::vector<double> direction_inhib = microenvironment.nearest_gradient_vector(x)[n_inhib]; 
			normalize( &direction ); 
			normalize( &direction_inhib);

			// std::vector<double> rand = UniformOnUnitCircle(); 
			double bias = 0.80;
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

				// force fast conversion of the prior gcone 
				set_single_behavior( pC, "transform to axon" , 10); 
				}
			}
		}

	}	

	t_next += interval; 
	return; 
}
