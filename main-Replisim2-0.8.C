#include <iostream>
#include <string.h>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <math.h>

#define PI 3.141592536

// commented out to disable profile integration from fibers
// #define WANGINT 1
#define WANG 1

// This project:
// do chr21 only with core origins
// Impose limit on active ORIs after t(origin) has beed randomly assigned
// 1 cell
// 1000 cells with randomized delays at specific times
// 1000 identical cells async

static bool async = false;
// static bool limit_origins_per_chr = true;
static bool limit_origins_per_chr = false;

// removed Wang data integration


// Simulate fibers with fluorescence
static double integration_time_step = 1. ;   // minutes

#ifdef WANG
class track // class for Wang data 
{
        public:
                unsigned int chr_id;
                int direction; // 0 is left, 1 is right
                long position;
                long parent_from;
                long parent_to;
                long parent_n;
                float time;
        // target function
        // requires reference global t_from_nearest or t_from_nearest_10.
                double delta_t();
};
#endif


// This project: time(origin) is variable.
//  * Origin didn't fire _because_ a fork reached it first
//  		( p(fire) per minute ).
//  		Can competences be derived from this model?
// vs
// ** Origin didn't fire because it was disabled a priori

static bool sqrt_tracks = false;
// static bool logarithmic_tracks = true;
static bool logarithmic_tracks = false;
static bool all_one_tracks = true;
static bool detrend_tracks = true;



static bool constant_vels = true;

static int hist_bins = 100000;

static int ori_integration_step = 1000;
static int profile_from_tracks_step = ori_integration_step;
// static int profile_from_tracks_step = 100;

//static std::string out_prefix = "chr21_async_10cells";
static std::string out_prefix = "chr21_CORE_flat_dist_30_minutes_100_cells";

static std::string oris_file = "Core_Origins_name_1_23_sorted_culled.bed";




static std::vector < long > chrl; // chromosome_len
static std::vector < long > chr_maxori; // chromosome_max origins
static std::vector < long > all_fibers; // fibers
static std::vector < long > all_fibers_joint; // fibers
static std::vector < long > met_fibers; // fibers
static std::vector < long > frk_fibers; // fibers
static std::vector <std::vector < long > > time_hists_chr; // partial histograms
static std::vector < long > time_hists; // final histograms
static std::vector < int > chrn; // chromosome_name
static std::vector <std::vector < long > > all_oris;
static std::vector <std::vector < float > > ttc;
static std::vector <std::vector < float > > velocities;
static std::vector <std::vector < double > > t_from_zero;
static std::vector <std::vector < float > > t_from_nearest;
static std::vector <std::vector < int > > id_of_nearest;
static std::vector <std::vector < float > > t_from_nearest_10;
static std::vector <std::vector < int > > id_of_nearest_10;
static std::vector <std::vector < int > > id_of_left_ori;

static std::vector <std::vector < int > > n_times;;
static std::vector <std::vector < double > > a_times;;
static std::vector <std::vector < double > > s_times;;


static std::vector <std::vector <std::vector < double > > > ori_delay;

static std::vector <std::vector < float > > profile_from_tracks;  // profile generated from Wang data
static std::vector <std::vector < float > > profile_from_tracks_differential;  // differential profile for Wang data


// for correlation
static std::vector <double> c_sxy;
static std::vector <double> c_sxx;
static std::vector <double> c_syy;
static std::vector <double> c_sx;
static std::vector <double> c_sy;
static std::vector <long> c_n;
static std::vector <long> c_n_total;

int main(int argc, char *argv[])
{

	if ( argc < 2 ) {
		std::cout <<
			"Usage " << 
			argv[0] << " " <<
			"1.7 0.5 30 27 0.0001 300 3000 " <<
			"chr21_CORE_flat_dist_30_minutes_100_cells " << 
			"Core_Origins_name_1_23_sorted_culled.bed " << 
			std::endl;
		exit(7);
	};

	long len, maxori;
	long orip, orip1, orip2;
	int name;
	name = 1;
	// int n_similar_cells = 10;
	// int n_similar_cells = 100;
	// int n_similar_cells = 10;
	int n_similar_cells = 3;

	// PARAMETERS
	float avg_speed = 1.7;
	avg_speed = atof(argv[1]);
	float sig_speed = 0.5;
	sig_speed = atof(argv[2]);
	float t_fibers = 30.0;
	t_fibers = atof(argv[3]);
	float max_delay = 28.;
	max_delay = atof(argv[4]);
	float p_dead = 0.;
	p_dead = atof(argv[5]);   // prob that dead origin
	float min_length = 100.;
	min_length = atof(argv[6]);   // power spectrum
	float rand_length = 3000.;
	rand_length = atof(argv[7]);   // power spectrum
	out_prefix = argv[8];
	oris_file = argv[9];

	static int v_step = 10;

	std::string out_percell_fibers = out_prefix + ".percell_fibers";
	std::string out_percell_trep = out_prefix + ".percell_trep";
	std::string out_percell_asrep = out_prefix + ".percell_asrep";

	// open out files
	std::ofstream out_time_bycell (out_percell_trep, std::ofstream::out);
	std::ofstream out_time_av_sig (out_percell_asrep, std::ofstream::out);
	std::ofstream out_fibers_bycell (out_percell_fibers, std::ofstream::out);

	int nchr = 0;
	int qqqqq;
	name = 0;
	// std::ifstream ifchrl ("Chromosome_lengths.txt", std::ifstream::in);
	// std::ifstream ifchrl ("hg19.chrom.sizes_23", std::ifstream::in);
	std::ifstream ifchrl ("hg19.chrom.sizes_23_with_maxoris", std::ifstream::in);
	while ( ifchrl >> qqqqq >> len >> maxori )
	{
		chrl.push_back(len);
		chr_maxori.push_back(maxori);
		chrn.push_back(name);
		++name;
		++nchr;
	}

	a_times.resize(nchr);
	s_times.resize(nchr);
	n_times.resize(nchr);

	all_oris.resize(nchr);
	id_of_left_ori.resize(nchr);
	ori_delay.resize(n_similar_cells);
	for (int c = 0; c < n_similar_cells ; ++c )
	{
		ori_delay[c].resize(nchr);
	}
	ifchrl.close();
	//std::ifstream iforis ("allOrisdNTPexp.wig", std::ifstream::in);
	long last_orip = 0;
	int  last_name = 999999;
	// for random delays
	// srand(17);
	srand(time(0));
	//
	//

	c_sxy.resize(nchr); c_sxx.resize(nchr); c_syy.resize(nchr); c_sx.resize(nchr); c_sy.resize(nchr); c_n.resize(nchr); c_n_total.resize(nchr);

//
//
//	Core Origins hg19
//
	// std::ifstream iforis ("Core_Origins.bed", std::ifstream::in);
	std::cerr << "READING ORIS " << std::endl;
	//this is good for core oris
	//std::ifstream iforis ("Core_Origins_name_1_23_sorted.bed", std::ifstream::in);
	//this is good for core oris only chr21
	//std::ifstream iforis ("Core_Origins_name_1_23_sorted_culled.bed", std::ifstream::in);
	std::ifstream iforis (oris_file.c_str(), std::ifstream::in);
	//     std::ifstream iforis ("open-chromatin-dnase-sorted-hg19.bed", std::ifstream::in);
	//double d1, d2;
	// while ( iforis >> name >> orip1 >> orip2  >> d1 >> d2 )
// 1	211812010	211812470	Valid_4770	1000	.	211812010	211812470
	std::string d1, d2, d3, d4, d5;
	int jjj;
	// while ( iforis >> name >> orip1 >> orip2  >> d1 >> d2 >> d3 >> d4 >> d5 )
	while ( iforis >> name >> orip1 >> orip2  >> d1 >> d2 )
	{
		std::cerr << "    ori " << name << " " << orip1 <<  " " <<orip2  <<  " " <<d1 <<  " " <<d2 <<  " " << ++ jjj << "\n"; 
		orip = int ( ( orip1+orip2 ) / 2 );
		all_oris[name-1].push_back(orip);
		if ( last_name == name ) 
		{
			// std::cerr << "NEW_LAST\t" << name << "\t"  << orip << " " <<  last_orip << std::endl;
			assert (orip >= last_orip ) ;
		}; // make sure origins are sorted
#pragma omp parallel for
		for (int c = 0; c < n_similar_cells ; ++c )
		{
			double delay;
			// delay = 0;// Zero delay 
			// delay = max_delay*(double)rand()/RAND_MAX - max_delay*0.5 ; // Flat delay , centered 
			delay = max_delay*(double)rand()/RAND_MAX ; // Flat delay 
			if ( ((double)rand()  / RAND_MAX ) < p_dead ) { delay = 1000000;};
			//
			//ori_delay[c][name-1].push_back( 0 );  
			ori_delay[c][name-1].push_back( delay ); 
		};
#pragma omp single
		last_name = name;
		last_orip = orip;
	};
	iforis.close();
	std::cerr << "READING ORIS ... DONE " << std::endl;
	
// disable origins other than max_origins.
	if (limit_origins_per_chr)
	{
	for (int c = 0; c < n_similar_cells ; ++c )
	{
		for ( int i=0; i<nchr; ++i )
		{
			long int nco = all_oris[i].size();
			std::vector<int> active_oris( nco , 0);
			long activated_oris = 0;
			do {
				// pick random
				int xoro = std::rand() / ((RAND_MAX + 1u) / nco );
				// check if 0 and { set to 1 and increment counter }
				if ( active_oris[xoro] == 0 )
				{
					++activated_oris;
					active_oris[xoro] = 1;
				};
			} 
			while ( activated_oris < chr_maxori[i] );
			for (int ior = 0; ior < nco; ++ior )
			{
				if ( ! active_oris[ior] ) { ori_delay[c][i][ior] = 2222222; };
			}
		};
	};
	};

// test output ori delays
	for (int c = 0; c < n_similar_cells ; ++c )
	{
		for ( int i=0; i<nchr; ++i )
		{
			std::cerr << c  << "\tNORIS\t" << i << "\t" << ori_delay[c][i].size() << "\n";
			// for ( auto del : ori_delay[c][i] ) { std::cerr << del << "\t"; } std::cerr <<  std::endl;
		};
	};

//
//	Late Origins DNAse HS  // later
//

	// std::cerr << "TRAP: \n";

	// make all vectors
	velocities.resize(nchr);
	ttc.resize(nchr);
	t_from_zero.resize(nchr);
	t_from_nearest.resize(nchr);
	id_of_nearest.resize(nchr);
	t_from_nearest_10.resize(nchr);
	id_of_nearest_10.resize(nchr);
	time_hists.resize(hist_bins); for (long q=0 ; q < hist_bins ; ++ q )  {time_hists[q] = 0 ;};
	time_hists_chr.resize(nchr);
	profile_from_tracks.resize(nchr); // profile from Wang tracks
	profile_from_tracks_differential.resize(nchr); // profile from Wang tracks
#pragma omp parallel for
	for ( int i=0; i<nchr; ++i )
	{
		long reduced_l = 1 + (long) ( chrl[i] / ori_integration_step ); // reduced  chromosome length
#pragma omp critical
		{
			std::cerr << "ALLOCATING : " << i << " with nchr=" << nchr << " \n";
			std::cerr << "ALLOCATING : " << i << " reduced_l=" << reduced_l << " \n";
		}
		a_times[i].resize(reduced_l); std::fill(a_times[i].begin(), a_times[i].end(), 0.);
		s_times[i].resize(reduced_l); std::fill(s_times[i].begin(), s_times[i].end(), 0.);
		n_times[i].resize(reduced_l); std::fill(n_times[i].begin(), n_times[i].end(), 0);
		velocities[i].resize(chrl[i]);
		ttc[i].resize(chrl[i]);
		time_hists_chr[i].resize(hist_bins); for (long q=0 ; q < hist_bins ; ++q ) {time_hists_chr[i][q] = 0 ;};
		t_from_zero[i].resize(chrl[i]);
		id_of_left_ori[i].resize(chrl[i]);
		// t_from_nearest[i].resize(chrl[i]);
		// id_of_nearest[i].resize(chrl[i]);
		t_from_nearest_10[i].resize(chrl[i]);
		id_of_nearest_10[i].resize(chrl[i]);
		// prof _ from _ trx
		profile_from_tracks[i].resize(chrl[i] / profile_from_tracks_step + 1 ) ;
		profile_from_tracks_differential[i].resize(chrl[i] / profile_from_tracks_step + 1 ) ;
		std::fill(profile_from_tracks[i].begin(), profile_from_tracks[i].end(), 0);
		std::fill(profile_from_tracks_differential[i].begin(), profile_from_tracks_differential[i].end(), 0);
	};
	std::cerr << "ALLOCATED ALL: \n";

#ifdef WANG
// read Wang data trx
        std::cerr << "Reading Wang " << "list_of_files.txt" << std::endl;
	int sample_id = 0;
	std::ifstream ifiles ("list_of_files.txt", std::ifstream::in);	
	std::vector < track > wang_data ;
	std::vector < float > wang_times ;
	std::vector < std::string > wang_files ;
	std::string fname;
	double stime;
	std::ifstream wang_file;
	//
	while ( ifiles >> fname >> stime )
	{
		//
		// NEED to add STIME for more times in one model
		//
		std::cerr << "reading from " << fname << std::endl;
                wang_files.push_back(fname);
                long start_bp;
                long end_bp;
                double strength;
                char direction;
                std::string track_id;
                std::string chr;
                int chrn;
                size_t pos = 3;
                wang_file.open (fname.c_str());
                while ( wang_file  >> chr >> start_bp >> end_bp >> track_id >> strength >> direction )
                {
                        std::string schr = chr.substr(3);
                        chrn = stoi ( schr );
			if ( chrn > 23 ) { continue ; } ;  
                        int chr_i = chrn -  1;
                        // std::cerr << "WANG\t" << chr << "\t" << chrn << std::endl;
                        // TUTAJ
			int i_strength;
			i_strength = 1;
			// std::cerr << "STR " <<  strength << "\n";
			//if ( all_one_tracks ) 
			//{
				//i_strength = 1;
			//}
			//else if ( sqrt_tracks ) 
			//{ 
				//i_strength =  int( .0001 + sqrt(strength) ) ; 
			//} 
			//else if ( logarithmic_tracks )
			//{
				//i_strength = int( 1.0001 + log(strength) ) ;
			//}
			//else
			//{
				//i_strength = int ( strength + 0.0001 );
			//};	
			//
			long int start_idx = int ( start_bp / profile_from_tracks_step );
			long int end_idx = int ( 1 + end_bp / profile_from_tracks_step );
			//
			//
			for (long int position = start_idx; position <= end_idx; ++position )
			{
				//std::cerr << "TRYING " << chr << "  (" << chr_i << ")  : " << position << std::endl ;
				profile_from_tracks[ chr_i ][ position ] += 1;
				if ( position % 1000 == 0 ) { std::cerr << "profile_from_tracks\tchr\t" << chr_i << "\tpos\t" << position << "\t" <<  profile_from_tracks[ chr_i ][ position ] << "\t" << fname.c_str() << std::endl; }
			};
			//
			//
                        // int step = ( end_bp - start_bp ) / i_strength;
                        // for (int ix = 0 ; ix < strength ; ++ ix )
                        // {
                                // track t;
                                // t.chr_id = chr_i;
                                // t.direction = ( direction == '+' ) ? 1 : 0 ;
                                // t.position = ix * step + start_bp + (int) step/2 ;
                                // t.parent_from = start_bp;
                                // t.parent_to   = end_bp;
                                // t.parent_n    = strength;
                                // t.time = stime;
                                // wang_data.push_back(t);
// 
                        // };

                };
                wang_file.close();
                std::cerr << "closing wang " << fname << std::endl;
	}
	ifiles.close();
// endif WANG
#endif

#pragma omp parallel for
	for ( int i=0; i<nchr; ++i )
	{
		int last_ori;
		long curr=0;
		for ( int q=0; q<all_oris[i].size(); ++q )
		{
			last_ori = q;
			while (curr < all_oris[i][q])
			{
				id_of_left_ori[i][curr] = q-1 ;   // -1 means we're left of zeroeth ori
				++curr;
			};
		};
		while (curr < chrl[i])
		{
			id_of_left_ori[i][curr] = last_ori ;   // between last ori ( all_oris[i].size()-1 ) and + end of chromosome.
			++curr;
		};

	};
 
	// now velocities
	double velscale = 1.0;
	std::cerr << "vels: \n";
	if (constant_vels)
	{
		std::cerr << "WARNING: constant vels: \n";
#pragma omp parallel for
		for ( int i=0; i<nchr; ++i )
		{
			for ( long p=0; p<chrl[i]; p++ )
			{
				velocities[i][p] = avg_speed;
				ttc[i][p] = 0.001 / velocities[i][p];
			}
		}
	} 
	else
	{
//
	double vel = avg_speed;
	double vsd = sig_speed;
	double vmin = 0.1;
	// power
	int n_freqs = 5;
	std::vector < double > freqs; freqs.resize(n_freqs);
	std::vector < double > amplitudes; amplitudes.resize(n_freqs);
	std::vector < double > phases; phases.resize(n_freqs);
	for ( int q=0; q<n_freqs; ++q )
	{
		// double length = 1000. + 30000. * ( (double)rand() ) / RAND_MAX;
		double length = min_length + rand_length * ( (double)rand() ) / RAND_MAX;
		// double length = 100. + 3000. * ( (double)rand() ) / RAND_MAX;
		freqs[q] = 2. * PI / length;
		amplitudes[q] = ( (double)rand() ) / RAND_MAX;
		phases[q] = 2. * PI * ( (double)rand() ) / RAND_MAX;
		std::cerr << q << "\t" 
			<< freqs[q] << "\t"
			<< amplitudes[q] << "\t"
			<< phases[q] << "\t" << std::endl;
	};
	double ssqvel = 0.;
	velscale = 0.;
	long allvel = 0;
#pragma omp parallel for
	for ( int i=0; i<nchr; ++i )
	{
		std::cerr << "Generate vels chr " << i+1 << std::endl;
		// for ( long p=0; p<chrl[i]; ++p )
		for ( long p=0; p<chrl[i]; p++ )
		{
			double lastvel;
			if ( p%v_step == 0 )
			{
				velocities[i][p] = 0;
				for ( int q=0; q<n_freqs; ++q )
				{
					velocities[i][p] += amplitudes[q]*sin(p*freqs[q]+phases[q]);
				};
				lastvel = velocities[i][p];
			}
			else
			{
				velocities[i][p] = lastvel;
			};
			++allvel;
			ssqvel+= velocities[i][p] * velocities[i][p];
			//std::cout << i << "\t" << p << "\t" << velocities[i][p] << std::endl;
		};
	};
	velscale = vsd / sqrt ( ssqvel / allvel );
#pragma omp parallel for
	for ( int i=0; i<nchr; ++i )
	{
		std::cerr << "Scaling chr " << i+1 << std::endl;
		for ( long p=0; p<chrl[i]; ++p )
		{
			velocities[i][p] *= velscale;
			velocities[i][p] += vel;
			if ( velocities[i][p] < vmin )
			{
				velocities[i][p] = vmin ;
			};
	 		// std::cout << i << "\t" << p << "\t" << velocities[i][p] << std::endl;
			ttc[i][p] = 0.001 / velocities[i][p];
		};
	}
	};
	for ( int i=0; i<nchr; ++i )
	{
		t_from_zero[i][0] = ttc[i][0];
		for ( long p=1; p<chrl[i]; ++p )
		{
			t_from_zero[i][p] =  t_from_zero[i][p-1] + ttc[i][p];
	 		// std::cout << "T_From_ZERO\t" << i << "\t" << p << "\t" << velocities[i][p] << "\t" << ttc[i][p] << "\t" << t_from_zero[i][p] << std::endl;
		};
	};
	// velocities generated
	// ttc generated
	std::cerr << " ttc generated " << std::endl; 
	
	

	for (int c = 0; c < n_similar_cells ; ++c )
	{
	std::cerr << "Doing CELL " << c << std::endl;
#pragma omp parallel for
	// for ( int i=0; i<nchr; ++i )
	for ( int i=21; i<22; ++i )
	{
		// create list of sorted ORIs
		// plus_sort  for t_from_0   + delay   for when we are LEFT  of ORI
		// minus_sort for t_from_end + delay   for when we are RIGHT of ORI
		//
		// REMOVED FOR NOW, 10 is OK
#if 0
		std::vector < std::pair < double, int > > plus_sort;
		std::vector < std::pair < double, int > > minus_sort;
		for ( int q=0; q<all_oris[i].size(); ++q ) // load ori vectors
		{
			double t_plus  = t_from_zero[i][ all_oris[i][q] ] + ori_delay[c][i][q] ;
			double t_minus = t_from_zero[i][chrl[i]-1] - t_from_zero[i][ all_oris[i][q] ] + ori_delay[c][i][q] ;
			plus_sort.push_back  (std::make_pair( t_plus,  q));
			minus_sort.push_back (std::make_pair( t_minus, q));
		};
		std::sort( plus_sort.begin(), plus_sort.end() );
		std::sort( minus_sort.begin(), minus_sort.end() );
		std::vector < int > pos_in_plus_sorted;
		std::vector < int > pos_in_minus_sorted;
		pos_in_plus_sorted.resize(all_oris[i].size());
		pos_in_minus_sorted.resize(all_oris[i].size());
		for (int q=0; q<all_oris[i].size(); ++q )
		{
			pos_in_plus_sorted[ plus_sort[q].second ] = q;
			pos_in_minus_sorted[ minus_sort[q].second ] = q;
		};
		for (int q=0; q<all_oris[i].size(); ++q )
		{
			// std::cerr << "RESORTED ORI chr " << i+1 << "\tq: " << q << "\tpos_in_plus_sorted: " << pos_in_plus_sorted[q] << std::endl;
		};
#endif
		//
#pragma omp critical
		{
		std::cerr << "nearest ORI chr " << i+1 << "\t" << all_oris[i].size() << std::endl;
		}
		//// for ( long p=0; p<chrl[i]; ++p )
		for ( long p=0; p<chrl[i]; p+= ori_integration_step )
		{
			long p_index = p/ori_integration_step;
			if ( p % 10000000 == 0 ) 
			{ 
#pragma omp critical
				{
					// std::cerr << "TRAPSF  C " << i << "  P " << p << std::endl ; 
				};
			};
			// the slow search
			//
			int min_o = 100000000;
			double min_t = 10000000000.;
#if 0
			for ( int q=0; q<all_oris[i].size(); ++q )
			{
				double dt = abs (t_from_zero[i][p] - t_from_zero[i][ all_oris[i][q] ] ) + ori_delay[c][i][q] ;
				if (dt < min_t)
				{
					min_t = dt;
					min_o = q;
				};
			};
			t_from_nearest[i][p] = min_t;
			id_of_nearest[i][p] = min_o;
#endif 
			//
			// the faster search
			//
			int min_o_10 = 100000000;
			double min_t_10 = 10000000000.;
			int qstart = id_of_left_ori[i][p] - 5;
			int qstop = id_of_left_ori[i][p] + 7;
			if ( qstart < 0 ) { qstart = 0 ; };
			if ( qstop > all_oris[i].size() ) { qstop = all_oris[i].size() ; };
			for ( int q=qstart; q<qstop; ++q )
			{
				double dt = abs (t_from_zero[i][p] - t_from_zero[i][ all_oris[i][q] ] ) + ori_delay[c][i][q] ;
				// std::cerr << "TFZERO" << c << " " << i << " " << p << " " << q << "\t" << t_from_zero[i][p] << "\t" << t_from_zero[i][ all_oris[i][q] ] << "\t ODELAY " << ori_delay[c][i][q]  << "\n";
				if (dt < min_t_10)
				{
					min_t_10 = dt;
					min_o_10 = q;
				};
			};
			t_from_nearest_10[i][p] = min_t_10;
			id_of_nearest_10[i][p] = min_o_10;
#if 0
			//
			// the fast search
			//
			int fmin_o = 100000000;
			double fmin_t = 10000000000.;
			// the left ori
			// the right ori
			int i_left  = id_of_left_ori[i][p];
			int i_left_in_sorted = pos_in_plus_sorted[i_left];
			int i_right = id_of_left_ori[i][p] + 1;
			int i_right_in_sorted = pos_in_plus_sorted[i_right];
			double dt_l = 100000000000.;
			int t_i_left;
			do
			{
				if (i_left_in_sorted) { --i_left_in_sorted; } // in plus-sorted
				t_i_left = plus_sort[i_left_in_sorted].second;	
			}
			while ( t_i_left > i_left );
			double fmin_t_left;
			double fmin_t_right;
			fmin_t_left  = abs (t_from_zero[i][p] - t_from_zero[i][ all_oris[i][i_left] ] ) + ori_delay[c][i][i_left] ;
			fmin_t_right = abs (t_from_zero[i][p] - t_from_zero[i][ all_oris[i][i_right] ] ) + ori_delay[c][i][i_right] ;
#endif 
			//
			//////dt_l = t_from_zero[i][p] - t_from_zero[i][ all_oris[i][q] ] 
			//
			// done both searches
			//
			// if ( id_of_nearest[i][p]  != id_of_nearest_10[i][p] )
			if ( 0 ) 
			{
#pragma omp critical
				{
	 		std::cout << "ALLBPS\t" << i << "\t" << p << "\t" << velocities[i][p] << "\t" << ttc[i][p] << "\t" << t_from_zero[i][p] << "\t" 
		//		<< "S\t" << t_from_nearest[i][p] << "\t" << id_of_nearest[i][p] << "\t"
				<< "F\t" << t_from_nearest_10[i][p] << "\t" << id_of_nearest_10[i][p] 
				<< std::endl;
				};
			};
			//
			// main bycell output 
#pragma omp critical
			{
			out_time_bycell 
				<<  c << "\t"
				<< i << "\t" << p << "\t" << t_from_nearest_10[i][p] 	
				<< std::endl;
			};
#pragma omp critical
			{
				++n_times[i][p_index];
				a_times[i][p_index] += t_from_nearest_10[i][p] ;
				s_times[i][p_index] += ( t_from_nearest_10[i][p] * t_from_nearest_10[i][p] );
				if ( p_index % 1000 == 0 ) { std::cerr << "MAKING_TIMES\t" << i << "\t" << p_index << "\t" << n_times[i][p_index] << "\t"  << a_times[i][p_index] << std::endl; };
			};
			//
			int bin = (long)(t_from_nearest_10[i][p] * 10 );
			if (bin < 0 ) { bin = 0;  };
			if (bin > hist_bins-1 ) { bin = hist_bins-1;  };
			time_hists_chr[i][bin] ++ ;
		}; // end loop over positions in chromosome
#pragma omp critical
		{
			std::cerr << "passed nearest ORI search  " << i << "\n";
		}

#ifdef FIBERS
	// FIBERS output
	//
	// check positions at t=t_fibers;

		// consecutive with t<t_fibers

		long last_off = 0;
		long last_on  = 0;
		for ( long p=1; p<chrl[i]-1; ++p )
		{
			if (
					( t_from_nearest_10[i][p] < t_fibers ) 
					&&
					(
						( t_from_nearest_10[i][p-1] > t_fibers )
						||
						( p == 1 )
					)
			   )
			{
				last_on = p;
			};
			if (
					( t_from_nearest_10[i][p] < t_fibers ) 
					&&
					(
						( t_from_nearest_10[i][p+1] > t_fibers )
						||
						( p == chrl[i]-2 )
					)
			   )
			{
				last_off = p;
				long length = abs ( last_off - last_on );
				all_fibers_joint.push_back(length);

			};
		};
		std::cerr << "passed all_fibers_joint.push_back " << i << "\n";
// #if 1
		for ( long p=0; p<chrl[i]-1; ++p )
		{
			// fork exists @ t=t_fibers
			if ( ( t_from_nearest_10[i][p] - t_fibers ) * ( t_from_nearest_10[i][p+1] - t_fibers ) < 0)
			{
				long length = abs ( p - all_oris[i][ id_of_nearest_10[i][p] ] );
				all_fibers.push_back(length);
				frk_fibers.push_back(length);
			};
		};
		std::cerr << "passed all_fibers_frk.push_back " << i << "\n";
		for ( long p=1; p<chrl[i]-1; ++p )
		{
			// forks met here at t < t_fibers
			if (
					( t_from_nearest_10[i][p] < t_fibers ) 
					&&
					( t_from_nearest_10[i][p] > t_from_nearest_10[i][p+1] )
					&&
					( t_from_nearest_10[i][p] > t_from_nearest_10[i][p-1] )
			   )
			{
				long length = abs ( all_oris[i][ id_of_nearest_10[i][p+1] ] - all_oris[i][ id_of_nearest_10[i][p-1] ] );
				all_fibers.push_back(length);
				met_fibers.push_back(length);
			};
		};
		std::cerr << "passed all_fibers_met.push_back " << i << "\n";
#endif

	}; // end loop over chromosomes 

	// HERE CORR
	
	}; // end loop  over similar cells 


#pragma omp single

	// output time profiles 
	// for ( int i=0; i<nchr; ++i )   
	// CHR21
	for ( int i=20; i<22; ++i )
	{
		for ( long pl = 0 ; pl < ( (chrl[i]-1) / ori_integration_step ) ; ++pl )
		{
			if ( pl % 1000 == 0 ) { std::cerr << "CHKING_TIMES\t" << i << "\t" << pl << "\t" << n_times[i][pl] << "\t"  << a_times[i][pl] << std::endl;};
			double t_avg = a_times[i][pl] /  n_times[i][pl] ;
			double t_sig = sqrt ( 
				s_times[i][pl] /  n_times[i][pl]  -  pow ( ( a_times[i][pl] / n_times[i][pl] ) , 2 ) 
				) ;
			out_time_av_sig << i+1 << "\t" << pl*ori_integration_step << "\t" << t_avg << "\t" << t_sig << "\t" << profile_from_tracks[ i ][ pl ] << std::endl;
		}
	}


#pragma omp single

	std::cerr << "ready to acc hists" << "\n";

	for ( int i=0; i<nchr; ++i )
	{
		for (long q =0; q<hist_bins; ++q )
		{
			time_hists[q] += time_hists_chr[i][q];
		};
	};
	for (long q =0; q<hist_bins; ++q )
	{
		// std::cout << "TH\t" << q/10. << "\t" << time_hists[q]  << std::endl;
	};

	for(long i : all_fibers )
	{
		std::cout << "FB\t" << i << std::endl;
	};
 	       for(long i : met_fibers )
        {
                std::cout << "FM\t" << i << std::endl;
        };
        for(long i : frk_fibers )
        {
                std::cout << "FF\t" << i << std::endl;
        };
        for(long i : all_fibers_joint )
        {
                std::cout << "FJ\t" << i << std::endl;
        };


#ifdef WANGINT
	// NOW CORR 

	for ( int qt=1; qt<300; ++qt )
	{
	double t_max = 2 * qt;

#pragma omp parallel for
	for ( int i=0; i<nchr; ++i )
	{
		c_sxy[i]=0; c_sxx[i]=0; c_syy[i]=0; c_sx[i]=0; c_sy[i]=0; c_n[i]=0; c_n_total[i]=0;
		//for ( long p = 0 ; p < chrl[i]; p += profile_from_tracks_step )
		for ( long itr = 1; itr < chrl[i] / profile_from_tracks_step ; itr +=  1  )
		{
			long p = itr * profile_from_tracks_step;
			// replication state from tracks
			double x = profile_from_tracks[i][itr] ;
			// replication state from replisim
			double t = t_from_nearest_10[i][p];
			double y = ( t < t_max / 10. ) ? 1 : 0 ;
			++ c_n_total[i];
			if (t > t_max) 
			{
				continue;
			};
			++ c_n[i];
			c_sxx[i] += x*x ;
			c_sxy[i] += x*y ;
			c_syy[i] += y*y ;
			c_sy[i] += y ;
			c_sx[i] += x ;
		}
#pragma omp critical
		{
			// std::cerr << "tmax: " << t_max << "\tchr " << i << "\t n= " << c_n[i] << std::endl;
		};
	}
	double ca_sxy=0; double ca_sxx=0; double ca_syy=0; double ca_sx=0; double ca_sy=0; long ca_n=0; long ca_n_total=0;
#pragma omp single
	for ( int i=0; i<nchr; ++i )
	{
		ca_sxx += c_sxx[i];
		ca_sxy += c_sxy[i];
		ca_syy += c_syy[i];
		ca_sx += c_sx[i];
		ca_sy += c_sy[i];
		ca_n  += c_n[i];
		ca_n_total  += c_n_total[i];
	};
	double corr = 
		( ca_n * ca_sxy - ca_sx * ca_sy ) 
		/ 
		( 
		 sqrt ( ca_n * ca_sxx - ca_sx * ca_sx )
			* 
		 sqrt ( ca_n * ca_syy - ca_sy * ca_sy )
		);

	std::cout << "t_max:\t" << t_max << "\tcorr:\t" << corr << "\tfrom n=" << ca_n << " ("  << int ( (100. * (double)ca_n )/ ca_n_total ) <<  "%)" << std::endl;
	};
#endif //WANGINT

	out_time_bycell.close();
	out_time_av_sig.close();
	out_fibers_bycell.close();

};

