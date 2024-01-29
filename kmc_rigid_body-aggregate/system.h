//NANOROD: system.h System Class (Revision Date: Oct 27,2023)
#ifndef _SYSTEM_H
#define _SYSTEM_H
#include "header.h"
#include "xyz.h"
#include "quarternion.h"
#include "molecule.h"
#include "utils.h"
#include "grid.h"
#include "aggregate.h"
class System
{
  public:
    //define lists/vectors of the systems
    vector<Molecule> M; //List of molecules
    vector<Grid> G;//grid list, should be constant throughout simulation
    list<hbond> H;//hbond list
    vector<Aggregate> Ag;//list of aggregates
    //define constants of systems
    int NGRID,NGRID3;//NGRID=#cell in one direction, NGRID3=total # of cells of grid, NCELL=#particles per cell
    double NCELL;
    double GRIDL;//Grid length in one direction
    const gsl_rng_type * gsl_T;
    gsl_rng * gsl_r;
    string Description;
    int NMOL; //Number of molecules
    double coresize=0.2;
    double arm_L=1.1;
    double bond_length=0.3;//minimum of fene potential, see energy.h
    double bond_extension=0.2;//max init energy is 6.4kT,see energy.h
    double searchl2_bond=pow(bond_length+bond_extension,2);//the length range of forming a bond
    double cm_L=arm_L*2+bond_length;
    double search2_cm=pow(cm_L+bond_extension,2);//the cm distance when it is possible to form a bond
    double search_angle=0.3;//max init energy is 13.5kT at Kalpha=300
    double search_xhi=0.3;
    double L; //Length of box
    int GSL_SEED; //Seed of random number generator
    int nsweep; //Number of MC sweeps
    double deltat; //Timestep
    double MCstep; //Step size of translation
    double MCstep_rotation;//Step size of rotation
    double E_1=10;//hbond dis enthalpy
    double free_bond_freeenergy=-1;//free bond entropy
    double omega_B=0.001;//arrhenius prefactor
    double omega_T=0.1;//frequency of change types
    void ReadInput(int argc, char *argv[])
    {
        double total_time;
        
        options_description desc("Usage:\nNANOROD <options>");
    
        desc.add_options()
        ("help,h", "print usage message")
        ("NGRID,G",value<int>(&NGRID)->default_value(20),"grids(default 20)")
        ("NMOL,N", value<int>(&NMOL)->default_value(400), "#molecules (default 400)")
        ("box_length,L", value<double>(&L)->default_value(40.0), "length of box (default 40.0)")
        ("time,s", value<double>(&total_time)->default_value(500.0), "total time in tau_0 units (default 500.0)")
        ("MCstep,m", value<double>(&MCstep)->default_value(0.1), "MC step size (default 0.1)")
        ("GSL_SEED,g", value<int>(&GSL_SEED)->default_value(10), "seed for the RNG (default 10)")
        ("Description,D", value<string>(&Description)->default_value("nanorod"), "Description (default nanorod)");

        MCstep_rotation=sqrt(3/4)*MCstep;
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);
        
        if (vm.count("help"))
        {
            cout << desc << "\n";
            exit(1);
        }
        
        gsl_rng_env_setup();
          
        gsl_T = gsl_rng_default;
        gsl_r = gsl_rng_alloc(gsl_T);
        gsl_rng_set(gsl_r, GSL_SEED);
        
        deltat=1.0/12.0*MCstep*MCstep;
        nsweep=int(ceil(total_time/deltat));
        NGRID3=NGRID*NGRID*NGRID;
        
    try
    {
      G.reserve(NGRID3);
    }
    catch (int e)
    {
	cout<<"Memory issues in cell list allocation.. exiting"<<endl;
	exit(1);
    }
        GRIDL=L/NGRID;
        if(GRIDL<cm_L*0.5)
        {
            cout<<"Error: Grid size too small"<<endl;
            exit(1);
        }
        NCELL=NMOL/NGRID3;
    }
    
    void Create();
    void WriteMol2(int timestep);
    void WriteDump(int timestep);
    void WriteBond(int timestep);
    void WriteOrientation(int timestep);
    void WriteGrid(int timestep);
    void WriteAggregate(int timestep);
    void UpdateGrid();
};
#endif
