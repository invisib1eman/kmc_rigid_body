//NANOROD: molecule.h Molecule Class (Revision Date:Oct 27, 2023)
#ifndef _MOLECULE_H
#define _MOLECULE_H
#include "header.h"
#include "xyz.h"
#include "quarternion.h"
#include "hbond.h"
#include "grid.h"
class Molecule
{
public:
    int MOL_ID;//molecule ID
    int N_VER; //Number of vertices
    int gID;//grid id
    int AID;//in which aggregate
    int AsubID;//id in the aggregate
    list<int> adjList;//AsubID of the neighbors, basically infer from hbond_list
    vector<XYZ> ver; //Coordinates of vertices
    vector<char> vertype;//type of vertices: A, B, C, D, I(inactive)
    int nbonds;//number of bonds
    vector<XYZ> aorigin; 
    
        
    
    vector<hbond> hbond_list;//list of hbonded neighbors and information
    XYZ centre; //Coordinates of centre
    quarternion orientation;
    Molecule() //HARD CODE BASED ON VERTEX INFORMATION,unit length=1nm,arm length=1.1nm
    {
        N_VER=6;
        nbonds=0;
        gID=0;
        centre.set(0.0,0.0,0.0);
        vertype.push_back('A');
        vertype.push_back('A');
        vertype.push_back('A');
        vertype.push_back('A');
        vertype.push_back('A');
        vertype.push_back('A');
        ver.push_back(XYZ(1.0,0.0,0.47));
        ver.push_back(XYZ(1.0,0.0,-0.47));
        ver.push_back(XYZ(-0.5,0.8660254,0.47));
        ver.push_back(XYZ(-0.5,0.8660254,-0.47));
        ver.push_back(XYZ(-0.5,-0.8660254,0.47));
        ver.push_back(XYZ(-0.5,-0.8660254,-0.47));
        aorigin.push_back(XYZ(1.0,0.0,0.47));
        aorigin.push_back(XYZ(1.0,0.0,-0.47));
        aorigin.push_back(XYZ(-0.5,0.866,0.47));
        aorigin.push_back(XYZ(-0.5,0.866,-0.47));
        aorigin.push_back(XYZ(-0.5,-0.866,0.47));
        aorigin.push_back(XYZ(-0.5,-0.866,-0.47));
        orientation=quarternion(1,0,0,0);
    }

        
    
    void UpdateVertices()
    {
        
        
        for(int i=0;i<N_VER;i++)
        {
            ver[i]=quarterrotation(aorigin[i],orientation)+centre;
        }
    }
};
#endif

