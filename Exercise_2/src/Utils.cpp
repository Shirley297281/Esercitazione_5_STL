#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <iomanip>

namespace PolygonalLibrary {
bool ImportMesh(PolygonalMesh& mesh)
{
    vector<double> edges_lenght(mesh.NumberCell1D);


    if(!ImportCell0Ds("Cell0Ds.csv",
                       mesh))
    {
        return false;
    }
    else
    {
        cout << "Cell0D markers -->" << endl;
        for(auto it = mesh.Cell0DMarkers.begin(); it != mesh.Cell0DMarkers.end(); it++)
        {
            cout << "marker:\t" << it -> first << "\t ids (points) for that marker:";
            for(const unsigned int id : it -> second)
                cout << "\t" << id;

            cout << endl;
        }
    }

    if(!ImportCell1Ds("Cell1Ds.csv",
                       mesh))
    {
        return false;
    }
    else
    {
        cout << "Cell1D markers-->" << endl;
        for(auto it = mesh.Cell1DMarkers.begin(); it != mesh.Cell1DMarkers.end(); it++)
        {
            cout << "marker:\t" << it -> first << "\t ids (points) for that marker:";
            for(const unsigned int id : it -> second)
                cout << "\t" << id;

            cout << endl;
        }

        cout << " \n \t √ all markers have been stored correctly \n  "<<endl;

    }

    // Test : the edges of polygons have non-zero length

/// switch to italian

    /// 1. controllo su 1D che i vertici non siano gli stessi
    /// 2. cotrollo incorciato con 0D la x e la y che non abbiano differenza minore della tolleranza

    for(unsigned int c = 0; c < mesh.NumberCell1D; c++)
    {

        Vector2i vertices = mesh.Cell1DVertices[c];

        // in Cell0D, cerco il punto e vedo se quello meno l'altro estremo sono di lungh > tol
        const unsigned int origin = vertices[0];
        const unsigned int end = vertices[1];

        double lunghezza_lato = distance(origin, end, mesh);

        edges_lenght.push_back(lunghezza_lato);

        // Verifica se la lunghezza del lato è maggiore della tolleranza
        if (lunghezza_lato < mesh.tolDefault) {
            cerr <<"Zero leght of edge "<< c<<endl;
            return false;
        }

    }

    cout << " \t √ all edges of polygons have non-zero length \n  "<<endl;




    int little = 0;
    if(!ImportCell2Ds("Cell2Ds.csv",
                       mesh))
    {
        return false;
    }
    else
    {
        double area = 0.0;


        // Test: the area of the polygons is non-zero:
        for(unsigned int a = 0; a < mesh.NumberCell2D; ++a)
        {
            int n = mesh.Cell2DNum[a];
            double area_poligonoTOT = 0.0;


            //devo fare il for e l'if per dividere il poligono in triangoli
            /// triangolo n = 3 --> 0 linee i = 0:0 , 1 ∆
            /// quadrilatero n = 4 --> 1 linea i = 0,1 , 2∆
            /// pentagono n = 5 --> 2 linee , 3∆...


            if( n == 3){

                double areat = 0.0;
                areat = calculate_area(mesh.Cell2DVertices[a], mesh);

                if (areat < mesh.tolDefault){

                    cerr<<a<<" Zero area of single triangle : "<<scientific<<setprecision(6)<<areat<<endl;
                    little ++;
                    continue;

                }

                cout<<scientific<<setprecision(6)<<"\n "<<a<<" triangle area: "<<areat<<endl<<endl;

                area = area + areat;

            }
            else //polygons with n>3 edges
            {
                area_poligonoTOT=calculate_area(mesh.Cell2DVertices[a], mesh);
                if (area_poligonoTOT < mesh.tolDefault){
                    cerr<<"Zero polygonal area";
                    little ++;
                }
                else{
                    cout << a << " polygonal area, "<< n << " edges : "<< scientific<<setprecision(8)<< area_poligonoTOT<<endl<<endl;
                    area = area + area_poligonoTOT;

                }
            }


            cout<<"total area : "<<area<<endl;
            cout << "--------------------------------------"<<endl;


        }// primo for

        if (little == 0 && abs(area - 1) < mesh.tolDefault)
            cout <<"\n\t √ all areas of polygons are non-zero with this method\n"<<endl;
        else
            cerr << "\n\n\t † There are "<<little<<" triangles with very small area. (or something wrong with the method) \n\n"<<endl;

    }
    return true;

}



// ***************************************************************************
bool ImportCell0Ds(const string &filename,
                   PolygonalMesh& mesh)
{

    ifstream file;
    file.open(filename);

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    listLines.pop_front();

    mesh.NumberCell0D = listLines.size();
    cout<<"Cell 0D \nNumber of cells 0D = "<<mesh.NumberCell0D<<endl;

    if (mesh.NumberCell0D == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }

    mesh.Cell0DId.reserve(mesh.NumberCell0D);
    mesh.Cell0DCoordinates.reserve(mesh.NumberCell0D);
    int counter=0;
    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2d coord;
        char sep=';';


        converter >>  id >> sep >> marker >> sep >> coord(0) >> sep >> coord(1);

        mesh.Cell0DId.push_back(id);
        mesh.Cell0DCoordinates.push_back(coord);

        if( marker != 0)
        {

            ///first method

            if (mesh.Cell0DMarkers.find(marker) == mesh.Cell0DMarkers.end())
                mesh.Cell0DMarkers.insert({marker, {id}});
            else{
                mesh.Cell0DMarkers[marker].push_back(id);
            }
        }
        else{counter ++;}


    }

    cout << "there are "<<counter<< " 0D Cells inside the square, not on edges (marker = 0)."<<endl;
    file.close();
    return true;
}
// ***************************************************************************
bool ImportCell1Ds(const string &filename,
                   PolygonalMesh& mesh)
{

    ifstream file;
    file.open(filename);

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    listLines.pop_front();

    mesh.NumberCell1D = listLines.size();
    cout<<"\nCell 1D \nNumber of cells 1D = "<<mesh.NumberCell1D<<endl;

    if (mesh.NumberCell1D == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DId.reserve(mesh.NumberCell1D);
    mesh.Cell1DVertices.reserve(mesh.NumberCell1D);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2i vertices;
        char sep=';';


        converter >>  id >> sep >> marker >> sep >> vertices(0) >> sep >> vertices(1);

        mesh.Cell1DId.push_back(id);
        mesh.Cell1DVertices.push_back(vertices);

        if( marker != 0)
        {
            ///second method

            auto ret = mesh.Cell1DMarkers.insert({marker, {id}});
            if(!ret.second)
                (ret.first)->second.push_back(id);
        }
    }

    file.close();

    return true;
}
// ***************************************************************************
bool ImportCell2Ds(const string &filename,
                   PolygonalMesh& mesh)
{

    ifstream file;
    file.open(filename);

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    listLines.pop_front();

    mesh.NumberCell2D = listLines.size();

    if (mesh.NumberCell2D == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DId.reserve(mesh.NumberCell2D);
    mesh.Cell2DVertices.reserve(mesh.NumberCell2D);
    mesh.Cell2DEdges.reserve(mesh.NumberCell2D);
    mesh.Cell2DNum.reserve(mesh.NumberCell2D);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int mark;
        unsigned int numV;
        unsigned int numE;

        char sep=';';
        converter >>  id >> sep >> mark>>sep;
        converter >> numV >> sep;

        mesh.Cell2DNum.push_back(numV);

        vector<unsigned int> vertices(numV);


        for(unsigned int i = 0; i < numV; i++)
            converter >> vertices[i]>>sep;
        converter >> numE >> sep;

        if (numV != numE){
            cerr<<"in Cell2D numbers of vertices and edges don't correspond!!"<<endl;
            return false;
        }

        vector<unsigned int> edges(numE);

        for(unsigned int i = 0; i < numE; i++)
            converter >> edges[i]>>sep;

        mesh.Cell2DId.push_back(id);
        mesh.Cell2DVertices.push_back(vertices);
        mesh.Cell2DEdges.push_back(edges);
    }


    file.close();
    return true;
}

////////////////////////////////////////////////////
//////////////////funzioni per operazioni ricorrenti



double distance(const unsigned int origin, const unsigned int end, PolygonalMesh& mesh)
{
    Vector2d origin_coo = mesh.Cell0DCoordinates[origin];
    Vector2d end_coo = mesh.Cell0DCoordinates[end];

    //cout << "coord x Origin: " << origin_coo[0] << " coord y Origin: " << origin_coo[1]<<endl;
    //cout<<"coord x End: " << end_coo[0]<<" coord y End: " << end_coo[1] << endl;

    double lunghezza_lato = sqrt(pow(end_coo[0] - origin_coo[0], 2) + pow(end_coo[1] - origin_coo[1], 2));
    return lunghezza_lato;
}

double calculate_area( vector<unsigned int>& Cell2DVertices, PolygonalMesh& mesh){

    double areatr = 0.0;
    int n = size(Cell2DVertices);

    for ( int i = 1; i<n-1 ; i++){

        // A = 1/2 * ∑ |(x_i-x0) * (y_i+1-y0) - (x_i+1-x0) * (y_i-y0)|

        areatr += abs((mesh.Cell0DCoordinates[Cell2DVertices[i]][0] - mesh.Cell0DCoordinates[Cell2DVertices[0]][0]) *
                          (mesh.Cell0DCoordinates[Cell2DVertices[i+1]][1] - mesh.Cell0DCoordinates[Cell2DVertices[0]][1])
                      - (mesh.Cell0DCoordinates[Cell2DVertices[i+1]][0] - mesh.Cell0DCoordinates[Cell2DVertices[0]][0]) *
                            (mesh.Cell0DCoordinates[Cell2DVertices[i]][1] - mesh.Cell0DCoordinates[Cell2DVertices[0]][1]));
    }

    areatr = 0.5 * areatr;


    return areatr;

}

/////
/////formula alternativa che crea problemi e non tiene conto di certe aree perchè troppo piccole

double alternatove_method( PolygonalMesh mesh, vector<double>& edges_lenght, int a, int n){
    //formula di Herone

    int punto_fisso = mesh.Cell2DVertices[a][0];
    double s=0.0;
    vector<double> lunghezze;
    lunghezze[0] = edges_lenght[mesh.Cell2DEdges[a][0]];
    lunghezze[1] = edges_lenght[mesh.Cell2DEdges[a][1]];
    s = lunghezze[0]+lunghezze[1];

    for (int t = 2; t<n; t++){
        lunghezze[t] = distance(punto_fisso, mesh.Cell2DVertices[a][t], mesh);
        s += lunghezze[t];

    }
    s = s/2.0;
    double areatriangolino = s;
    double areapoli = 0.0;

    for (int m = 0; m<n-2; m++){
    for (int i = 0; i<n; i++){
        //Herone formula
        areatriangolino *= (abs(s - lunghezze[i]));
        if (areatriangolino < mesh.tolDefault){
            cerr<<a<<" Zero area of triangle in polygon: "<<scientific<<setprecision(6)<<areatriangolino<<endl;

        }else
            cout << "area triangolino : "<<areatriangolino<<endl;
    }
        areapoli += areatriangolino;
    }
    areapoli = sqrt( areapoli );


    if (areapoli < mesh.tolDefault){
        cerr<<"Zero polygonal area";
    }
    else{
        cout << a << " area poligono "<< n << " lati : "<< scientific<<setprecision(8)<< areapoli<<endl<<endl;
    }

    return areapoli;

}
////////////////////////////////////////////////////

}
