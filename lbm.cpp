#include <iostream>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <vector>
#include <math.h>

#include "imageClass/GrayScaleImage.h"
#include "imageClass/lodepng.h"

#define TIMESTEP 0.0015

//#define SPACEING 0.0005

using namespace std;

static bool scenario1 = false;
static double viscosity_l = 0.0;
static double timeToSimulate = 0.0;
static double acceleration = 0.0;
static double acceleration_l = 0.0;
static int resolution = 0;
static double timestep = 0.0;
static double spaceing = 0.0;
static double omega = 0.0;
static int numCellsX = 0;
static int numCellsY = 0;
static int steps = 0;


// template class for a 3D grid
template<typename T> class grid_lattice {

private:
    size_t lengthInX;
    size_t lengthInY;
    size_t lengthInF;
    std::vector<T> data;

public:
    grid_lattice() {
        lengthInX = 0;
        lengthInY = 0;
        lengthInF = 0;
    }
    grid_lattice(size_t xDim, size_t yDim, size_t f=9) {
        data = std::vector<T>((xDim+2) * (yDim+2) * f, (T) (0.0));
        lengthInX = xDim+2;
        lengthInY = yDim+2;
        lengthInF = f;
    }		// Standart-constructor
    grid_lattice(size_t xDim, size_t yDim, T value, size_t f=9) {
        data = std::vector<T>((xDim+2) * (yDim+2) * f, value);
        lengthInX = xDim+2;
        lengthInY = yDim+2;
        lengthInF = f;
    }	// initalisation-constructor
    virtual ~grid_lattice() {
    }		// Destructor

    size_t lengthX() {
        return lengthInX;
    }	// returns number of elements in x direction
    size_t lengthY() {
        return lengthInY;
    }	// returns number of elements in y direction
    size_t lengthF() {
        return lengthInF;
    }	// returns number of elements in f direction

    T& operator()(size_t i, size_t j, size_t f=0) {
        assert(i < lengthInX);
        assert(j < lengthInY);
        assert(f < lengthInF);
        return data[j * ((int)lengthInX*(int)lengthInF) + i*(int)lengthInF +(int)f];
    }

    void getCopy( grid_lattice &source ){
        for( int i=0; i<(int)lengthInX; i++){
            for( int j=0; j<(int)lengthInY; j++){
                for( int f=0; f<(int)lengthInF; f++){
                    data[j * ((int)lengthInX*(int)lengthInF) + i*(int)lengthInF +f] = source( i,j,f);
                }
            }
        }
    }

};

static grid_lattice<double> grid;
static grid_lattice<double> gridCopy;
static grid_lattice<int> isBoundary;
static grid_lattice<double> equilibrium;

void streamStep(){
    gridCopy.getCopy( grid );

    int lengthX = (int)grid.lengthX();
    int lengthY = (int)grid.lengthY();

    // inner grid points without cylinder
    for( int i=1; i<lengthX-1; i++ ){
        for( int j=1; j<lengthY-1; j++ ){
            if( isBoundary(i,j) == 0){
                grid(i+1, j,   1) = gridCopy( i,j,1);
                grid(i+1, j+1, 2) = gridCopy( i,j,2);
                grid(i,   j+1, 3) = gridCopy( i,j,3);
                grid(i-1, j+1, 4) = gridCopy( i,j,4);
                grid(i-1, j,   5) = gridCopy( i,j,5);
                grid(i-1, j-1, 6) = gridCopy( i,j,6);
                grid(i,   j-1, 7) = gridCopy( i,j,7);
                grid(i+1, j-1, 8) = gridCopy( i,j,8);
            }
        }
    }
    
    //Es stehen noch keine RB im gridCopy!!!! Deshalb wird 0 kopiert

    // periodic RB
    for (int i=1; i<lengthY-1; i++){
        // linke RB
        grid(lengthX-2, i, 4) = grid( 0,i,4);
        grid(lengthX-2, i, 5) = grid( 0,i,5);
        grid(lengthX-2, i, 6) = grid( 0,i,6);
        // rechte RB
        grid(1, i, 1) = grid( lengthX-1,i,1);
        grid(1, i, 2) = grid( lengthX-1,i,2);
        grid(1, i, 8) = grid( lengthX-1,i,8);
    }
    
    //no slip RB
    for (int i=1; i<lengthX-1; i++){
        // obere RB
        grid(i-1, lengthY-2, 6) = grid( i,lengthY-1,2);
        grid(i, lengthY-2, 7) = grid( i,lengthY-1,3);
        grid(i+1, lengthY-2, 8) = grid( i,lengthY-1,4);
        // untere RB
        grid(i-1, 1, 4) = grid( i,0,8);
        grid(i, 1, 3) = grid( i,0,7);
        grid(i+1, 1, 2) = grid( i,0,6);
    }
    grid(1, 1, 2) = grid( 0,0, 6);
    grid(1, lengthY-2, 8) = grid( 0,lengthY-1, 4);
    grid(lengthX-2, 1, 4) = grid( lengthX-1,0, 8);
    grid(lengthX-2,lengthY-2, 6) = grid( lengthX-1,lengthY-1, 2);

    // Cylinder cells
    for( int i=2; i<lengthX-2; i++ ){
        for( int j=2; j<lengthY-2; j++ ){
            if( isBoundary(i,j)==1 ){
                grid(i-1, j,   5) = grid( i,j,1);
                grid(i, j+1,   3) = grid( i,j,7);
                grid(i, j-1,   7) = grid( i,j,3);
                grid(i+1, j,   1) = grid( i,j,5);
                grid(i-1, j-1, 6) = grid( i,j,2);
                grid(i+1, j-1, 8) = grid( i,j,4);
                grid(i+1, j+1, 2) = grid( i,j,6);
                grid(i-1, j+1, 4) = grid( i,j,8);
            }
        }
    }
}

void initBoundBoolean(){

    //double magic = sqrt(0.5)*SPACEING;
    //magic = 0.5*SPACEING;

    int lengthX = (int)grid.lengthX();
    int lengthY = (int)grid.lengthY();

    for(int i=0; i<lengthX; ++i){
        isBoundary(i,0) = 1;
        isBoundary(i,lengthY-1) = 1;
    }
    for(int i=0; i<lengthY; ++i){
        isBoundary(0,i) = 1;
        isBoundary(lengthX-1,i) = 1;
    }
    for(int i=1; i<lengthY-1; ++i){
        for(int j=1; j<lengthX-1; ++j){
            if(    (sqrt((i*spaceing-0.008)*(i*spaceing-0.008) + (j*spaceing-0.02)*(j*spaceing-0.02))) <= 0.0025 + 0.5*spaceing  && (sqrt((i*spaceing-0.008)*(i*spaceing-0.008) + (j*spaceing-0.02)*(j*spaceing-0.02))) >= 0.0025 - 0.5*spaceing ){
                isBoundary(j,i) = 1;
            }else if((sqrt((i*spaceing-0.008)*(i*spaceing-0.008) + (j*spaceing-0.02)*(j*spaceing-0.02))) <= 0.0025 -0.5*spaceing){
                isBoundary(j,i) = -1;
            }
            else{
                //isBoundary(j,i) = 0;
            }
        }
    }
    for(int i=1; i<lengthY-1; ++i){
        for(int j=1; j<lengthX-1; ++j){
            if( isBoundary(j,i) == -1 ){
                    if( isBoundary(j-1,i-1) == 0 || isBoundary(j-1,i) == 0 ||isBoundary(j-1,i+1) == 0 || isBoundary(j,i-1) == 0 || isBoundary(j,i) == 0 || isBoundary(j,i+1) == 0 || isBoundary(j+1,i-1) == 0 || isBoundary(j+1,i) == 0 || isBoundary(j+1,i+1) == 0 ){
                        isBoundary(j,i) = 1;
                    }
            }
        }
    }
    
}

void initLaticeGrid(){
    int lengthX = (int)grid.lengthX();
    int lengthY = (int)grid.lengthY();

    for(int i=1; i<lengthX-1; ++i){
        for(int j=1; j<lengthY-1; ++j){
            if(isBoundary(i,j) == 0){
                    grid(i,j,0) = 4.0/9.0;
                    grid(i,j,1) = 1.0/9.0;
                    grid(i,j,3) = 1.0/9.0;
                    grid(i,j,5) = 1.0/9.0;
                    grid(i,j,7) = 1.0/9.0;
                    grid(i,j,2) = 1.0/36.0;
                    grid(i,j,4) = 1.0/36.0;
                    grid(i,j,6) = 1.0/36.0;
                    grid(i,j,8) = 1.0/36.0;
            }
        }
    }
}

void collideStep(){

    int lengthX = (int)grid.lengthX();
    int lengthY = (int)grid.lengthY();
    double density = 0.0;
    double velosity[2];
    double skalar = 0.0;
    double skalar_u = 0.0;
    double vier_neun = 4.0/9.0;
    double ein_neun = 1.0/9.0;
    double ein_sechs = 1.0/36.0;

    for( int i=1; i<lengthY-1; i++ ){
        for( int j=1; j<lengthX-1; j++ ){
            if(isBoundary(j,i) == 0){

                for(int k=0; k<9; k++){
                    density += grid(j,i,k);
                }

                velosity[0] = grid(j,i,1)+ grid(j,i,2) - grid(j,i,4) - grid(j,i,5) - grid(j,i,6) + grid(j,i,8);
                velosity[1] = grid(j,i,2)+ grid(j,i,3) + grid(j,i,4) - grid(j,i,6) - grid(j,i,7) - grid(j,i,8);
                velosity[0] = velosity[0] / density;
                velosity[1] = velosity[1] / density;

                // Compute equilibrium distribution
                skalar_u = (velosity[0] * velosity[0]) + (velosity[1] * velosity[1]);

                skalar = 0;
                equilibrium(j,i,0) = vier_neun * density * (1.0 + 3.0*skalar + 4.5*skalar*skalar - 1.5* skalar_u);

                skalar = 1*velosity[0];
                equilibrium(j,i,1) = ein_neun * density * (1.0 + 3.0*skalar + 4.5*skalar*skalar - 1.5* skalar_u);

                skalar = 1*velosity[0] + 1*velosity[1];
                equilibrium(j,i,2) = ein_sechs * density * (1.0 + 3.0*skalar + 4.5*skalar*skalar - 1.5* skalar_u);

                skalar = 1*velosity[1];
                equilibrium(j,i,3) = ein_neun * density * (1.0 + 3.0*skalar + 4.5*skalar*skalar - 1.5* skalar_u);

                skalar = -1*velosity[0] + 1*velosity[1];
                equilibrium(j,i,4) = ein_sechs * density * (1.0 + 3.0*skalar + 4.5*skalar*skalar - 1.5* skalar_u);

                skalar = -1*velosity[0];
                equilibrium(j,i,5) = ein_neun * density * (1.0 + 3.0*skalar + 4.5*skalar*skalar - 1.5* skalar_u);

                skalar = -1*velosity[0] + -1*velosity[1];
                equilibrium(j,i,6) = ein_sechs * density * (1.0 + 3.0*skalar + 4.5*skalar*skalar - 1.5* skalar_u);

                skalar = -1*velosity[1];
                equilibrium(j,i,7) = ein_neun * density * (1.0 + 3.0*skalar + 4.5*skalar*skalar - 1.5* skalar_u);

                skalar = 1*velosity[0] + -1*velosity[1];
                equilibrium(j,i,8) = ein_sechs * density * (1.0 + 3.0*skalar + 4.5*skalar*skalar - 1.5* skalar_u);

                // Update rule for grid
                grid(j,i,0) = grid(j,i,0) - omega * (grid(j,i,0) - equilibrium(j,i,0)) + 3 * vier_neun * density * 0 * acceleration_l;
                grid(j,i,1) = grid(j,i,1) - omega * (grid(j,i,1) - equilibrium(j,i,1)) + 3 * ein_neun * density * 1 * acceleration_l;
                grid(j,i,2) = grid(j,i,2) - omega * (grid(j,i,2) - equilibrium(j,i,2)) + 3 * ein_sechs * density * 1 * acceleration_l;
                grid(j,i,3) = grid(j,i,3) - omega * (grid(j,i,3) - equilibrium(j,i,3)) + 3 * ein_neun * density * 0 * acceleration_l;
                grid(j,i,4) = grid(j,i,4) - omega * (grid(j,i,4) - equilibrium(j,i,4)) + 3 * ein_sechs * density * -1 * acceleration_l;
                grid(j,i,5) = grid(j,i,5) - omega * (grid(j,i,5) - equilibrium(j,i,5)) + 3 * ein_neun * density * -1 * acceleration_l;
                grid(j,i,6) = grid(j,i,6) - omega * (grid(j,i,6) - equilibrium(j,i,6)) + 3 * ein_sechs * density * -1 * acceleration_l;
                grid(j,i,7) = grid(j,i,7) - omega * (grid(j,i,7) - equilibrium(j,i,7)) + 3 * ein_neun * density * 0 * acceleration_l;
                grid(j,i,8) = grid(j,i,8) - omega * (grid(j,i,8) - equilibrium(j,i,8)) + 3 * ein_sechs * density * 1 * acceleration_l;

                //cout << "density: " << density << endl;
                density = 0.0;
            }
        }
    }
}

void testStream(){
    int lengthX = (int)grid.lengthX();
    int lengthY = (int)grid.lengthY();

    for(int i=1; i<lengthX-1; ++i){
        for(int j=1; j<lengthY-1; ++j){
            if(isBoundary(i,j) == 0){
                cout << "cell: " << i << "," << j << endl;
                for(int f=0; f<9; ++f){
                    if(f==0){
                        if(grid(i,j,f) == 4.0/9.0){
                            cout << f << " Richtig " << grid(i,j,f) << endl;
                        }else{
                            cout << f << " Falsch " << grid(i,j,f) << endl;
                        }
                    }

                    if(f==1 || f==3 || f==5 || f==7){
                        if(grid(i,j,f) == 1.0/9.0){
                            cout << f << " Richtig " << grid(i,j,f) << endl;
                        }else{
                            cout << f << " Falsch " << grid(i,j,f) << endl;
                        }
                    }

                    if(f==2 || f==4 || f==6 || f==8){
                        if(grid(i,j,f) == 1.0/36.0){
                            cout << f << " Richtig " << grid(i,j,f) << endl;
                        }else{
                            cout << f << " Falsch " << grid(i,j,f) << endl;
                        }
                    }

                }
                cout << endl;
            }
        }
    }
}

void printDomain(){
    int lengthX = (int)grid.lengthX();
    int lengthY = (int)grid.lengthY();

        for(int i=lengthY-1; i>=0; --i){
            for(int j=0; j<lengthX; ++j){
                if(isBoundary(j,i) == 0)
                    cout << "c" << " ";
                else if(isBoundary(j,i) == 1)
                    cout << "b" << " ";
                else if(isBoundary(j,i) == -1)
                    cout << "z" << " ";
                else
                    cout << "!" << " ";
            }
            cout << endl;
        }
}

void ausgabeBild(){



    double velosity [2];
    double density= 0.0;
    double speed = 0.0;
    double max_speed = 0.0;

    int lengthX = (int)grid.lengthX();
    int lengthY = (int)grid.lengthY();

    GrayScaleImage image = GrayScaleImage(numCellsX, numCellsY);

    for(int j = 1; j<lengthY-1; ++j){
        for(int i = 1; i<lengthX-1; ++i){
            if(isBoundary(i,j) == 0){
                for(int k=0; k<9; k++){
                    density += grid(i,j,k);
                }
                //cout << "density bild: " << density << endl;
                velosity[0] = grid(i,j,1)+ grid(i,j,2) - grid(i,j,4) - grid(i,j,5) - grid(i,j,6) + grid(i,j,8);
                velosity[1] = grid(i,j,2)+ grid(i,j,3) + grid(i,j,4) - grid(i,j,6) - grid(i,j,7) - grid(i,j,8);
                velosity[0] = velosity[0] / density;
                velosity[1] = velosity[1] / density;
                speed = sqrt(velosity[0]*velosity[0]+velosity[1]*velosity[1]);
                //cout << "speed bild: " << speed << endl;
                if(speed > max_speed){
                    max_speed = speed;
                }
                density = 0.0;
            }
        }
    }

    cout << "max speed bild: " << max_speed << endl;

    for(int j = 1; j<lengthY-1; ++j){
        for(int i = 1; i<lengthX-1; ++i){
            if(isBoundary(i,j) == 0){
                for(int k=0; k<9; k++){
                    density += grid(i,j,k);
                }
                velosity[0] = grid(i,j,1)+ grid(i,j,2) - grid(i,j,4) - grid(i,j,5) - grid(i,j,6) + grid(i,j,8);
                velosity[1] = grid(i,j,2)+ grid(i,j,3) + grid(i,j,4) - grid(i,j,6) - grid(i,j,7) - grid(i,j,8);
                velosity[0] = velosity[0] / density;
                velosity[1] = velosity[1] / density;
                speed = sqrt(velosity[0]*velosity[0]+velosity[1]*velosity[1]);

                image.setElement(i-1, j-1, (speed/max_speed));

                density=0.0;
            }
        }
     }
    if(scenario1){
            image.save("scenario1.png");
    }else{
        image.save("scenario2.png");
    }


}


int main( int args, char** argv ){

    if( args != 2 ){
        cout << "USAGE: lbm scenario" << endl;
        exit( EXIT_SUCCESS );
    }

    if( strcmp( argv[1], "scenario1") == 0 ){
        timestep = TIMESTEP;
        spaceing = 0.005/30;
        scenario1 = true;
        viscosity_l = 1e-6 * (timestep/(spaceing*spaceing));
        timeToSimulate = 3.0;
        acceleration = 0.01;
        acceleration_l = (acceleration*timestep*timestep)/spaceing;
        resolution = 30;
        omega = 1.0 / (3.0*viscosity_l + 0.5 );
        numCellsX = (int)(0.06 / spaceing);
        numCellsY = (int)(0.02 / spaceing);
        grid = grid_lattice<double>( numCellsX, numCellsY );
        gridCopy = grid_lattice<double>( numCellsX, numCellsY );
        isBoundary = grid_lattice<int>( numCellsX, numCellsY, 0, 1 );
        equilibrium = grid_lattice<double>( numCellsX, numCellsY );
        steps = timeToSimulate/timestep;
    }else if( strcmp( argv[1], "scenario2") == 0){
        timestep = TIMESTEP;
        spaceing = 0.005/60;
        scenario1 = false;
        viscosity_l = 1e-6 * (timestep/(spaceing*spaceing));
        timeToSimulate = 5.0;
        acceleration = 0.016;
        acceleration_l = (acceleration*timestep*timestep)/spaceing;
        resolution = 60;
        omega = 1.0 / (3.0*viscosity_l + 0.5 );
        numCellsX = 0.06 / spaceing;
        numCellsY = 0.02 / spaceing;
        grid = grid_lattice<double>( numCellsX, numCellsY );
        gridCopy = grid_lattice<double>( numCellsX, numCellsY );
        isBoundary = grid_lattice<int>( numCellsX, numCellsY, 0, 1 );
        equilibrium = grid_lattice<double>( numCellsX, numCellsY );
        steps = timeToSimulate/timestep;
    }else{
        cout << "scenario has to be 1 or 2" << endl;
        exit( EXIT_SUCCESS );
    }
    cout << "X " << numCellsX << endl;
    cout << "Y " << numCellsY << endl;
    cout << "timestep " << timestep << endl;
    cout << "spaceing " << spaceing << endl;
    cout << "accel_l " << acceleration_l << endl;
    cout << "omega " << omega << endl;
    cout << "viscosity_l " << viscosity_l << endl;



    cout << "Start Inits" << endl;
    initLaticeGrid();
    initBoundBoolean();

    /*
    //Test Stream Step
    streamStep();
    testStream();
    */
    /*
    //Test Collide Step
    collideStep();
    testStream();
    */


    cout << "Start Berechnung" << endl;
    for(int i=0; i<steps; ++i){
        streamStep();
        collideStep();
    }


    cout << "Start Ausgabe Bild" << endl;
    ausgabeBild();


    exit( EXIT_SUCCESS );
}

















