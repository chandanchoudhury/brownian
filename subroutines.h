#include"constants_MP.h"

// function to return a random number between 'start' to 'end'
double rand(const double start, const double end) {
    return (end-start) * (double(rand()) / RAND_MAX) + start;
}

double rand_normal(const double mean, const double stddev) { //Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached) {
        double x, y, r;
        do {
            x = ((double) rand() / (RAND_MAX)) * 2 - 1;
            y = ((double) rand() / (RAND_MAX)) * 2 - 1;
            // x = (end-start) * (double(rand()) / RAND_MAX) + start;
            // y = (end-start) * (double(rand()) / RAND_MAX) + start;

            r = x*x + y*y;            
        } while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
        
    }
    else {
        n2_cached = 0; return n2*stddev + mean;
    }
}

// Subroutine to apply periodic boundaries
double make_periodic(double x, const double box) {
    while (x < -0.5*box) {
        x = x + box;
    }

    while (x > 0.5*box) {
        x = x - box;
    }

    return x;
}

// Subroutine to wrap the coordinates into a box
double wrap_into_box(double x, double box) {
    while (x > box) {
        x = x - box;
    }

    while (x < 0) {
        x = x + box;
    }

    return x;
}

// Subroutine to print a PDB of the coordinates: Equilibrium Stage
void print_pdbeq(double **coords, const int n_atoms, const int move, const int nsites, const double sigmaI) {
    char filename[128];
    int i = 0;
    
    snprintf(filename, 128, "Eq%000008d.pdb", move);

    FILE *f = fopen(filename, "w");

    fprintf(f, "TITLE   move = %10d, sigmaI = %8.4f\n",move, sigmaI);
    fprintf(f, "CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00\n", 
               box_size[0], box_size[1], box_size[2]);

    for (i = 0; i < n_atoms-nsites; i = i + 1) {
        fprintf(f, "ATOM  %5d  Kr   Kr     1    %8.3f%8.3f%8.3f  1.00  0.00          Kr\n",
                       i+1, coords[i][0], coords[i][1], coords[i][2]); 
        fprintf(f, "TER\n");
    }

    for (int j = n_atoms-nsites; j < n_atoms; j = j + 1) {  
        fprintf(f, "ATOM  %5d  Pb   Pb     1    %8.3f%8.3f%8.3f  1.00  0.00          Pb\n",\
        i+1, coords[j][0], coords[j][1], coords[j][2]);

        fprintf(f, "TER\n"); i++;
    }

    fclose(f);
}
// Subroutine to print a PDB of the coordinates
void print_pdb(double **coords, const int n_atoms, const int move, const int nsites, const double sigmaI) {
    char filename[128];
    int i = 0;
    
    snprintf(filename, 128, "output%000008d.pdb", move);

    FILE *f = fopen(filename, "w");

    fprintf(f, "TITLE   move = %10d, sigmaI = %8.4f\n",move, sigmaI);
    fprintf(f, "CRYST1 %8.3f %8.3f %8.3f  90.00  90.00  90.00\n", 
               box_size[0], box_size[1], box_size[2]);

    for (i = 0; i < n_atoms-nsites; i = i + 1) {
        fprintf(f, "ATOM  %5d  Kr   Kr     1    %8.3f%8.3f%8.3f  1.00  0.00          Kr\n",
                       i+1, coords[i][0], coords[i][1], coords[i][2]); 
        fprintf(f, "TER\n");
    }

    for (int j = n_atoms-nsites; j < n_atoms; j = j + 1) {  
        fprintf(f, "ATOM  %5d  Pb   Pb     1    %8.3f%8.3f%8.3f  1.00  0.00          Pb\n",\
        i+1, coords[j][0], coords[j][1], coords[j][2]);

        fprintf(f, "TER\n"); i++;
    }

    fclose(f);
}

// Subroutine that calculates the energies of the atoms
double calculate_energy(double **coords, const int n_atoms, const int atom, const double *box_size,
                        const double sigmaP, const double epsilon) {
    // Loop over all pairs of atoms and calculate
    // the LJ energy
    double total_energy = 0;

    for (int i = 0; i < n_atoms; i = i + 1) {
        if(i != atom) {
            double delta_x = coords[atom][0] - coords[i][0];
            double delta_y = coords[atom][1] - coords[i][1];
            double delta_z = coords[atom][2] - coords[i][2];

            // Apply periodic boundaries
            delta_x = make_periodic(delta_x, box_size[0]);
            delta_y = make_periodic(delta_y, box_size[1]);
            delta_z = make_periodic(delta_z, box_size[2]);

            const double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                              (delta_z*delta_z);

            const double r = sqrt(r2);
            if(r <= (sigmaP + sigmaP + rand(0, 0.1)))
                total_energy += 1.0;
            else
                total_energy += 0.0;
            /*
            E_LJ = 4*epsilon[ (sigma/r)^12 - (sigma/r)^6 ]
            const double sig2_over_r2 = (sigmaP*sigmaP) / r2;
            const double sig6_over_r6 = sig2_over_r2*sig2_over_r2*sig2_over_r2;
            const double sig12_over_r12 = sig6_over_r6 * sig6_over_r6;

            const double e_lj = 4.0 * epsilon * ( sig12_over_r12 - sig6_over_r6 );

            total_energy = total_energy + e_lj;
            */
        }
    }
    // return the total energy of the atoms
    return total_energy;
}

// Subroutine that predicts the overlap of the particles
//P-P
//P-I
bool calculate_overlap(double **coords, const int atom, const int n_atoms, const int nsites, const double sigmaP, const double sigmaI, int *overlap, int *Ioverlap)
{
   *Ioverlap = 0; 
    // Loop over all atoms and calculate
    for(int i = 0; i < n_atoms; i = i + 1) {
        if(i != atom) {
        double delta_x = coords[i][0] - coords[atom][0];
        double delta_y = coords[i][1] - coords[atom][1];
        double delta_z = coords[i][2] - coords[atom][2];

        // Apply periodic boundaries
        delta_x = make_periodic(delta_x, box_size[0]);
        delta_y = make_periodic(delta_y, box_size[1]);
        delta_z = make_periodic(delta_z, box_size[2]);

        const double r2 = (delta_x*delta_x) + (delta_y*delta_y) +
                           (delta_z*delta_z);
        const double r  = sqrt(r2);

//      atom is never a part of Ice.
        if(i < (n_atoms-nsites)) // i is a particle
        {
            if(r <= (sigmaP + sigmaP)) {
                *overlap = 0; 
                return true; //immediately transfers the control back to the calling function.
            }
            else {
                *overlap = 1;
            }
        }
        else // if i is an Ice
        {
            if(r <= (sigmaP + sigmaI)) {
                *overlap = 0; *Ioverlap = i;
                return true;
            }
            else
                *overlap = 1;
        }
        }
    }
    return 0;
}
