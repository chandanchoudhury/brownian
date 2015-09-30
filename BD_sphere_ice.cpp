
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <unistd.h>

#include "subroutines.h"

using namespace std;

int main(int argc, const char **argv) {
    double **coords = new double*[n_atoms];
    int i = 0, Site[10] = {0}, count = 0; 
    // bool overlap = false;
    int overlap = 1;

    // Randomly generate the coordinates of the atoms in the box
    srand(time(NULL));
    for (int i = 0; i < n_atoms; i = i + 1) {
        coords[i] = new double[3];

        // Note "rand(0,x)" would generate a random number
        // between 0 and $x
        coords[i][0] = rand(0, box_size[0]);
        coords[i][1] = rand(0, box_size[1]);
        coords[i][2] = rand(0, box_size[2]);
    }

    // The total number of accepted moves
    int naccept = 0;

    // The total number of rejected moves
    int nreject = 0;

    // Print the initial PDB file
    print_pdbeq(coords, n_atoms, 0, nsites, sigmaI);
    printf("\n");
    printf("Equilibration Steps . . . . . \n");

    for (int move = 1; move <= eq_moves; move = move + 1) {
        for(int atom = 0; atom < n_atoms; atom++) {
        // int atom = int(rand(0, n_atoms));
        // save the old coordinates
        const double old_coords[3] = { coords[atom][0], coords[atom][1],
                                       coords[atom][2] };

        // Make the move - translate by a delta in each dimension
        const double delta_x = rand_normal(0.0, 1.0) * sqrt(Ddt);
        const double delta_y = rand_normal(0.0, 1.0) * sqrt(Ddt);
        const double delta_z = rand_normal(0.0, 1.0) * sqrt(Ddt);
/*      const double d       = rand(-max_translate, max_translate);
        const double theta   = rand(0, M_PI);
        const double phi     = rand(0, 2*M_PI);
        const double delta_x = d*sin(theta)*cos(phi);
        const double delta_y = d*sin(theta)*sin(phi);
        const double delta_z = d*cos(theta);
*/
        coords[atom][0] = coords[atom][0] + delta_x;
        coords[atom][1] = coords[atom][1] + delta_y;
        coords[atom][2] = coords[atom][2] + delta_z;

        // wrap the coordinates back into the box
        coords[atom][0] = wrap_into_box(coords[atom][0], box_size[0]);
        coords[atom][1] = wrap_into_box(coords[atom][1], box_size[1]);
        coords[atom][2] = wrap_into_box(coords[atom][2], box_size[2]);

        const double new_energy = calculate_energy(coords, n_atoms, atom, box_size, sigmaP, epsilon);
        // Automatically accept if the energy goes down
        if (new_energy == 0) {
            naccept = naccept + 1;
        }
        else {
            nreject = nreject + 1;
            coords[atom][0] = old_coords[0];
            coords[atom][1] = old_coords[1];
            coords[atom][2] = old_coords[2];
        }
        }

        if (move % 10 == 0) {
            printf("%8d:  %8d  %6d  %8.3f\n", move, naccept, nreject, sigmaI);
        }
        if (move % 10 == 0) {
            print_pdbeq(coords, n_atoms, move, nsites, sigmaI);
        }
    }

    printf("Post Equilibration Steps . . . . . \n");
    print_pdb(coords, n_atoms, 0, nsites, sigmaI);

    for (int move = 1; move <= num_moves; move = move + 1) {
        // Pick a random Particle atom. We do not pick the Ice particles.
        // int atom = int(rand(0, n_atoms-nsites));
        for(int atom = 0; atom < n_atoms-nsites; atom++) {
            overlap = 1; // False
            int Ioverlap = 0;
            // Save the old coordinates
            const double old_coords[3] = { coords[atom][0], coords[atom][1],
                                           coords[atom][2] };
            // Make the move - translate by a delta in each dimension
            const double delta_x = rand_normal(0.0, 1.0) * sqrt(Ddt);
            const double delta_y = rand_normal(0.0, 1.0) * sqrt(Ddt);
            const double delta_z = rand_normal(0.0, 1.0) * sqrt(Ddt);

/*          const double delta_x = rand_normal(-max_translate, max_translate, 0.0, 1.0);
            const double delta_y = rand_normal(-max_translate, max_translate, 0.0, 1.0);
            const double delta_z = rand_normal(-max_translate, max_translate, 0.0, 1.0);
*/
            coords[atom][0] = coords[atom][0] + delta_x;
            coords[atom][1] = coords[atom][1] + delta_y;
            coords[atom][2] = coords[atom][2] + delta_z;

            // Wrap the coordinates back into the box
            coords[atom][0] = wrap_into_box(coords[atom][0], box_size[0]);
            coords[atom][1] = wrap_into_box(coords[atom][1], box_size[1]);
            coords[atom][2] = wrap_into_box(coords[atom][2], box_size[2]);

            // Particle - Particle or Particle-Ice distance
            calculate_overlap(coords, atom, n_atoms, nsites, sigmaP, sigmaI, &overlap, &Ioverlap);

            if(overlap == 0) { // Overlap True
                if(Ioverlap >= n_atoms-nsites) {
                    // Restore the old coords.
                    coords[atom][0] = old_coords[0];
                    coords[atom][1] = old_coords[1];
                    coords[atom][2] = old_coords[2];
                    // Change the direction of the move.
                    coords[atom][0] = coords[atom][0] - delta_x;
                    coords[atom][1] = coords[atom][1] - delta_y;
                    coords[atom][2] = coords[atom][2] - delta_z;

                    // Wrap the coordinates back into the box
                    coords[atom][0] = wrap_into_box(coords[atom][0], box_size[0]);
                    coords[atom][1] = wrap_into_box(coords[atom][1], box_size[1]);
                    coords[atom][2] = wrap_into_box(coords[atom][2], box_size[2]);

                    overlap = 1; Ioverlap = 0;

                    calculate_overlap(coords, atom, n_atoms, nsites, sigmaP, sigmaI, &overlap, &Ioverlap);

                    if(overlap == 0) { // Overlap True for Ice - particle
                        nreject = nreject + 1;
                        coords[atom][0] = old_coords[0];
                        coords[atom][1] = old_coords[1];
                        coords[atom][2] = old_coords[2];
                    }
                    else {
                    naccept = naccept + 1;
                    }
                }
                else { // Particle - particle overlap
                    nreject = nreject + 1;
                    coords[atom][0] = old_coords[0];
                    coords[atom][1] = old_coords[1];
                    coords[atom][2] = old_coords[2];
                }
            }
            else
                naccept = naccept + 1;
        }
        // Print the energy every 1000 moves
        if (move % 1000 == 0) {
            printf("%8d:  %12d  %12d  %8.3f\n", move, naccept, nreject, sigmaI);
        }
        // print the coordinates every 10000 moves
        if (move % 1000 == 0) {
            print_pdb(coords, n_atoms, move, nsites, sigmaI);
        }
        sigmaI += IceGrowth; //Radially increasing the Ice volume
    }
    return 0;
}
