#include "Checker.h"
# include <iostream>
# include <fstream>
void Checker::run(){
  
    int k_first_event, cell_first_event;
    double displacement_max =
    min(min(Param::box[0], Param::box[1])/2, min(Param::cell_size[0], Param::cell_size[1])) - 2*Param::sigma;
    ofstream oszero ("NumberUnvisitedConstraints.dat");
    double distance_to_go = 50000000*Param::lambda_0 *Param::factor; //standard run length

    time_start = clock();
  
    while (true) {   //loop over collisions
        double displacement_first_event = distance_to_go;
        expl_cell( displacement_first_event, k_first_event, cell_first_event);

        displacement_first_event = max(displacement_first_event, .0);
        double new_position =
                cell[cell_cur][k_cur][0] + min(min(displacement_first_event, distance_to_go), displacement_max);

        refresh_cell( new_position );

        if (displacement_max < min(displacement_first_event,distance_to_go)) {
            distance_to_go -= displacement_max;
            continue; //move limited by cell size, the move continues with the same disk and the same direction
        }
        else if (displacement_first_event < distance_to_go) {
            distance_to_go = distance_to_go - displacement_first_event;
            k_cur = k_first_event;
            cell_cur = cell_first_event;
            number_liftings++;
            if(number_liftings % Param::outcount == 0) countzero(number_liftings, oszero);
            continue; //a regular lifting move
        }
        else {
            break; //end of the chain
        }
    }
    time_end = clock();
}
