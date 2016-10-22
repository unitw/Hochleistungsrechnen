#include <stdio.h>

// Definieren Sie ein enum cardd

typedef enum {
    N = 1, S = 2, E = 4, W = 8
} cardd;


// Definieren Sie ein 3x3-Array namens map, das Werte vom Typ cardd enthält
cardd map[3][3];


//prüft ob die übergebene Richtung gültig ist oder nicht

int check_cardd(cardd dir) {
    return (dir >= N && dir <= 10 && dir != 3);
}

// Die Funktion set_dir soll an Position x, y den Wert dir in das Array map eintragen
// Überprüfen Sie x und y um mögliche Arrayüberläufe zu verhindern
// Überprüfen Sie außerdem dir auf Gültigkeit

void set_dir(int x, int y, cardd dir) {

    if ((x >= 0 && x < 3) && (y >= 0 && y < 3) && check_cardd(dir)) {
        map[x][y] = dir;

    }

}





// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben

void show_map(void) {


    for (int x = 0; x < 3; x++) {
        for (int y = 0; y < 3; y++) {

            switch (map[x][y]) {

                case N:
                    printf("|N|");
                    break;
                case E:
                    printf("|E|");
                    break;
                case S:
                    printf("|S|");
                    break;
                case W:
                    printf("|W|");
                    break;
                case N | W:
                    printf("|NW|");
                    break;
                case N | E:
                    printf("|NE|");
                    break;
                case S | W:
                    printf("|SW|");
                    break;
                case S | E:
                    printf("|SE|");
                    break;
                default:
                    printf("|0|");
                    break;


            }
            if (y == 2) {
                printf("\n");
            }
        }
    }
    printf("\n");


}

/*
int main(void) {
    // In dieser Funktion darf nichts verändert werden!
    set_dir(0, 1, N);
    set_dir(1, 0, W);
    set_dir(1, 4, W);
    set_dir(1, 2, E);
    set_dir(2, 1, S);

    show_map();

    set_dir(0, 0, N | W);
    set_dir(0, 2, N | E);
    set_dir(0, 2, N | S);
    set_dir(2, 0, S | W);
    set_dir(2, 2, S | E);
    set_dir(2, 2, E | W);

    show_map();

    return 0;
}
*/
