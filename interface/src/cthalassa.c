#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <cthalassa/cthalassa.h>

THALASSAPathStruct* THALASSAPathStruct_new(char *phys_path, char *earth_path, char *kernel_path) {
    // Allocate pointer
    THALASSAPathStruct *paths = malloc(sizeof(THALASSAPathStruct));

    // Physical constants path
    strcpy(paths->phys_path, phys_path);
    paths->phys_path_len = strlen(phys_path);

    // Earth model path
    strcpy(paths->earth_path, earth_path);
    paths->earth_path_len = strlen(earth_path);

    // SPICE kernel path
    strcpy(paths->kernel_path, kernel_path);
    paths->kernel_path_len = strlen(kernel_path);

    // Return pointer
    return paths;
}
