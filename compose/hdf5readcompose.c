/*
 *  This file contains some sample routines to read data from the HDF5 file.
 */

#define H5_USE_16_API
#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RANK1   1
#define file_name "eoscompose.h5"


void hdf5_read_dimensions_(int* m_nb, int*n_temp, int* o_y_q, int * n_qty_rd,int * n_add_rd,int * n_p_rd,int * n_q_rd,int * n_m_rd,int * n_err_rd,int * i_tab_rd)
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    herr_t      status;

    /* open existing file */    
    file = H5Fopen (file_name, H5F_ACC_RDONLY,H5P_DEFAULT);

    /* points for n_b */

    /*
     * Open an existing data set
     */
    dataset = H5Dopen2(file, "pointsnb", H5P_DEFAULT);
    /*
     * Read the data from the dataset using default transfer properties.
     */
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);
    
/*      * Close/release resources.
     */
    H5Dclose(dataset);

    /* points for T */
    dataset = H5Dopen2(file, "pointst", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_temp);
    H5Dclose(dataset);

    /* points for y_q */

    dataset = H5Dopen2(file, "pointsyq", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, o_y_q);
    H5Dclose(dataset);


    /* read number of different thermo quantities */

    dataset = H5Dopen2(file, "pointsqty", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_qty_rd);
    H5Dclose(dataset);

    /* read itab */

    dataset = H5Dopen2(file, "itab", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, i_tab_rd);
    H5Dclose(dataset);

    /* read number of different additional thermo quantities */

    dataset = H5Dopen2(file, "pointsadd", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_add_rd);
    H5Dclose(dataset);

    /* read number of pairs */

    dataset = H5Dopen2(file, "pointspairs", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_p_rd);
    H5Dclose(dataset);


    /* read number of quadruples */

    dataset = H5Dopen2(file, "pointsav", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_q_rd);
    H5Dclose(dataset);

    /* read number of microscopic quantities */

    dataset = H5Dopen2(file, "pointsmicro", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_m_rd);
    H5Dclose(dataset);

    /* read number of error quantities */

    dataset = H5Dopen2(file, "pointserr", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_err_rd);
    H5Dclose(dataset);



    /* finally closing the file */
    H5Fclose(file);

}


void hdf5_read_(double * nb, double * t, double * y_q, double * quantity, int * index_quantity,int * openfile, char * qty_name, char* index_name)
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    herr_t      status;



    /* open existing file */    
    file = H5Fopen (file_name, H5F_ACC_RDONLY,H5P_DEFAULT);

    if(*openfile == 1) { 

    /* read points for density */

    dataset = H5Dopen2(file, "nb", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nb);
    H5Dclose(dataset);

    /* read points for temperature */

    dataset = H5Dopen2(file, "t", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
    H5Dclose(dataset);

    /* read points for charge fraction */

    dataset = H5Dopen2(file, "yq", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y_q);
    H5Dclose(dataset);

    }



    /* read points for quantities */
    dataset = H5Dopen2(file, qty_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,quantity);
    H5Dclose(dataset);

    /* read points for indices */

    dataset = H5Dopen2(file, index_name, H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_quantity);
    H5Dclose(dataset);

    /* finally closing the file */
    H5Fclose(file);

}

/* for reading the quadruples we need more data */

void hdf5_read_av_(double * nb, double * t, double * y_q, double * yav,double * aav, double * zav, int * index_quantity,int * openfile)
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    herr_t      status;

    /* open existing file */    
    file = H5Fopen (file_name, H5F_ACC_RDONLY,H5P_DEFAULT);

    if(*openfile == 1) { 

    /* read points for density */

    dataset = H5Dopen2(file, "nb", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nb);
    H5Dclose(dataset);

    /* read points for temperature */

    dataset = H5Dopen2(file, "t", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
    H5Dclose(dataset);

    /* read points for charge fraction */

    dataset = H5Dopen2(file, "yq", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y_q);
    H5Dclose(dataset);

    }



    /* read points for yav */
    dataset = H5Dopen2(file, "yav", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,yav);
    H5Dclose(dataset);

    /* read points for aav */
    dataset = H5Dopen2(file, "aav", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,aav);
    H5Dclose(dataset);

    /* read points for zav */
    dataset = H5Dopen2(file, "zav", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,zav);
    H5Dclose(dataset);

    /* read points for indices */

    dataset = H5Dopen2(file, "index_av", H5P_DEFAULT);
    status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_quantity);
    H5Dclose(dataset);

    /* finally closing the file */
    H5Fclose(file);

}
