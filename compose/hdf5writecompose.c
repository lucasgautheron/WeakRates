/*
 *  This function writes data to the HDF5 file.
 */

#define H5_USE_16_API
#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RANK1   1
#define file_name "eoscompose.h5"

void hdf5_write_thermo_(int* m_nb, int*n_temp, int* o_y_q, int * n_qty,double * nb, double * t, double * y_q, double * thermo_hdf5, int * index_thermo,int * openfile, int* itab )
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hsize_t     dimsf[1];/* dataset dimensions */
    hsize_t     dimsfthermo[1];/* dataset dimensions */
    hsize_t     dimsfindex[1];/* dataset dimensions */
    herr_t      status;
    int mrho,nt,oye;

    mrho = *m_nb;
    nt = *n_temp;
    oye = *o_y_q;
    dimsf[0] = 1;
    dimsfindex[0] = *n_qty;
    dimsfthermo[0] = mrho*nt*oye*(*n_qty);


    if(*openfile == 1) { 
    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */     
      file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* points for n_b */
    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsnb", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);
    
/*      * Close/release resources.
     */
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* store the value of itab for later reading purposes */

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "itab", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, itab);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* points for T */
    /* if itab =0, then data points are common for nb, t, yq*/

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointst", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab ==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_temp);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* points for y_q */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsyq", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, o_y_q);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for density */

    dimsf[0] = *m_nb;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "nb", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nb);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for temperature */

    if(*itab==0) {
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *n_temp;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "t", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for charge fraction */

    if(*itab ==0 ){
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *o_y_q;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "yq", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y_q);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    }

    if(*openfile == 2) { 
      /*
	/* open existing file */
      file = H5Fopen (file_name, H5F_ACC_RDWR,H5P_DEFAULT);}


    /* write number of different thermodynamic quantities */

    dimsf[0] = 1;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsqty", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_qty);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    if(*n_qty != 0) {

    /* write points for thermo quantities */
    dataspace = H5Screate_simple(RANK1, dimsfthermo, NULL);
    dataset = H5Dcreate2(file, "thermo", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, thermo_hdf5);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for thermo indices */

    dataspace = H5Screate_simple(RANK1, dimsfindex, NULL);
    dataset = H5Dcreate2(file, "index_thermo", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_thermo);
    H5Dclose(dataset);
    H5Sclose(dataspace);}

    /* finally closing the file */
    H5Fclose(file);

}


void hdf5_write_thermo_add_(int* m_nb, int*n_temp, int* o_y_q, int * n_qty,double * nb, double * t, double * y_q, double * thermo_hdf5, int * index_thermo,int * openfile, int* itab )
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hsize_t     dimsf[1];/* dataset dimensions */
    hsize_t     dimsfthermo[1];/* dataset dimensions */
    hsize_t     dimsfindex[1];/* dataset dimensions */
    herr_t      status;
    int mrho,nt,oye;

    mrho = *m_nb;
    nt = *n_temp;
    oye = *o_y_q;
    dimsf[0] = 1;
    dimsfindex[0] = *n_qty;
    dimsfthermo[0] = mrho*nt*oye*(*n_qty);


    if(*openfile == 1) { 
    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */     
      file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* points for n_b */
    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsnb", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);
    
/*      * Close/release resources.
     */
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* store the value of itab for later reading purposes */

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "itab", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, itab);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* points for T */
    /* if itab =0, then data points are common for nb, t, yq*/

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointst", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab ==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_temp);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* points for y_q */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsyq", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, o_y_q);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for density */

    dimsf[0] = *m_nb;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "nb", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nb);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for temperature */

    if(*itab==0) {
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *n_temp;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "t", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for charge fraction */

    if(*itab ==0 ){
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *o_y_q;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "yq", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y_q);
    H5Dclose(dataset);
    H5Sclose(dataspace);

}

    if(*openfile == 2) { 
      /*
	/* open existing file */
      file = H5Fopen (file_name, H5F_ACC_RDWR,H5P_DEFAULT);}


    /* write number of different thermodynamic quantities */

    dimsf[0] = 1;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsadd", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_qty);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    if(*n_qty !=0) {
    /* write points for thermo quantities */
    dataspace = H5Screate_simple(RANK1, dimsfthermo, NULL);
    dataset = H5Dcreate2(file, "thermo_add", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, thermo_hdf5);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for thermo indices */

    dataspace = H5Screate_simple(RANK1, dimsfindex, NULL);
    dataset = H5Dcreate2(file, "index_thermo_add", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_thermo);
    H5Dclose(dataset);
    H5Sclose(dataspace);}

    /* finally closing the file */
    H5Fclose(file);

}



void hdf5_write_compo_p_(int* m_nb, int*n_temp, int* o_y_q, int * n_qty,double * nb, double * t, double * y_q, double * thermo_hdf5, int * index_thermo,int * openfile, int* itab )
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hsize_t     dimsf[1];/* dataset dimensions */
    hsize_t     dimsfthermo[1];/* dataset dimensions */
    hsize_t     dimsfindex[1];/* dataset dimensions */
    herr_t      status;
    int mrho,nt,oye;

    mrho = *m_nb;
    nt = *n_temp;
    oye = *o_y_q;
    dimsf[0] = 1;
    dimsfindex[0] = *n_qty;
    dimsfthermo[0] = mrho*nt*oye*(*n_qty);


    if(*openfile == 1) { 
    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */     
      file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* points for n_b */
    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsnb", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);
    
/*      * Close/release resources.
     */
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* store the value of itab for later reading purposes */

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "itab", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, itab);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* points for T */
    /* if itab =0, then data points are common for nb, t, yq*/

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointst", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab ==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_temp);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* points for y_q */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsyq", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, o_y_q);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for density */

    dimsf[0] = *m_nb;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "nb", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nb);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for temperature */

    if(*itab==0) {
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *n_temp;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "t", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for charge fraction */

    if(*itab ==0 ){
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *o_y_q;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "yq", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y_q);
    H5Dclose(dataset);
    H5Sclose(dataspace);

}

    if(*openfile == 2) { 
      /*
	/* open existing file */
      file = H5Fopen (file_name, H5F_ACC_RDWR,H5P_DEFAULT);}


    /* write number of different pairs */

    dimsf[0] = 1;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointspairs", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_qty);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    
    if(*n_qty != 0 ) {

    /* write points for particle fractions*/
    dataspace = H5Screate_simple(RANK1, dimsfthermo, NULL);
    dataset = H5Dcreate2(file, "yi", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, thermo_hdf5);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for particle indices */

    dataspace = H5Screate_simple(RANK1, dimsfindex, NULL);
    dataset = H5Dcreate2(file, "index_yi", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_thermo);
    H5Dclose(dataset);
    H5Sclose(dataspace);}

    /* finally closing the file */
    H5Fclose(file);

}


void hdf5_write_compo_q_(int* m_nb, int*n_temp, int* o_y_q, int * n_q,double * nb, double * t, double * y_q, double * yav,double * aav, double * zav, int * index_av,int * openfile, int* itab )
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hsize_t     dimsf[1];/* dataset dimensions */
    hsize_t     dimsfthermo[1];/* dataset dimensions */
    hsize_t     dimsfindex[1];/* dataset dimensions */
    herr_t      status;
    int mrho,nt,oye;

    mrho = *m_nb;
    nt = *n_temp;
    oye = *o_y_q;
    dimsf[0] = 1;
    dimsfindex[0] = *n_q;
    dimsfthermo[0] = mrho*nt*oye*(*n_q);


    if(*openfile == 1) { 
    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */     
      file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* points for n_b */
    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsnb", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);
    
/*      * Close/release resources.
     */
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* store the value of itab for later reading purposes */

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "itab", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, itab);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* points for T */
    /* if itab =0, then data points are common for nb, t, yq*/

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointst", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab ==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_temp);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* points for y_q */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsyq", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, o_y_q);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for density */

    dimsf[0] = *m_nb;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "nb", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nb);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for temperature */

    if(*itab==0) {
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *n_temp;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "t", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for charge fraction */

    if(*itab ==0 ){
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *o_y_q;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "yq", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y_q);
    H5Dclose(dataset);
    H5Sclose(dataspace);

}

    if(*openfile == 2) { 
      /*
	/* open existing file */
      file = H5Fopen (file_name, H5F_ACC_RDWR,H5P_DEFAULT);}


    /* write number of different thermodynamic quantities */

    dimsf[0] = 1;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsav", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_q);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    if(*n_q != 0 ) {

    /* write points for average fractions*/
    dataspace = H5Screate_simple(RANK1, dimsfthermo, NULL);
    dataset = H5Dcreate2(file, "yav", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, yav);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for average baryon number*/
    dataspace = H5Screate_simple(RANK1, dimsfthermo, NULL);
    dataset = H5Dcreate2(file, "aav", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, aav);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for average charge */
    dataspace = H5Screate_simple(RANK1, dimsfthermo, NULL);
    dataset = H5Dcreate2(file, "zav", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, zav);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for average indices */

    dataspace = H5Screate_simple(RANK1, dimsfindex, NULL);
    dataset = H5Dcreate2(file, "index_av", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_av);
    H5Dclose(dataset);
    H5Sclose(dataspace);}

    /* finally closing the file */
    H5Fclose(file);

}


void hdf5_write_micro_(int* m_nb, int*n_temp, int* o_y_q, int * n_qty,double * nb, double * t, double * y_q, double * thermo_hdf5, int * index_thermo,int * openfile, int* itab )
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hsize_t     dimsf[1];/* dataset dimensions */
    hsize_t     dimsfthermo[1];/* dataset dimensions */
    hsize_t     dimsfindex[1];/* dataset dimensions */
    herr_t      status;
    int mrho,nt,oye;

    mrho = *m_nb;
    nt = *n_temp;
    oye = *o_y_q;
    dimsf[0] = 1;
    dimsfindex[0] = *n_qty;
    dimsfthermo[0] = mrho*nt*oye*(*n_qty);


    if(*openfile == 1) { 
    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */     
      file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* points for n_b */
    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsnb", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);
    
/*      * Close/release resources.
     */
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* store the value of itab for later reading purposes */

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "itab", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, itab);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* points for T */
    /* if itab =0, then data points are common for nb, t, yq*/

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointst", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab ==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_temp);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* points for y_q */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsyq", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, o_y_q);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for density */

    dimsf[0] = *m_nb;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "nb", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nb);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for temperature */

    if(*itab==0) {
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *n_temp;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "t", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for charge fraction */

    if(*itab ==0 ){
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *o_y_q;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "yq", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y_q);
    H5Dclose(dataset);
    H5Sclose(dataspace);

}

    if(*openfile == 2) { 
      /*
	/* open existing file */
      file = H5Fopen (file_name, H5F_ACC_RDWR,H5P_DEFAULT);}


    /* write number of different microscopic quantities*/

    dimsf[0] = 1;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsmicro", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_qty);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    if(*n_qty != 0 ) {

    /* write points for microscopic quantities*/
    dataspace = H5Screate_simple(RANK1, dimsfthermo, NULL);
    dataset = H5Dcreate2(file, "micro", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, thermo_hdf5);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for indices for microscopic quantities */

    dataspace = H5Screate_simple(RANK1, dimsfindex, NULL);
    dataset = H5Dcreate2(file, "index_micro", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_thermo);
    H5Dclose(dataset);
    H5Sclose(dataspace);}

    /* finally closing the file */
    H5Fclose(file);

}



void hdf5_write_err_(int* m_nb, int*n_temp, int* o_y_q, int * n_err,double * nb, double * t, double * y_q, double * err_hdf5, int * index_err,int * openfile, int* itab )
{
    hid_t       file, dataset;         /* file and dataset handles */
    hid_t       datatype, dataspace;   /* handles */
    hsize_t     dimsf[1];/* dataset dimensions */
    hsize_t     dimsfthermo[1];/* dataset dimensions */
    hsize_t     dimsfindex[1];/* dataset dimensions */
    herr_t      status;
    int mrho,nt,oye;

    mrho = *m_nb;
    nt = *n_temp;
    oye = *o_y_q;
    dimsf[0] = 1;
    dimsfindex[0] = *n_err;
    dimsfthermo[0] = mrho*nt*oye*(*n_err);


    if(*openfile == 1) { 
    /*
     * Create a new file using H5F_ACC_TRUNC access,
     * default file creation properties, and default file
     * access properties.
     */     
      file = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* points for n_b */
    /*
     * Create a new dataset within the file using defined dataspace and
     * datatype and default dataset creation properties.
     */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsnb", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Write the data to the dataset using default transfer properties.
     */
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);
    
/*      * Close/release resources.
     */
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* store the value of itab for later reading purposes */

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "itab", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, itab);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* points for T */
    /* if itab =0, then data points are common for nb, t, yq*/

    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointst", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab ==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_temp);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* points for y_q */
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointsyq", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(*itab==0){
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m_nb);}
    else {
      status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, o_y_q);}
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for density */

    dimsf[0] = *m_nb;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "nb", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nb);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for temperature */

    if(*itab==0) {
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *n_temp;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "t", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, t);
    H5Dclose(dataset);
    H5Sclose(dataspace);


    /* write points for charge fraction */

    if(*itab ==0 ){
      dimsf[0] = *m_nb;}
    else {
      dimsf[0] = *o_y_q;}
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "yq", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y_q);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    }

    if(*openfile == 2) { 
      /*
	/* open existing file */
      file = H5Fopen (file_name, H5F_ACC_RDWR,H5P_DEFAULT);}


    /* write number of different thermodynamic quantities */

    dimsf[0] = 1;
    dataspace = H5Screate_simple(RANK1, dimsf, NULL);
    dataset = H5Dcreate2(file, "pointserr", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, n_err);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    if(*n_err != 0 ) {

    /* write points for error quantities */
    dataspace = H5Screate_simple(RANK1, dimsfthermo, NULL);
    dataset = H5Dcreate2(file, "error", H5T_NATIVE_DOUBLE, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, err_hdf5);
    H5Dclose(dataset);
    H5Sclose(dataspace);

    /* write points for error indices */

    dataspace = H5Screate_simple(RANK1, dimsfindex, NULL);
    dataset = H5Dcreate2(file, "index_err", H5T_NATIVE_INT, dataspace,
			H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, index_err);
    H5Dclose(dataset);
    H5Sclose(dataspace);}

    /* finally closing the file */
    H5Fclose(file);

}


