#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <ctime>
#include <malloc.h>
#include </usr/include/xdrfile/xdrfile_trr.h> // xdr include file 
#include </usr/include/xdrfile/xdrfile_xtc.h> // xdr include file 

    using namespace std;

int main(int argc, char **argv){
	
	/*variables*/
	// XTC variables
	XDRFILE *xd;		// the xtc file
	int natoms;	// number of total atoms
	int step;		// the step counter for the simulation
	float time;		// simulation time
	matrix box;		// box coordinates in a 3x3 matrix
	rvec *coor;		// atom coordinates in a 2-D matrix
	rvec *vel;		// atom coordinates in a 2-D matrix
	rvec *force;		// atom coordinates in a 2-D matrix
//	float prec;	
	float lambda;
	
	int firstframe = 500;				
	int lastframe ;		// ps

	//other variables
	
	char f[40];
	FILE *data;
	int i,j,k;
	
	double *rhoX,*rhoXtmp;
	double *rhoY,*rhoYtmp;
	double *rhoZ,*rhoZtmp;

	double z;
	//these are average values	
	double dx = 0.1;
	int xbin,ybin,zbin;
	
	double zs;
	double comZ;
	double integral;
	double indicator1,indicator2;

	/*read xtc files*/
	read_trr_natoms(argv[1],&natoms);

	coor = (rvec *)malloc(natoms*sizeof(rvec));
	vel = (rvec *)malloc(natoms*sizeof(rvec));
	force = (rvec *)malloc(natoms*sizeof(rvec));


	sprintf(f,"rho_z_com_new.dat");
	data = fopen(f,"w");
	fprintf(data,"# z\tRhoz\tint Rhoz\n");

	//open xtc file and loop through each frame
	xd=xdrfile_open(argv[1],"r");
	k = 0;
	while( ! read_trr(xd, natoms, &step, &time, &lambda, box, coor, vel, force)){
		if(coor == 0){
			printf("Insufficient memory to load .trr file. \n");
			return 1;
		}
		if(step == 0){
			xbin = int(box[0][0]/dx);
			ybin = int(box[1][1]/dx);
			zbin = int(box[2][2]/dx);

			printf("# of atoms = %d\n# Xbin = %d\n# Ybin = %d\n# Zbin = %d\n",natoms,xbin,ybin,zbin);

			rhoX = (double *)malloc(sizeof(double)*xbin);
			rhoY = (double *)malloc(sizeof(double)*ybin);
			rhoZ = (double *)malloc(sizeof(double)*zbin);

			rhoXtmp = (double *)malloc(sizeof(double)*xbin);
			rhoYtmp = (double *)malloc(sizeof(double)*ybin);
			rhoZtmp = (double *)malloc(sizeof(double)*zbin);

			for(i=0;i<xbin;i++) rhoX[i] = 0.0;
			for(i=0;i<ybin;i++) rhoY[i] = 0.0;
			for(i=0;i<zbin;i++) rhoZ[i] = 0.0;
			printf("# natoms / zbin = %d\n",natoms/zbin);

		}

	    if(time >= firstframe ){

			// get Zcom
			comZ = 0.0;
			for(i=0;i<natoms;i++) comZ += coor[i][2]/double(natoms);
			if(abs(comZ - box[2][2]/2.0) > box[2][2]/double(zbin)){
				// get rho'(z)
				for(i=0;i<zbin;i++) rhoZtmp[i] = 0.0;
				for(i=0;i<natoms;i++){
					rhoZtmp[int(coor[i][2]/dx)] += 1.0/double(natoms)/dx;
				}
				
				indicator2 = 0.0;	
				for(i=0;i< (zbin*10 );i++){
					zs = i*dx/10.0;
					integral = 0.0;
					for(j=0;j *dx < zs;j++)	integral +=	(rhoZtmp[j]+rhoZtmp[j+1])/2.0*dx;
					integral += rhoZtmp[j]*(zs-j*dx);
					indicator1 = integral - zs / box[2][2] + comZ / box[2][2] - 0.50;
					if( i > 0 && indicator2*indicator1 <= 0 && (zs-box[2][2]/2.0)*(comZ - box[2][2]/2.0) <= 0 && rhoZtmp[int(zs/dx)] < 1.0/double(zbin)/dx ){
						zs = dx/10.0 / (abs(indicator1)+abs(indicator2))*abs(indicator2) + (i-1)*dx/10.0;
						break;
					}
					indicator2 = indicator1;
			//		printf("%d\t%e\t%e\t%e\n",i,zs,integral,indicator1);
					if(i == zbin*10 - 1 && indicator2*indicator1 > 0 ) printf("# Error when remove c.o.m...%f.(%f)..\n",comZ,box[2][2]/2.0);
				}
				//printf("# zs = %e\n# ComZ = %e\n# rhoZ[zs] = %f\n",zs,comZ,rhoZtmp[int(zs/dx)]);
			}
			else zs = 0;
		
			for(i=0;i<natoms;i++){
				if(coor[i][2] - zs >= 0) rhoZ[int((coor[i][2] - zs )/dx)] += 1.0;
				else rhoZ[int((coor[i][2] - zs + box[2][2])/dx)] += 1.0;
			}
			k++;
	    }
	}
	for(i=0;i<zbin;i++) rhoZ[i] = rhoZ[i]/double(k);
	
	printf("# Finish reading .trr file...\n");

	double nCheck;	// integral of rhoZ to check the area under the curve	

	double rhoMin,rhoMax;	// get the max and min value of rho
	nCheck = 0.0;
	rhoMin = rhoZ[0];
	rhoMax = rhoZ[0];
	for(i=0;i<zbin;i++){
		if(i<zbin-1) nCheck += (rhoZ[i] + rhoZ[i+1])/2.0*double(dx);
		fprintf(data,"%.3f\t%e\t%e\n",i*dx,rhoZ[i],nCheck);
		if(rhoMin > rhoZ[i]) rhoMin = rhoZ[i];
		if(rhoMax < rhoZ[i] ) rhoMax = rhoZ[i];
	}

	fprintf(data,"# rhoMax = %f\n# rhoMin = %f\n",rhoMax,rhoMin);

	/*	get the average max and min
 *	average from max to 0.99 max is the average of max
 *	average from min to 1.01 min is the average of min
 */
	double rhoMaxAverage = 0.0;
	double rhoMinAverage = 0.0;
	int count1 = 0;
	int count2 = 0;

	for(i=0;i<zbin;i++){
		if(rhoZ[i] >= 0.99* rhoMax){
			rhoMaxAverage += rhoZ[i];
			count1++;
		}
		if(rhoZ[i] <= 1.01 * rhoMin){
			rhoMinAverage += rhoZ[i];
			count2++;
		}
	}

	fprintf(data,"# <rhoMax> = %f\n# <rhoMin> = %f\n",rhoMaxAverage/double(count1),rhoMinAverage/double(count2));
	fprintf(data,"# average from max - 0.99 max, min - 1.01 min\n");

	fclose(data);

	return 0;
}
