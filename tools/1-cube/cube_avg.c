#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct{
 double cell[3];
 int nAtoms;
 int nGrid[3];
 int* atomN;
 double* atomX;
 double* atomY;
 double* atomZ;
 double* grid;
}cube;

void cubeAvg(cube input, int dir);

main(int argc, char *argv[]){
 int i, j, k, dir;
 char cubeFile[200];

 if(argc==3){
  sprintf(cubeFile,"%s",argv[1]);
  if(access(cubeFile,F_OK) == -1){
   printf("File does not exist.\n");
   return 0;
  }
 }
 else{
  printf("Usage: ./cube_avg *.cube avg_dir[x=1|y=2|z=3]\n");
  return 0;
 }
 dir=atoi(argv[2]);

 cube input;

 FILE* fp;
 char str[2048], *str_ptr;
 int flag=0, nline=0, cnt=0, nAtoms=0;

 fp=fopen(cubeFile, "r");

 while(fgets(str,sizeof(str),fp)!=NULL){
  nline++;
  if(nline==3){
   sscanf(str,"%d %*f %*f %*f",&input.nAtoms);
   input.atomN = (int *) malloc (input.nAtoms*sizeof(int));
   input.atomX = (double *) malloc (input.nAtoms*sizeof(double));
   input.atomY = (double *) malloc (input.nAtoms*sizeof(double));
   input.atomZ = (double *) malloc (input.nAtoms*sizeof(double));
  }
  if(nline==4){
   sscanf(str,"%d %lf %*f %*f",&input.nGrid[0], &input.cell[0]);
   input.cell[0]=input.nGrid[0]*input.cell[0];
  }
  if(nline==5){
   sscanf(str,"%d %*f %lf %*f",&input.nGrid[1], &input.cell[1]);
   input.cell[1]=input.nGrid[1]*input.cell[1];
  }
  if(nline==6){
   sscanf(str,"%d %*f %*f %lf",&input.nGrid[2], &input.cell[2]);
   input.cell[2]=input.nGrid[2]*input.cell[2];
  }
  if(nline>=7 && nline<7+input.nAtoms){
    sscanf(str,"%d %*f %lf %lf %lf",&input.atomN[nline-7],&input.atomX[nline-7],&input.atomY[nline-7],&input.atomZ[nline-7]);
  }
  if(nline==input.nAtoms+6){
   printf("### Cell parameters(Bohr): %f %f %f\n", input.cell[0], input.cell[1], input.cell[2]);
   printf("### Number of Atoms: %d\n", input.nAtoms);
   printf("### Grid Dimension: %d %d %d\n", input.nGrid[0], input.nGrid[1], input.nGrid[2]);
   input.grid=(double *)malloc(input.nGrid[0]*input.nGrid[1]*input.nGrid[2]*sizeof(double));
  }
  if(nline>=input.nAtoms+7 && cnt<input.nGrid[0]*input.nGrid[1]*input.nGrid[2]){
   str_ptr=strtok(str," ");
   for(;str_ptr!=NULL;cnt++){
     sscanf(str_ptr,"%lf",&input.grid[cnt]);
    str_ptr=strtok(NULL," ");
   }
  }
 }//parsing cube file end

 cubeAvg(input,dir);

 free(input.atomN);
 free(input.atomX);
 free(input.atomY);
 free(input.atomZ);
 free(input.grid);

}//main END

void cubeAvg(cube input, int dir){ // in-plane average along the selected direction
 FILE* fp;
 int  i, j, k;
 char fname[512];
 double avg=0;

 if(dir==1){
   sprintf(fname,"cube.x.avg");
   fp=fopen(fname,"w");
   for(i=0;i<input.nGrid[0];i++){
    for(j=0;j<input.nGrid[1];j++){
     for(k=0;k<input.nGrid[2];k++){
      avg+=input.grid[k+j*input.nGrid[2]+i*input.nGrid[2]*input.nGrid[1]];
     }
    }
    avg/=input.nGrid[1]*input.nGrid[2];
    fprintf(fp,"%lf %lf\n",i*input.cell[0]/input.nGrid[0],avg);
    avg=0;
   }
  fclose(fp);
 }
 if(dir==2){
   sprintf(fname,"cube.y.avg");
   fp=fopen(fname,"w");
   for(j=0;j<input.nGrid[1];j++){
    for(i=0;i<input.nGrid[0];i++){
     for(k=0;k<input.nGrid[2];k++){
      avg+=input.grid[k+j*input.nGrid[2]+i*input.nGrid[2]*input.nGrid[1]];
     }
    }
    avg/=input.nGrid[0]*input.nGrid[2];
    fprintf(fp,"%lf %lf\n",j*input.cell[1]/input.nGrid[1],avg);
    avg=0;
   }
  fclose(fp);
 }
 if(dir==3){
   sprintf(fname,"cube.z.avg");
   fp=fopen(fname,"w");
   for(k=0;k<input.nGrid[2];k++){
    for(j=0;j<input.nGrid[1];j++){
     for(i=0;i<input.nGrid[0];i++){
      avg+=input.grid[k+j*input.nGrid[2]+i*input.nGrid[2]*input.nGrid[1]];
     }
    }
    avg/=input.nGrid[0]*input.nGrid[1];
    fprintf(fp,"%lf %lf\n",k*input.cell[2]/input.nGrid[2],avg);
    avg=0;
   }
  fclose(fp);
 }
 
 printf("Averaged values were saved: %s\n",fname);
}

