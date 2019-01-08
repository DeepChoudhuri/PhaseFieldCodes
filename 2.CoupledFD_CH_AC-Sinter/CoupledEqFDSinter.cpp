#include "CoupledEqFDSinter.h"
#include "PFUtilities.h"

//Funtions related to microstrucutre creation
void prepareInitialMicrostructure(double **Cxy, std::vector<double **> Evec,\
                                  const computeParam *cp, \
                                  const materialParam *mp)
{
  std::cout << "Preparing initial microstructutre ..." << std::endl;

  //Coordinates of particle - dynamically assigned
  double *xc;
  double *yc;

  // used when more than two etas or particles
  double RN = 10.0;
  double RNnew = RN;

  // used when more than two etas or particles
  double R1 = 20.0;
  double R2 = 0.5*R1;

  // For evaluating is the point falls within a particle
  double rxy1 = 0;
  double rxy2 = 0; // This one is used for bi-crystal configuration


  if(mp->etaNum !=2){ // Here etaNum is treated as number of particles
    xc = new double[mp->etaNum];
    yc = new double[mp->etaNum];

    //Equilaterial Triangle - 3 particles of equal size
    RNnew = 20.0;
    xc[0] = 30.0; yc[0] = 30.0;
    xc[1] = 70.0; yc[1] = 30.0;
    xc[2] = 50.0; yc[2] = 64.64;

    /*//Square - 9 particles
    xc[0] = 29.0; yc[0] = 50.0;
    xc[1] = 50.0; yc[1] = 50.0;
    xc[2] = 71.0; yc[2] = 50.0;
    xc[3] = 50.0; yc[3] = 29.0;
    xc[4] = 50.0; yc[4] = 71.0;
    xc[5] = 39.0; yc[5] = 39.0;
    xc[6] = 61.0; yc[6] = 39.0;
    xc[7] = 39.0; yc[7] = 61.0;
    xc[8] = 61.0; yc[8] = 61.0;*/

    for(int ei = 0; ei < mp->etaNum; ++ei){
      /*if(ei > 4){
        RNnew = 0.5*RN;
      }*/

      for(int i = 0; i < cp->Nx; ++i){
        for(int j = 0; j < cp->Ny; ++j){
          rxy1 = std::sqrt(std::pow((i-xc[ei]),2.0) + std::pow((j-yc[ei]),2.0));
          if(rxy1 <= RNnew){ //If the distance rx1 is less than the Radius
            Cxy[i][j] = 0.999; // Assign max. values for concentration
            Evec[ei][i][j] = 0.999; // Assign max. values for etas for each particle
          }
        }
      }
    }
  } // This ends microstrucutre with more than 2 particles

  if(mp->etaNum == 2){ // When we have a Bi-crystal, i.e. only TWO particles
    xc = new double[mp->etaNum];
    yc = new double[mp->etaNum];
    //Coordinates assignment
    xc[0] = cp->Nx*0.5;
    yc[0] = 40.0;
    yc[1] = 70.0;

    for(int i = 0; i < cp->Nx; ++i){
      for (int j = 0; j < cp->Ny; ++j) {
        rxy1 = std::sqrt(std::pow((i-xc[0]),2.0) + std::pow((j-yc[0]),2.0));
        rxy2 = std::sqrt(std::pow((i-xc[0]),2.0) + std::pow((j-yc[1]),2.0));

        if(rxy1 <= R1){ //If the distance rx1 is less than the Radius
          Cxy[i][j] = 0.999; // Assign max. values for concentration
          Evec[0][i][j] = 0.999; // Assign max. values for etas for each particle
        }

        if(rxy2 <= R2){ //If the distance rx1 is less than the Radius
          Cxy[i][j] = 0.999; // Assign max. values for concentration
          Evec[1][i][j] = 0.999; // Assign max. values for etas for each particle
        }
      }
    }
    //xc = NULL; yc = NULL;
  } //End bi-crystal configuration
  std::cout << "Intial microstructure created!" << std::endl;
  //Garbage cleaning
  delete[] xc;
  delete[] yc;
}//END function



//This function does the following:
// 1. Takes all the Eta data from the Vector and compresses them
//    to a single writable matrix, and
// 2. Writes the concetration file also to a the same file as the Etas.
void writeAllDataToVTKFile(double ** cxy, std::vector<double **> evec,\
                           const computeParam *cp, const materialParam *mp,\
                           long istep)
{
  double **exy2 = create2DField(cp);

  for(int i = 0; i < cp->Nx; ++i){
    for(int j = 0; j < cp->Ny; ++j){
      for(int ne = 0; ne < mp->etaNum; ++ne){
        exy2[i][j] = exy2[i][j] + std::pow(evec[ne][i][j], 2.0);
      }
    }
  }
  write2DVTKfile(cxy, exy2, cp, istep);
}//END function



//Derivative of the Free energy with respect to the concentration.
////Takes in one concentration value,
//and MULTIPLE non-conservative (Eta) parameter fields
double dFdCon(double cxy, std::vector<double **> evec,\
              const materialParam *mp,long p, long q)
{
  double sumei2 = 0.0;
  double sumei3 = 0.0;
  for(int ei = 0; ei < mp->etaNum; ++ei){
    sumei2 = sumei2 + std::pow(evec[ei][p][q],2.0);
    sumei3 = sumei3 + std::pow(evec[ei][p][q],3.0);
  }
  return 2.0*mp->A*(cxy*(1.0-cxy)*(1.0-cxy) - cxy*cxy*(1.0-cxy)) + \
         mp->B*(2.0*cxy-6.0*sumei2+4.0*sumei3);
//         2.0*mp->A*(1.0 - cxy)*(1.0 - 2.0*cxy);
//         2.0*mp->A*(cxy*std::pow(1.0-cxy,2.0) - std::pow(cxy,2.0)*(1.0-cxy));

} //END function



//Derivative of the Free energy with respect
//to the non-conserved parameter(s) eta
////Takes in one concentration value,
//and MULTIPLE non-conservative (Eta) parameter fields
double dFdEtai(double cxy, double eixy, std::vector<double **> evec,\
               const materialParam *mp, long p, long q)
{
  double sumei2 = 0.0;
  for(int i=0; i < mp->etaNum; ++i){
    sumei2 = sumei2 + std::pow(evec[i][p][q],2.0);
  }
  return 12.0*mp->B*((1.0-cxy)*eixy -\
         (2.0-cxy)*std::pow(eixy,2.0) +\
         eixy*sumei2);//Expression was analytically derived
}//END function



//Total  bulk free energy at a given time step
//Takes in ONE concentration field,
//and MULTIPLE non-conservative (Eta) parameter fields
double BulkFreeEnergy(double **cxy, std::vector<double **> evec,\
                      const computeParam *cp,const materialParam *mp)
{
  double totFreeEng = 0.0;
  double sumei2, sumei3;

  for(int i = 0; i < cp->Nx; ++i){
    for(int j = 0; j < cp->Ny; ++j){

      sumei2 = 0.0;sumei3 = 0.0;
      for(int ei = 0; ei < mp->etaNum; ++ei){
        sumei2 = sumei2 + pow(evec[ei][i][j],2.0);
        sumei3 = sumei3 + pow(evec[ei][i][j],3.0);
      }

      totFreeEng = totFreeEng + mp->A*std::pow((cxy[i][j]*(1.0-cxy[i][j])),2.0)\
                  + mp->B*(std::pow(cxy[i][j], 2.0) + 6.0*(1.0-cxy[i][j])*sumei2\
                  - 4.0*(2.0-cxy[i][j])*sumei3 + 3.0*std::pow(sumei2,2.0));
    }
  }
  return totFreeEng;
}//END function



// Concentration dependent mobility
double mobility(double cxy, std::vector<double **> evec,\
              const materialParam *mp,long p, long q)
{
  double mob = 0.0;
  double phi = std::pow(cxy,3.0)*(10.0-15.0*cxy+6.0*std::pow(cxy,2.0));
  double sumProdij = 0.0;
  for(int i = 0; i < mp->etaNum; ++i){
    for(int j = 0; j < mp->etaNum; ++j){
      if(i!=j){
        sumProdij = sumProdij + evec[i][p][q]*evec[j][p][q];
      }
    }
  }
  mob = mp->dvol*phi + \
        mp->dvap*(1.0-phi) + \
        mp->dsur*cxy*(1.0-cxy) +\
        mp->dgb*sumProdij;

  return mob;
}



// Solve differential equations to evolve concentration and the Eta fiels
void evolveMicrostrucutre(double **Cxy, std::vector<double **>Evec,\
                          const computeParam *cp,const materialParam *mp)
{
  std::cout << "Solving CH and AC coupled differential equations ... "<< std::endl;

  // Stores bulk Free Energy at each time step
  double *FBulk = new double[cp->nstep];
  // For storing time data
  double *step = new double[cp->nstep];
  double *time = new double[cp->nstep]; // To be used for plotting F vs. t
  double tt = 0; // time variable

  //For concentration fields for CH equations
  double **lapConc; // Laplacian of concetration
  double **idFdc = create2DField(cp); // Matrix evaluating inner derivative
  double **lapIdFdc; // Laplacian of the inner derivative

  //For non-conservative ith eta field for AC equations
  double **lapEtai; // Laplacian of ith eta
  double idFdetai; // Evaluation of inner derivative
  double **etaxyi  = create2DField(cp); // ith eta matrix for temporary storage


  // START time integration
  for(int t = 0; t < cp->nstep; ++t){
    step[t] = t;
    time[t] = tt;

  /*-------------------Evolve concentration fields-----------------------------*/
    lapConc = get2DLaplacian(Cxy, cp);// Laplacian of the concentration field

    //START loop for evaluating DF/Dcon - k.Laplacian(con)
    for(int i = 0; i < cp->Nx; i++){
      for(int j = 0; j < cp->Ny; j++){
        idFdc[i][j] =  dFdCon(Cxy[i][j],Evec,mp,i,j) - 0.5*mp->Kc*lapConc[i][j];
      }
     }//END loop for evaluating DF/Dcon - k.Laplacian(con)

     lapIdFdc = get2DLaplacian(idFdc,cp); //Laplacian of th einner derivative
     //START time integration loop for concentration
     for(int i = 0; i < cp->Nx; ++i){
       for(int j = 0; j < cp->Ny; ++j){

         Cxy[i][j] = Cxy[i][j] + mobility(Cxy[i][j],Evec,mp,i,j)\
                                 *cp->dtime*lapIdFdc[i][j];

         if(Cxy[i][j] >= 0.9999){ //Set MAX limits for small variations
           Cxy[i][j] = 0.9999;
         }

         if(Cxy[i][j] < 0.00001){ //Set MIN limits for small variations
           Cxy[i][j] = 0.00001;
         }
       }
      }//END time integration loop for concentration
  /*-------------------END evolving concentration fields-------------------*/


  /*-------------------Evolve non-consevatinve Eta fields-------------------*/
    for(int ei = 0; ei < mp->etaNum; ei++){//Run though number of Etas

      for(int i = 0; i < cp->Nx; i++){
        for(int j = 0; j < cp->Ny; j++){
          etaxyi[i][j] = Evec[ei][i][j];
        }
      }

      lapEtai = get2DLaplacian(etaxyi, cp);
      //START time integration loop for etai
      for(int i = 0; i < cp->Nx; i++){
        for(int j = 0; j < cp->Ny; j++){
          idFdetai = dFdEtai(Cxy[i][j],etaxyi[i][j],Evec,mp,i,j) - \
                           0.5*mp->Ke*lapEtai[i][j]; // DF/Detai - k.Laplacian(etai)
          etaxyi[i][j] = etaxyi[i][j] - cp->dtime*mp->L*idFdetai; // time step
          if(etaxyi[i][j] >= 0.9999){ //Set MAX limits for small variations
             etaxyi[i][j] = 0.9999;
          }
          if(etaxyi[i][j] < 0.0001){ //Set MIN limits for small variations
            etaxyi[i][j] = 0.0001;
          }
        }
      } //END time integration loop

      //HARD copying of temporary eta values to the Parent eta data strucutre
      for(int i = 0; i < cp->Nx; i++){
        for(int j = 0; j < cp->Ny; j++){
          Evec[ei][i][j] = etaxyi[i][j];
        }
      }
      //displayArray(Evec[0], cp);
    } // END etas run
  /*-------------------END evolving non-consevatinve Eta fields--------------*/

  if((t%250 == 0) && (t < 20001)){
    writeAllDataToVTKFile(Cxy,Evec,cp,mp,t);
  }
  FBulk[t] = BulkFreeEnergy(Cxy,Evec,cp,mp);
  tt = tt + cp->dtime;
}//END MAIN time integration loop
write1DTXTfile(time, FBulk, cp->nstep, "energy-time.txt");
write1DTXTfile(step, FBulk, cp->nstep, "energy-step.txt");
//writeAllDataToVTKFile(Cxy,Evec,cp,mp,cp->nstep);
std::cout << "Finished solving CH and AC coupled differential equations "<< std::endl;

//Garbage cleaning
delete[] FBulk;
delete[] time;
delete2DArray(lapConc,cp->Nx);
delete2DArray(lapIdFdc,cp->Nx);
delete2DArray(etaxyi,cp->Nx);
}//END function










































//
