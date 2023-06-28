#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <algorithm>
#include "BM.h"

#define a1 1.0
#define a2 1.4
#define N1 512
#define N2 512
#define Nl 50
#define phi0 0.70
#define phi_max 1.20
#define dphi_fixed 1.0e-4
#define dphi_a 1.0e-4
#define dphi_g 1.0e-4
#define dt0 1.0e-3
#define dt_max 0.1
#define zeta 1.0
#define dim 2
#define polydispersity 0.0
#define cut 1.4
#define skin 1.0

void ini_coord(double (*x)[dim],double L){
  srand((unsigned int)time(NULL));
  for(int i=0;i<N1+N2;i++){
    x[i][0] = (double)rand()/RAND_MAX*L;
    x[i][1] = (double)rand()/RAND_MAX*L;
  }
}

void set_diameter(double *a){
  for(int i=0;i<N1;i++){
    a[i]=a1+polydispersity*gaussian_rand();
  }
  for(int i=0;i<N2;i++){
    a[i+N1]=a2+polydispersity*gaussian_rand();
  }
}

void LE_boundary(double (*x)[dim],double gamma,double L){
  for(int i=0;i<N1+N2;i++){
    if(x[i][0]<x[i][1]*gamma){
      x[i][0]+=L;
    }else if(x[i][0]>x[i][1]*gamma+L){
      x[i][0]-=L;
    }

    if(x[i][1]<0){
      x[i][1]+=L;
      x[i][0]+=L*gamma;
    }else if(x[i][1]>L){
      x[i][1]-=L;
      x[i][0]-=L*gamma;
    }
  }
}

void ini_array(double (*x)[dim]){
  for(int i=0;i<N1+N2;i++)
    for(int j=0;j<dim;j++)
      x[i][j]=0.0;
}

void affine_arithmetic(double (*x)[dim],double *L,double *phi,double dphi){
  int i,j=0;
  // update phi
  *phi+=dphi;
  // transform L and x by using new phi
  *L*=sqrt((*phi-dphi)/(*phi));
  for(i=0;i<N1+N2;i++){
    for(j=0;j<dim;j++){
      x[i][j] *= sqrt((*phi-dphi)/(*phi));
    }
  }
}

void affine_geometric(double (*x)[dim],double *L,double *phi){
  int i,j=0;
  // update phi
  *phi*=(1.0+dphi_g);
  // transform L and x by using old phi
  *L*=sqrt(1.0/(1.0+dphi_g));
  for(i=0;i<N1+N2;i++){
    for(j=0;j<dim;j++){
      x[i][j]*=sqrt(1.0/(1.0+dphi_g));
    }
  }
}

void list_verlet(int (*list)[Nl],double (*x)[dim],double gamma,double L){
  double dx,dy,dr2;
  int i,j;
  for(i=0;i<N1+N2;i++){
    for(j=0;j<Nl;j++){
      list[i][j]=0;
    }
  }
  
  for(i=0;i<N1+N2;i++){
    for(j=0;j<i;j++){
      dx=x[i][0]-x[j][0];
      dy=x[i][1]-x[j][1];
      if(dy>0.5*L){
        dy-=L;
        dx-=L*gamma;
      }else if(dy<-0.5*L){
        dy+=L;
        dx+=L*gamma;
      }
      if(dx>0.5*L){
        dx-=L;
      }else if(dx<-0.5*L){
        dx+=L;
      }
      dr2=dx*dx+dy*dy;
      if(dr2<(cut+skin)*(cut+skin)){
        list[i][0]++;
        list[i][(int)list[i][0]]=j;
      }
    }
  }
}

void update(double (*x_update)[dim],double (*x)[dim]){
  for(int i=0;i<N1+N2;i++)
    for(int j=0;j<dim;j++)
      x_update[i][j]=x[i][j];
}

void calc_disp_max(double *disp_max,double (*x)[dim],double (*x_update)[dim],double gamma,double L){
  double dx,dy;
  double disp;
  for(int i=0;i<N1+N2;i++){
    dx=x[i][0]-x_update[i][0];
    dy=x[i][1]-x_update[i][1];
    //dx-=L*floor((dx+0.5*L)/L);
    //dy-=L*floor((dy+0.5*L)/L);
    if(dy>0.5*L){
        dy-=L;
        dx-=L*gamma;
      }else if(dy<-0.5*L){
        dy+=L;
        dx+=L*gamma;
      }
      if(dx>0.5*L){
        dx-=L;
      }else if(dx<-0.5*L){
        dx+=L;
      }
    disp = dx*dx+dy*dy;
    if(disp > *disp_max){
      *disp_max =disp;
    }
  }
  //std::cout<<"disp_max is"<<"\t"<<*disp_max<<"\n";
}

void auto_list_update(double *disp_max,double (*x)[dim],double (*x_update)[dim],int (*list)[Nl],double gamma,double L){
  static int count=0;
  count++;
  calc_disp_max(&(*disp_max),x,x_update,gamma,L);
  if(*disp_max > 0.5*skin*skin*0.25){
    list_verlet(list,x,gamma,L);
    update(x_update,x);
    //std::cout<<"update"<<*disp_max<<" "<<count<<std::endl;
    *disp_max=0.0;
    count=0;
  }
}

void calc_force(int (*list)[Nl],double (*x)[dim],double (*f)[dim],double *a,double *U,double *p,double *sxy,double gamma,double L){
  double dx,dy,dr,dr2,dUr,aij;
  int i,j;

  *U=0.0;
  *p=0.0;
  *sxy=0.0;
  ini_array(f);
  
  // calculate force
  for(i=0;i<N1+N2;i++){
    for(j=1;j<=list[i][0];j++){
      dx=x[i][0]-x[list[i][j]][0];
      dy=x[i][1]-x[list[i][j]][1];
      if(dy>0.5*L){
        dy-=L;
        dx-=L*gamma;
      }else if(dy<-0.5*L){
        dy+=L;
        dx+=L*gamma;
      }
      if(dx>0.5*L){
        dx-=L;
      }else if(dx<-0.5*L){
        dx+=L;
      }
      dr2=dx*dx+dy*dy;
      dr=sqrt(dr2);
      aij=0.5*(a[i]+a[list[i][j]]);

      if(dr<aij){
        dUr=-(1.0-dr/aij)/aij;
        f[i][0]-=dUr*dx/dr;
        f[list[i][j]][0]+=dUr*dx/dr;
        f[i][1]-=dUr*dy/dr;
        f[list[i][j]][1]+=dUr*dy/dr;

        *U+=(1.0-dr/aij)*(1.0-dr/aij)/2.0/(double)(N1+N2);  // alpha=2 at Hertzian potential
        *p-=dr*dUr/(2.0*L*L);
        *sxy-=dx*dy*dUr/(2.0*L*L);
      }
    //std::cout<<i<<"\t"<<list[i][j]<<"\t"<<aij-dr<<"\t"<<f[i][0]<<"\t"<<f[i][1]<<std::endl;
    }
  }
}

void apply_shear(int (*list)[Nl],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *sxy,double *disp_max,double gamma,double dgamma,double L){
  for(int i=0;i<N1+N2;i++){
    x[i][0]+=x[i][1]*dgamma;
  }
  auto_list_update(disp_max,x,x_update,list,gamma,L);
  calc_force(list,x,f,a,&(*U),&(*p),&(*sxy),gamma,L);
}

void FIRE(int (*list)[Nl],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *sxy,double *f_tot,double *disp_max,double gamma,double L){
  double v[N1+N2][dim],P,v_norm,f_norm;
  double dt=dt0,A=0.1;
  int count=0;

  // initialize
  ini_array(v);
  list_verlet(list,x,gamma,L);
  calc_force(list,x,f,a,&(*U),&(*p),&(*sxy),gamma,L);

  for(;;){
    // initialize
    P=0.0;
    *f_tot=0.0;
    // 1. MD integration with Velocity Verlet
    for(int i=0;i<N1+N2;i++){
      for(int j=0;j<dim;j++){
        x[i][j]+=v[i][j]*dt+0.5*f[i][j]*dt*dt;
        v[i][j]+=0.5*f[i][j]*dt;
      }
    }
    calc_force(list,x,f,a,&(*U),&(*p),&(*sxy),gamma,L);
    for(int i=0;i<N1+N2;i++){
      for(int j=0;j<dim;j++){
        v[i][j]+=0.5*f[i][j]*dt;
      }
    }
    LE_boundary(x,gamma,L);
    auto_list_update(&(*disp_max),x,x_update,list,gamma,L);

    // 3. calculate power P
    for(int i=0;i<N1+N2;i++){
      v_norm=sqrt(v[i][0]*v[i][0]+v[i][1]*v[i][1]);
      f_norm=sqrt(f[i][0]*f[i][0]+f[i][1]*f[i][1]);
      *f_tot+=f_norm/(N1+N2);
      P+=f[i][0]*v[i][0]+f[i][1]*v[i][1];
      // 4. velocity modification
      for(int j=0;j<dim;j++){
        v[i][j]=(1.0-A)*v[i][j]+A*f[i][j]*v_norm/(f_norm+DBL_EPSILON);
      }
    }

    // converge criterion
    if(*f_tot<1.0e-12){
      break;
    }

    // 5. evaluate power P
    if(P>=0.0){
      // 6. P>0
      count++;
      if(count>5){
        count=0;
        dt=std::min(1.1*dt,dt_max);
        A*=0.99;
      }
    }else{
      // 7. P<0
      count=0;
      ini_array(v);
      dt*=0.5;
      A=0.1;
    }
  }
}

void output_coord(double (*x)[dim],double *a,double dphi,double gamma){
  char filename[128];
  std::ofstream file;
  if(dphi>0.0){
    sprintf(filename,"coord_AQSshear_ajam_dphi%1.1e_strain%1.4e.dat",fabs(dphi),gamma);
  }else{
    sprintf(filename,"coord_AQSshear_bjam_dphi%1.1e_strain%1.4e.dat",fabs(dphi),gamma);
  }
  file.open(filename);
  for(int i=0;i<N1+N2;i++){
    file <<x[i][0]<<"\t"<<x[i][1]<<"\t"<<a[i]<<std::endl;
  }
  file.close();
}

void output_phiJ(double phi){
  char filename[128];
  std::ofstream file;
  if(dphi_fixed>0.0){
    sprintf(filename,"phiJ_AQSshear_ajam_AQStrain_dphi%1.1e_N%d.dat",fabs(dphi_fixed),N1+N2);
  }else{
    sprintf(filename,"phiJ_AQSshear_bjam_AQStrain_dphi%1.1e_N%d.dat",fabs(dphi_fixed),N1+N2);
  }
  
  file.open(filename,std::ios::app);
  file << phi << std::endl;
  file.close();
}

void output_pressure(double p,double gamma){
  char filename[128];
  std::ofstream file;
  if(dphi_fixed>0.0){
    sprintf(filename,"pressure_AQSshear_ajam_AQStrain_dphi%1.1e_N%d.dat",fabs(dphi_fixed),N1+N2);
  }else{
    sprintf(filename,"pressure_AQSshear_bjam_AQStrain_dphi%1.1e_N%d.dat",fabs(dphi_fixed),N1+N2);
  }
  file.open(filename,std::ios::app);
  file << gamma << "\t" << p << std::endl;
  file.close();
}

void output_stress(double sxy0,double sxy,double gamma){
  char filename[128];
  std::ofstream file;
  if(dphi_fixed>0.0){
    sprintf(filename,"stress_AQSshear_ajam_AQStrain_dphi%1.1e_N%d.dat",fabs(dphi_fixed),N1+N2);
  }else{
    sprintf(filename,"stress_AQSshear_bjam_AQStrain_dphi%1.1e_N%d.dat",fabs(dphi_fixed),N1+N2);
  }
  file.open(filename,std::ios::app);
  file << gamma << "\t" << fabs(sxy-sxy0) << std::endl;
  file.close();
}

void output_stress_exam(double sxy,double gamma){
  char filename[128];
  std::ofstream file;
  if(dphi_fixed>0.0){
    sprintf(filename,"stress_AQStrain_ajam_exam_phi0.85_N%d.dat",N1+N2);
  }else{
    sprintf(filename,"stress_AQStrain_bjam_exam_phi0.85_N%d.dat",N1+N2);
  }
  file.open(filename,std::ios::app);
  file << gamma << "\t" << sxy << std::endl;
  file.close();
}

void initialize(int (*list)[Nl],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *sxy,double *f_tot,double *phi,double *disp_max,double L){
  *phi=phi0;
  ini_array(x);
  ini_coord(x,L);
  list_verlet(list,x,0.0,L);
  FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),0.0,L);
}

void make_jamming(int (*list)[Nl],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *sxy,double *f_tot,double *phi,double dphi,double *disp_max,double *L){
  while(*U<1.0e-16){
    affine_arithmetic(x,&(*L),&(*phi),dphi);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),0.0,*L);
  }
  while(*U>1.0e-12){
    affine_arithmetic(x,&(*L),&(*phi),-dphi);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),0.0,*L);
  }
  while(*U>1.0e-14){
    affine_arithmetic(x,&(*L),&(*phi),-dphi*1.0e-1);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),0.0,*L);
  }
  while(*U>1.0e-16){
    affine_arithmetic(x,&(*L),&(*phi),-dphi*1.0e-2);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),0.0,*L);
  }

  //output_coord(x,a,0.0,0.0);
  //output_phiJ(*phi);
}

void make_jamming_compress_training(int (*list)[Nl],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *sxy,double *f_tot,double *phi,double dphi,double *disp_max,double *L){
  while(*phi<=phi_max){
    affine_arithmetic(x,&(*L),&(*phi),dphi);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),0.0,*L);
  }
  while(*phi>=phi0){
    affine_arithmetic(x,&(*L),&(*phi),-dphi);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),0.0,*L);
  }
  while(*U<1.0e-16){
    affine_arithmetic(x,&(*L),&(*phi),dphi);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),0.0,*L);
  }
  while(*U>1.0e-12){
    affine_arithmetic(x,&(*L),&(*phi),-dphi);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),0.0,*L);
  }
  while(*U>1.0e-14){
    affine_arithmetic(x,&(*L),&(*phi),-dphi*1.0e-1);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),0.0,*L);
  }
  while(*U>1.0e-16){
    affine_arithmetic(x,&(*L),&(*phi),-dphi*1.0e-2);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&(*sxy),&(*f_tot),&(*disp_max),0.0,*L);
  }

  //output_coord(x,a,0.0,0.0);
  //output_phiJ(*phi);
}

void make_jamming_shear_training(int (*list)[Nl],double (*x)[dim],double (*x_update)[dim],double (*f)[dim],double *a,double *U,double *p,double *f_tot,double *phi,double dphi,double *disp_max,double *L){
  double sxy,sxy0;                  // local variable
  double gamma=0.0,dgamma=1.0e-3;   // local variable
  while(*phi<0.85){
    affine_arithmetic(x,&(*L),&(*phi),dphi);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*f_tot),&(*disp_max),0.0,*L);
    std::cout<<"phi="<<*phi<<std::endl;
  }
  std::cout<<"phi=0.85"<<"\t"<<"pressure="<<*p<<std::endl;
  while(gamma<=1.0){  // apply shear until strain=1.0
    gamma+=dgamma;
    apply_shear(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*disp_max),gamma,dgamma,*L);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*f_tot),&(*disp_max),gamma,*L);
    std::cout<<"gamma="<<gamma<<std::endl;
  }
  std::cout<<"gamma=1.0"<<"\t"<<"sxy0="<<sxy<<std::endl;
  dgamma=1.0e-5;
  while(sxy>1.0e-6){
    gamma-=dgamma;
    apply_shear(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*disp_max),gamma,-dgamma,*L);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*f_tot),&(*disp_max),gamma,*L);
    //output_stress_exam(sxy,gamma);
    std::cout<<std::fixed<<std::setprecision(8)<<"gamma="<<gamma<<"\t"<<"sxy="<<sxy<<std::endl;
  }
  dgamma=1.0e-7;
  while(sxy>1.0e-7){
    gamma-=dgamma;
    apply_shear(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*disp_max),gamma,-dgamma,*L);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*f_tot),&(*disp_max),gamma,*L);
    //output_stress_exam(sxy,gamma);
     std::cout<<std::fixed<<std::setprecision(8)<<"gamma="<<gamma<<"\t"<<"sxy="<<sxy<<std::endl;
  }
  dgamma=1.0e-9;
  while(sxy>0.0){
    gamma-=dgamma;
    apply_shear(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*disp_max),gamma,-dgamma,*L);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*f_tot),&(*disp_max),gamma,*L);
    //output_stress_exam(sxy,gamma);
    std::cout<<std::fixed<<std::setprecision(8)<<"gamma="<<gamma<<"\t"<<"sxy="<<sxy<<std::endl;
  }
  std::cout<<"sxy=0"<<std::endl;
  while(*U<1.0e-16){
    affine_arithmetic(x,&(*L),&(*phi),dphi);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*f_tot),&(*disp_max),0.0,*L);
  }
  while(*U>1.0e-12){
    affine_arithmetic(x,&(*L),&(*phi),-dphi);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*f_tot),&(*disp_max),0.0,*L);
  }
  while(*U>1.0e-14){
    affine_arithmetic(x,&(*L),&(*phi),-dphi*1.0e-1);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*f_tot),&(*disp_max),0.0,*L);
  }
  while(*U>1.0e-16){
    affine_arithmetic(x,&(*L),&(*phi),-dphi*1.0e-2);
    FIRE(list,x,x_update,f,a,&(*U),&(*p),&sxy,&(*f_tot),&(*disp_max),0.0,*L);
  }
  std::cout<<"phiJ="<<*phi<<std::endl;
  //output_coord(x,a,0.0,0.0);
  //output_phiJ(*phi);
}


int main(){
  double x[N1+N2][dim],x_update[N1+N2][dim],v[N1+N2][dim],f[N1+N2][dim],a[N1+N2];
  double U,p,sxy,sxy0,phi,phi_J,f_tot;
  double L=sqrt((N1*a1/2.0*a1/2.0+N2*a2/2.0*a2/2.0)*M_PI/phi0),disp_max=0.0,gamma=0.0,dgamma=1.0e-7;
  char filename[128];
  int list[N1+N2][Nl];
  int k=0,out=1;

  set_diameter(a);
  initialize(list,x,x_update,f,a,&U,&p,&sxy0,&f_tot,&phi,&disp_max,L);

  // calculate phi_J
  //make_jamming(list,x,x_update,f,a,&U,&p,&sxy0,&f_tot,&phi,dphi_a,&disp_max,&L);
  //make_jamming_compress_training(list,x,x_update,f,a,&U,&p,&sxy0,&f_tot,&phi,dphi_a,&disp_max,&L);
  make_jamming_shear_training(list,x,x_update,f,a,&U,&p,&f_tot,&phi,dphi_a,&disp_max,&L);
  phi_J=phi;
  /*
  // increase phi
  affine_arithmetic(x,&L,&phi,dphi_fixed);
  list_verlet(list,x,0.0,L);
  FIRE(list,x,x_update,f,a,&U,&p,&sxy0,&f_tot,&disp_max,0.0,L);
  //output_coord(x,a,dphi_fixed,0.0);
  
  // apply shear
  while(gamma<=1.0e-3){
    dgamma *= pow(10.0,0.1);  // logalithmic increasing (if linear plot is needed, let comment-out this line)
    gamma+=dgamma;
    apply_shear(list,x,x_update,f,a,&U,&p,&sxy,&disp_max,gamma,dgamma,L);
    FIRE(list,x,x_update,f,a,&U,&p,&sxy,&f_tot,&disp_max,gamma,L);
    //output_coord(x,a,dphi_fixed,gamma);
    output_pressure(p,gamma);
    output_stress(sxy0,sxy,gamma);
  }
  
  while(gamma<=1.0){
    dgamma = 1.0e-3;
    gamma+=dgamma;
    apply_shear(list,x,x_update,f,a,&U,&p,&sxy,&disp_max,gamma,dgamma,L);
    FIRE(list,x,x_update,f,a,&U,&p,&sxy,&f_tot,&disp_max,gamma,L);
    //output_coord(x,a,dphi_fixed,gamma);
    output_pressure(p,gamma);
    output_stress(sxy0,sxy,gamma);
  }
  */
  return 0;
}
